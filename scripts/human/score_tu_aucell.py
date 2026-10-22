#!/usr/bin/env python3
"""Actual AUCell via pinned decoupler 2.2.0, with explicit effective coverage.

Uses the official sparse AnnData interface. Frozen signatures remain unchanged;
the official method prunes unavailable/empty targets and normalises each AUC by
its effective target set, which is recorded rather than called a 25-gene UCell.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
from importlib import import_module
import inspect
import json
import os
from pathlib import Path
import platform
import sys

PROGRAMS = [f"GEP{i}" for i in range(1, 6)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--bsize", type=int, default=256)
    parser.add_argument("--n-up", type=int, default=1500)
    parser.add_argument("--tmin", type=int, default=5)
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    if not 1 <= args.bsize <= 256:
        raise ValueError("Use sparse AUCell bsize <= 256")
    # Configure the optional project add-on library without modifying shared envs.
    addons = args.work / "envs" / "python_addons"
    if addons.exists():
        sys.path.insert(0, str(addons))
    os.environ.setdefault("NUMBA_NUM_THREADS", str(args.threads))
    os.environ.setdefault("OMP_NUM_THREADS", str(args.threads))
    import numpy as np
    import pandas as pd
    from scipy import sparse
    import decoupler as dc
    import anndata as ad
    import numba
    version = importlib.metadata.version("decoupler")
    if version != "2.2.0":
        raise RuntimeError(f"Expected pinned decoupler 2.2.0; found {version}")
    numba.set_num_threads(args.threads)
    np.random.seed(20261021)
    out = args.work / "01_human_programs"
    mapping = pd.read_csv(out / "PROGRAM_GENE_MAPPING.tsv", sep="\t").fillna({"mapped_gene": ""})
    frozen = pd.read_csv(out / "FROZEN_K5_PROGRAMS.tsv", sep="\t")
    if not frozen.groupby("program").size().reindex(PROGRAMS).eq(25).all():
        raise ValueError("Exactly five original frozen Top25 programs required")
    all_scores, effective_rows, sample_records, coverage_rows = [], [], [], []
    # Bind the run to the actual installed implementation and its documented tie rule.
    au_module = import_module("decoupler.mt._aucell")
    data_module = import_module("decoupler.pp.data")
    implementation = {}
    for name, module in [("aucell", au_module), ("data_extract", data_module)]:
        source = Path(inspect.getfile(module))
        implementation[name] = {"path": str(source), "sha256": hashlib.sha256(source.read_bytes()).hexdigest()}
    for sample in range(1, 9):
        gsm = "GSM" + str(5048707 + sample)
        sample_dir = args.work / "data_cache" / "Tu_sparse" / gsm
        mat = sparse.load_npz(sample_dir / "expression_gene_by_cell.npz").T.tocsr()
        genes = pd.read_csv(sample_dir / "genes.tsv", sep="\t").gene.astype(str).to_numpy()
        cells = pd.read_csv(sample_dir / "cells.tsv", sep="\t").cell_id.astype(str).to_numpy()
        if mat.shape != (len(cells), len(genes)) or len(set(genes)) != len(genes):
            raise ValueError(f"Sparse dimensions/feature identities invalid: {gsm}")
        if not np.isfinite(mat.data).all() or (mat.data < 0).any():
            raise ValueError(f"Invalid expression: {gsm}")
        nonempty_features = np.asarray(mat.getnnz(axis=0) > 0).ravel()
        effective_features = set(genes[nonempty_features])
        input_features = set(genes)
        if len(effective_features) < args.n_up:
            raise ValueError(f"{gsm}: fewer nonempty genes than the fixed n_up={args.n_up}; do not silently tune")
        selected = mapping[mapping.gsm.eq(gsm)].copy()
        net_rows = []
        for row in selected.itertuples():
            mapped = str(row.mapped_gene)
            included = bool(mapped and mapped in effective_features)
            reason = "shared_nonempty_target" if included else "not_mapped" if not mapped else "zero_feature_removed_by_official_empty_true"
            effective_rows.append({"gsm": gsm, "program": row.program,
                "original_frozen_gene": row.original_gene, "mapped_gene": mapped,
                "present_in_input_features": bool(mapped and mapped in input_features),
                "present_nonempty_feature": included, "effective_target_status": reason,
                "score_method": "AUCell_decoupler_2.2.0"})
            if mapped:
                net_rows.append({"source": row.program, "target": mapped})
        net = pd.DataFrame(net_rows).drop_duplicates(["source", "target"])
        effective_counts = net[net.target.isin(effective_features)].groupby("source").size()
        available_programs = effective_counts[effective_counts >= args.tmin].index.tolist()
        for program in PROGRAMS:
            count = int(effective_counts.get(program, 0))
            coverage_rows.append({"gsm": gsm, "program": program, "n_original_frozen_genes": 25,
                "n_AUCell_shared_nonempty_targets": count, "program_retained_after_tmin": count >= args.tmin,
                "score_method": "AUCell_decoupler_2.2.0", "tmin": args.tmin,
                "algorithm_denominator": "effective shared nonempty target set; official maximum AUC"})
        # If one program is technically not assessable, retain NA for it; never
        # fabricate zero enrichment. Other programs can still complete.
        if not available_programs:
            raise RuntimeError(f"No program has tmin={args.tmin} shared nonempty targets: {gsm}")
        net.to_csv(sample_dir / "AUCell_input_network.tsv", sep="\t", index=False)
        data = ad.AnnData(X=mat, obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=genes))
        print(json.dumps({"start_AUCell": gsm, "shape": mat.shape,
            "nonempty_genes": len(effective_features), "effective_program_targets": effective_counts.to_dict(),
            "method": "decoupler.mt.aucell", "n_up": args.n_up, "tmin": args.tmin}), flush=True)
        returned = dc.mt.aucell(data=data, net=net, tmin=args.tmin, empty=True,
                               bsize=args.bsize, n_up=args.n_up, verbose=True)
        scored = returned if isinstance(returned, ad.AnnData) else data
        if "score_aucell" not in scored.obsm:
            raise RuntimeError("Official decoupler did not produce score_aucell")
        scores = scored.obsm["score_aucell"].copy()
        if not isinstance(scores, pd.DataFrame):
            raise TypeError("Unexpected decoupler score output type")
        if not np.isfinite(scores.to_numpy()).all():
            raise ValueError(f"Nonfinite AUCell output: {gsm}")
        if not set(scores.index).issubset(set(cells)) or scores.index.duplicated().any():
            raise ValueError(f"Output cell identity mismatch: {gsm}")
        for program in PROGRAMS:
            if program not in scores:
                scores[program] = np.nan
        scores = scores[PROGRAMS]
        scores.columns = [p + "_top25_AUCell" for p in PROGRAMS]
        scores.insert(0, "cell_id", scores.index)
        scores["score_method"] = "AUCell_decoupler_2.2.0"
        scores.to_csv(sample_dir / "AUCell_scores.tsv", sep="\t", index=False)
        all_scores.append(scores)
        sample_records.append({"gsm": gsm, "n_input_cells": len(cells), "n_scored_nonempty_cells": len(scores),
            "nonempty_features": len(effective_features),
            "effective_target_counts": {p: int(effective_counts.get(p, 0)) for p in PROGRAMS},
            "programs_pruned_by_tmin": [p for p in PROGRAMS if p not in available_programs]})
        del data, scored, mat
        print(json.dumps({"complete_AUCell": gsm, "scored_cells": len(scores)}), flush=True)
    combined = pd.concat(all_scores, ignore_index=True)
    combined.to_csv(out / "TU_PROGRAM_AUCELL_CELL_SCORES.tsv.gz", sep="\t", index=False)
    pd.DataFrame(effective_rows).to_csv(out / "AUCELL_EFFECTIVE_PROGRAM_TARGETS.tsv", sep="\t", index=False)
    pd.DataFrame(coverage_rows).to_csv(out / "AUCELL_EFFECTIVE_PROGRAM_COVERAGE.tsv", sep="\t", index=False)
    (out / "AUCELL_RUN_METADATA.json").write_text(json.dumps({
        "status": "DONE" if all(not r["programs_pruned_by_tmin"] for r in sample_records) else "PARTIAL_PROGRAM_COVERAGE",
        "algorithm": "AUCell", "implementation": "official decoupler.mt.aucell", "version": version,
        "data_interface": "official sparse AnnData; pinned 2.2.0 tuple/list sparse extraction is unsuitable",
        "python": platform.python_version(), "n_up": args.n_up, "tmin": args.tmin, "bsize": args.bsize,
        "threads": args.threads, "external_numpy_seed": 20261021,
        "ties": "Official pp.data._break_ties uses default_rng(seed=0) to permute feature order; scipy.stats.rankdata(method='ordinal') then resolves ties. No external tie seed argument exists.",
        "denominator": "Official AUCell maximum AUC uses each shared nonempty target set after feature/network pruning; not a forced 25-member denominator.",
        "normalization": "Published processed expression used directly; no rank-changing feature scaling or extra normalization.",
        "frozen_signature_boundary": "Original five Top25 lists retained unchanged; mapping/coverage/effective target pruning tabulated separately.",
        "implementation_files": implementation, "samples": sample_records,
        "p_values": "AUCell provides enrichment scores, not inferential P values.",
        "method_choice": "Fixed Tu AUCell implementation; keep its numerical scale distinct from discovery UCell. No R UCell run is claimed for Tu.",
    }, indent=2), encoding="utf-8")
    (out / "AUCELL_METHOD_ACTUALLY_RUN.txt").write_text(
        "Official decoupler 2.2.0 AUCell on sparse AnnData; n_up=1500, tmin=5, bsize<=256.\n"
        "This is AUCell, not UCell. Historical Swahn UCell values are retained with their original method label.\n"
        "Official internal tie permutation seed=0; ordinal ranks. Effective target denominators follow official pruning.\n"
        "Scores from different methods are not numerically comparable. No group P values are produced.\n", encoding="utf-8")
    print(json.dumps({"AUCell_cells_total": len(combined), "version": version,
                      "output": str(out / "TU_PROGRAM_AUCELL_CELL_SCORES.tsv.gz")}), flush=True)


if __name__ == "__main__":
    main()
