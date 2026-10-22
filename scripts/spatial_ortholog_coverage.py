#!/usr/bin/env python3
"""Offline MGI 1:1 ortholog mapping and real sparse spatial-gene coverage.

No title-case mapping, imputation, ROI construction or genotype significance
testing. The MGI report defines orthology; targeted Ensembl lookup supplies
stable mouse IDs so old mm10 symbols can match by gene identity.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd
from scipy import sparse

CANDIDATES = "FN1 ITGB1 PRSS23 CCN2 PLOD2 THBS1 LOX TGFBI GJA1 CD44 IGFBP3".split()
SOURCE_URL = "https://www.informatics.jax.org/downloads/reports/HOM_ProteinCoding.rpt"
MGI_COLUMNS = ["mgi_id", "mouse_symbol", "mouse_ncbi_gene_id", "hgnc_id", "human_symbol", "human_ncbi_gene_id"]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def no_version(value):
    return re.sub(r"\.\d+$", "", str(value).strip())


def prepare_mapping(programs, mgi_file, lookup_file=None):
    frozen = pd.read_csv(programs, sep="\t", keep_default_na=False)
    assert {"program", "rank", "gene"}.issubset(frozen.columns)
    assert not frozen.duplicated(["program", "gene"]).any()
    assert frozen.gene.ne("").all()
    raw = pd.read_csv(mgi_file, sep="\t", header=None, names=MGI_COLUMNS, dtype=str, keep_default_na=False)
    # MGI may list the same MGI/HGNC gene pair with multiple NCBI IDs. Preserve
    # those IDs together; never mistake alternate accessions for extra genes.
    report = raw.drop_duplicates().groupby(
        ["mgi_id", "mouse_symbol", "hgnc_id", "human_symbol"], sort=False, as_index=False
    ).agg(mouse_ncbi_gene_id=("mouse_ncbi_gene_id", lambda z: ";".join(sorted(set(z)))),
          human_ncbi_gene_id=("human_ncbi_gene_id", lambda z: ";".join(sorted(set(z)))))
    if report.human_symbol.duplicated().any() or report.mouse_symbol.duplicated().any():
        raise ValueError("MGI one-to-one report has nonidentical repeated human or mouse symbols")
    lookup = json.loads(Path(lookup_file).read_text(encoding="utf8")) if lookup_file else {}
    rows = []
    for gene in sorted(set(frozen.gene) | set(CANDIDATES)):
        found = report.loc[report.human_symbol == gene]
        row = {"human_gene": gene, "is_prespecified_candidate": gene in CANDIDATES,
               "frozen_programs": ";".join(frozen.loc[frozen.gene == gene, "program"].astype(str)),
               "mapping_status": "NOT_IN_STRICT_PROTEIN_CODING_ONE_TO_ONE_RESOURCE",
               "orthology_class": "UNMAPPED", "source_url": SOURCE_URL,
               "mouse_gene_symbol": "", "mouse_ensembl_gene_id": "",
               "mgi_id": "", "hgnc_id": "", "mouse_ncbi_gene_id": "", "human_ncbi_gene_id": "",
               "ensembl_current_biotype": "", "mouse_ID_lookup_status": "NOT_APPLICABLE"}
        if len(found) == 1:
            r = found.iloc[0]
            gene_lookup = lookup.get(r.mouse_symbol) or {}
            row.update(mapping_status="MGI_ONE_TO_ONE", orthology_class="one_to_one_protein_coding",
                       mouse_gene_symbol=r.mouse_symbol, mgi_id=r.mgi_id, hgnc_id=r.hgnc_id,
                       mouse_ncbi_gene_id=r.mouse_ncbi_gene_id, human_ncbi_gene_id=r.human_ncbi_gene_id,
                       mouse_ensembl_gene_id=gene_lookup.get("id", ""),
                       ensembl_current_biotype=gene_lookup.get("biotype", ""),
                       mouse_ID_lookup_status="ENSEMBL_SYMBOL_LOOKUP" if gene_lookup.get("id") else "NO_ENSEMBL_ID_IN_SNAPSHOT")
            if row["mouse_ensembl_gene_id"] and not row["mouse_ensembl_gene_id"].startswith("ENSMUSG"):
                raise ValueError(f"Unexpected mouse stable ID for {gene}: {row['mouse_ensembl_gene_id']}")
        rows.append(row)
    return frozen, pd.DataFrame(rows), len(raw), len(report)


def capture_coverage(folder, mapping):
    folder = Path(folder)
    var = pd.read_csv(folder / "MOUSE_SYMBOL_FEATURES.tsv", sep="\t", keep_default_na=False)
    obs = pd.read_csv(folder / "BIN48_METADATA_QC.tsv.gz", sep="\t", keep_default_na=False)
    x = sparse.load_npz(folder / "counts48_symbols_qc.npz").tocsr()
    assert x.shape == (len(obs), len(var)) and sparse.issparse(x)
    assert var.mouse_gene_symbol.is_unique and obs.spot_id.is_unique
    ensembl_indices = {}
    for i, text in enumerate(var.mouse_ensembl_ids):
        for eid in set(no_version(z) for z in str(text).split(";") if z):
            ensembl_indices.setdefault(eid, []).append(i)
    symbol_indices = dict(zip(var.mouse_gene_symbol, range(len(var))))
    library = np.asarray(x.sum(axis=1)).ravel().astype(float)
    rows = []
    for r in mapping.to_dict("records"):
        result = dict(r)
        result.update(capture_id=str(obs.capture_id.iloc[0]) if len(obs) else folder.name,
                      source_GSM=str(obs.source_GSM.iloc[0]) if len(obs) else "",
                      genotype=str(obs.genotype.iloc[0]) if len(obs) else "",
                      n_qc_bins=len(obs), native_mouse_feature_present=False,
                      matched_native_mouse_symbol="", matched_native_ensembl_ids="",
                      feature_match_method="NOT_MAPPED", n_detected_bins=0,
                      detected_bin_fraction=0.0, total_raw_UMI=0,
                      mean_log1p_per_10000=0.0, assessability="NO_ONE_TO_ONE_MAPPING")
        if r["mapping_status"] != "MGI_ONE_TO_ONE":
            rows.append(result)
            continue
        indices = ensembl_indices.get(no_version(r["mouse_ensembl_gene_id"]), []) if r["mouse_ensembl_gene_id"] else []
        method = "STABLE_MOUSE_ENSEMBL_ID"
        if not indices and r["mouse_gene_symbol"] in symbol_indices:
            indices = [symbol_indices[r["mouse_gene_symbol"]]]
            method = "EXACT_MGI_MOUSE_SYMBOL"
        if len(indices) > 1:
            result.update(assessability="AMBIGUOUS_NATIVE_FEATURE_MATCH", feature_match_method="MULTIPLE_FEATURES_SAME_ENSEMBL_ID")
        elif len(indices) == 0:
            result.update(assessability="ORTHOLOG_NOT_IN_RELEASED_MATRIX_FEATURES", feature_match_method="NO_EXACT_STABLE_ID_OR_SYMBOL_MATCH")
        else:
            index = indices[0]
            column = x[:, index].tocoo()
            normalized_nonzero = np.log1p(column.data.astype(float) / np.maximum(library[column.row], 1) * 10000)
            n_detected = int(np.count_nonzero(column.data))
            result.update(native_mouse_feature_present=True,
                          matched_native_mouse_symbol=str(var.iloc[index].mouse_gene_symbol),
                          matched_native_ensembl_ids=str(var.iloc[index].mouse_ensembl_ids),
                          feature_match_method=method, n_detected_bins=n_detected,
                          detected_bin_fraction=n_detected / len(obs) if len(obs) else 0.0,
                          total_raw_UMI=int(column.data.sum()),
                          mean_log1p_per_10000=float(normalized_nonzero.sum() / len(obs)) if len(obs) else 0.0,
                          assessability="DETECTED_IN_QC_BINS" if n_detected else "FEATURE_PRESENT_BUT_ZERO_IN_QC_BINS")
        rows.append(result)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--programs", required=True, type=Path)
    source = ap.add_mutually_exclusive_group(required=True)
    source.add_argument("--mgi-report", type=Path, help="Explicit alternate current source; may change mappings")
    source.add_argument("--targeted-mapping", type=Path,
                        help="Frozen, source-attributed limited MGI/Ensembl mapping used in the completed analysis")
    ap.add_argument("--expected-mapping-sha256")
    ap.add_argument("--ensembl-lookup", type=Path)
    ap.add_argument("--source-metadata", type=Path)
    ap.add_argument("--prepared48", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()
    out = a.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if a.targeted_mapping:
        if a.expected_mapping_sha256 and sha(a.targeted_mapping) != a.expected_mapping_sha256:
            raise ValueError("Frozen targeted mapping hash mismatch")
        frozen = pd.read_csv(a.programs, sep="\t", keep_default_na=False)
        mapping = pd.read_csv(a.targeted_mapping, sep="\t", keep_default_na=False)
        if mapping.human_gene.duplicated().any() or not (set(frozen.gene) | set(CANDIDATES)).issubset(set(mapping.human_gene)):
            raise ValueError("Targeted mapping is ambiguous or does not cover the frozen program and candidate query")
        one = mapping.mapping_status.eq("MGI_ONE_TO_ONE")
        if not mapping.loc[one, "orthology_class"].eq("one_to_one_protein_coding").all():
            raise ValueError("Inconsistent strict one-to-one classification")
        if mapping.loc[one, "mouse_gene_symbol"].eq("").any():
            raise ValueError("Mapped human query lacks mouse ortholog symbol")
        nraw = nuniq = None
    else:
        frozen, mapping, nraw, nuniq = prepare_mapping(a.programs, a.mgi_report, a.ensembl_lookup)
    mapping.to_csv(out / "HUMAN_MOUSE_ONE_TO_ONE_TARGETED_MAPPING.tsv", sep="\t", index=False)
    frozen.merge(mapping, how="left", left_on="gene", right_on="human_gene", validate="many_to_one").to_csv(
        out / "FROZEN_K5_WITH_MOUSE_ORTHOLOGS.tsv", sep="\t", index=False)
    records = []
    for program, f in frozen.groupby("program", sort=True):
        q = mapping.set_index("human_gene").loc[f.gene]
        records.append({"program": program, "frozen_gene_count": len(f),
                        "strict_one_to_one_mapped": int(q.mapping_status.eq("MGI_ONE_TO_ONE").sum()),
                        "strict_one_to_one_fraction": float(q.mapping_status.eq("MGI_ONE_TO_ONE").mean()),
                        "unmapped_genes": ";".join(q.index[q.mapping_status.ne("MGI_ONE_TO_ONE")])})
    pd.DataFrame(records).to_csv(out / "PROGRAM_ORTHOLOGY_COVERAGE_ONLY.tsv", sep="\t", index=False)
    status = {"status": "MAPPING_PREPARED_NO_REAL_MATRIX_READ", "finished_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
              "mgi_source": SOURCE_URL, "mgi_report_sha256": sha(a.mgi_report) if a.mgi_report else None,
              "targeted_mapping_sha256": sha(a.targeted_mapping) if a.targeted_mapping else None,
              "mapping_mode": "frozen_targeted_snapshot" if a.targeted_mapping else "explicit_current_MGI_report",
              "mgi_raw_rows": nraw, "mgi_unique_gene_pairs": nuniq,
              "alternate_accession_report_rows_consolidated": nraw - nuniq if nraw is not None else None,
              "frozen_programs_sha256": sha(a.programs), "query_genes": len(mapping),
              "query_one_to_one_mapped": int(mapping.mapping_status.eq("MGI_ONE_TO_ONE").sum()),
              "titlecase_or_inferred_orthology": False, "ROIs_read_or_changed": False,
              "source_metadata": json.loads(a.source_metadata.read_text()) if a.source_metadata else None}
    if a.prepared48:
        folders = sorted(p.parent for p in a.prepared48.glob("*/counts48_symbols_qc.npz"))
        if not folders:
            raise ValueError("No actual prepared 48-um capture sparse matrices found")
        coverage = pd.concat([capture_coverage(p, mapping) for p in folders], ignore_index=True)
        coverage.to_csv(out / "HUMAN_MOUSE_ORTHOLOG_COVERAGE.tsv", sep="\t", index=False)
        detailed = frozen[["program", "gene", "rank"]].merge(coverage, left_on="gene", right_on="human_gene", how="left", validate="one_to_many")
        detailed["detected_mouse_symbol"] = detailed.matched_native_mouse_symbol.where(detailed.n_detected_bins > 0, "")
        detailed.to_csv(out / "PROGRAM_GENE_SPATIAL_COVERAGE.tsv", sep="\t", index=False)
        summary = detailed.groupby(["capture_id", "source_GSM", "genotype", "program"], sort=True).agg(
            frozen_gene_count=("gene", "size"), n_one_to_one_mapped=("mapping_status", lambda z: int(z.eq("MGI_ONE_TO_ONE").sum())),
            n_native_features_present=("native_mouse_feature_present", "sum"), n_genes_detected=("n_detected_bins", lambda z: int((z > 0).sum())),
            detected_gene_list=("detected_mouse_symbol", lambda z: ";".join(v for v in z if v)),
        ).reset_index()
        summary.to_csv(out / "PROGRAM_ORTHOLOG_COVERAGE_SUMMARY.tsv", sep="\t", index=False)
        status.update(status="DONE_MAPPING_AND_REAL_MATRIX_COVERAGE", actual_capture_count=len(folders))
    (out / "ORTHOLOG_MAPPING_STATUS.json").write_text(json.dumps(status, indent=2) + "\n", encoding="utf8")
    print(json.dumps(status, indent=2), flush=True)


if __name__ == "__main__":
    main()
