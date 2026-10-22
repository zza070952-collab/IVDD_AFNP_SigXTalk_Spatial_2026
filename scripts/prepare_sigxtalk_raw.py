#!/usr/bin/env python3
"""Prepare donor-resolved SigXTalk input from genuine public GEO integer counts.

No expression matrix is inferred from program scores. Raw mode runs the sibling
human discovery interface plus the fixed AF identity pass; cached mode reads real
cell-level QC/identity/UCell tables and independently resamples their barcodes.
"""
from __future__ import annotations

import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
import argparse
import datetime as dt
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
import urllib.request

import numpy as np
import pandas as pd
from scipy import io, sparse

REPO = Path(__file__).resolve().parents[1]
SEED = 20261017
MIN_TAIL = 20
CAP = 150
ALIASES_FOR_HISTORICAL_NP_QC = {
    "CTGF": "CCN2", "CYR61": "CCN1", "WISP1": "CCN4", "WISP2": "CCN5", "WISP3": "CCN6"}
SUFFIXES = ["matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"]


class EligibilityError(ValueError):
    def __init__(self, message, eligibility):
        super().__init__(message)
        self.eligibility = eligibility


def sha(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def dump(value, path):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def truth(values):
    text = values.astype(str).str.lower()
    bad = ~text.isin(["true", "false", "1", "0"])
    if bad.any():
        raise ValueError("Unrecognised boolean metadata values: " + str(text[bad].unique().tolist()))
    return text.isin(["true", "1"])


def data_path(root, value):
    path = Path(value)
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def resource_path(value):
    path = Path(value)
    return path.resolve() if path.is_absolute() else (REPO / path).resolve()


def writable(root, path):
    if not path.resolve().is_relative_to(root.resolve()):
        raise ValueError("Writable output must remain in the new project root: " + str(path))
    return path


def read_cells(path, columns):
    frame = pd.read_csv(path, sep="\t", dtype={"cell_id": str, "gsm": str, "donor_id": str})
    if not set(columns).issubset(frame):
        raise ValueError(f"{path} missing actual cell-level fields: {sorted(set(columns)-set(frame))}")
    if frame.cell_id.isna().any() or frame.cell_id.duplicated().any():
        raise ValueError("Missing or duplicate cell ID in " + str(path))
    return frame.set_index("cell_id", drop=False)


def validate_manifest(path):
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {"series", "gsm", "donor_id", "tissue", "age", "sex", "disc", "condition", "matrix_url"}
    if not required.issubset(frame) or frame.gsm.duplicated().any():
        raise ValueError("Invalid public sample manifest")
    if len(frame) != 16 or frame.donor_id.nunique() != 8 or set(frame.series) != {"GSE230808"}:
        raise ValueError("The declared raw-selection input requires eight GSE230808 AF/NP pairs")
    if set(frame.condition) != {"diseased"}:
        raise ValueError("Do not mix healthy/external donors into this disease-internal selection")
    for donor, group in frame.groupby("donor_id", sort=True):
        if len(group) != 2 or set(group.tissue) != {"AF", "NP"}:
            raise ValueError("Unpaired donor: " + donor)
        for field in ["age", "sex", "disc"]:
            if group[field].nunique() != 1 or not group[field].iloc[0]:
                raise ValueError("Unverified donor pairing field: " + donor + " " + field)
    if not frame.matrix_url.str.startswith("https://ftp.ncbi.nlm.nih.gov/geo/samples/").all():
        raise ValueError("The distributed raw manifest must use official GEO URLs")
    return frame


def validate_cell_sources(frame, manifest, tissue):
    samples = manifest.set_index("gsm")
    for gsm, group in frame.groupby("gsm", sort=False):
        if gsm not in samples.index or samples.loc[gsm, "tissue"] != tissue:
            raise ValueError("Cell table has wrong tissue/accession: " + gsm)
        if set(group.donor_id) != {samples.loc[gsm, "donor_id"]}:
            raise ValueError("Cell donor does not match the verified public sample: " + gsm)
        if not group.cell_id.str.startswith(gsm + ":").all():
            raise ValueError("Expected an actual GSM:barcode cell identifier")
    if frame.donor_id.nunique() != 8:
        raise ValueError("Both actual metadata inputs must retain all eight donors")


def select_cells(np_metadata, scores, af_metadata, minimum=MIN_TAIL, cap=CAP, seed=SEED,
                 expected_eligible=6):
    """The predeclared selection, independent of any historical selected-cell list."""
    metadata = np_metadata.copy()
    if not metadata.index.isin(scores.index).all():
        raise ValueError("Discovery cell-level UCell scores are incomplete")
    metadata["GEP2"] = pd.to_numeric(scores.loc[metadata.index, "GEP2_top25_UCell"], errors="raise").to_numpy()
    if not np.isfinite(metadata.GEP2.to_numpy()).all():
        raise ValueError("Non-finite GEP2 UCell input; no missing-value replacement is allowed")
    high_confidence = truth(metadata.high_confidence_NP)
    if "pass_qc" in metadata and (high_confidence & ~truth(metadata.pass_qc)).any():
        raise ValueError("High-confidence NP input contains a failed-QC cell")
    broad_af = truth(af_metadata.broad_stromal)
    if "pass_qc" in af_metadata and (broad_af & ~truth(af_metadata.pass_qc)).any():
        raise ValueError("Broad-stromal AF input contains a failed-QC cell")
    if "predicted_doublet" in af_metadata and (broad_af & truth(af_metadata.predicted_doublet)).any():
        raise ValueError("Broad-stromal AF input contains a predicted doublet")
    rng = np.random.default_rng(seed)
    parts, eligibility = [], []
    for donor, group in metadata.groupby("donor_id", sort=True, observed=True):
        subset = group[high_confidence.loc[group.index]]
        quantiles = subset.GEP2.quantile([.25, .75])
        q25, q75 = float(quantiles.iloc[0]), float(quantiles.iloc[1])
        low, high = subset[subset.GEP2 <= q25], subset[subset.GEP2 >= q75]
        eligible = len(low) >= minimum and len(high) >= minimum and q75 > q25
        eligibility.append({"donor_id": donor, "n_NP": len(subset), "n_low": len(low), "n_high": len(high),
            "q25": q25, "q75": q75, "eligible": bool(eligible),
            "eligibility_reason": "eligible" if eligible else "insufficient_tail_support_or_nonseparating_quartiles"})
        if eligible:
            for name, candidates in [("NP_GEP2_low", low), ("NP_GEP2_high", high)]:
                selected = sorted(rng.choice(sorted(candidates.index), min(cap, len(candidates)), replace=False))
                part = candidates.loc[selected].copy()
                part["bridge_group"], part["tissue"] = name, "NP"
                parts.append(part)
    eligible_donors = [row["donor_id"] for row in eligibility if row["eligible"]]
    if len(eligible_donors) != expected_eligible:
        raise EligibilityError(f"Actual fixed rules yielded {len(eligible_donors)} eligible pairs; expected {expected_eligible}. Preserve the changed result and investigate, do not change thresholds.", pd.DataFrame(eligibility))
    for donor in eligible_donors:
        candidates = af_metadata[(af_metadata.donor_id == donor) & broad_af]
        if len(candidates) < minimum:
            raise ValueError("Insufficient matched broad-stromal AF cells: " + donor)
        selected = sorted(rng.choice(sorted(candidates.index), min(cap, len(candidates)), replace=False))
        part = candidates.loc[selected].copy()
        part["bridge_group"], part["tissue"] = "AF_stromal", "AF"
        parts.append(part)
    selected = pd.concat(parts, axis=0)
    if not selected.index.is_unique:
        raise ValueError("A cell entered more than one selected group")
    selected["GEP2_top25_UCell"] = selected.get("GEP2")
    return selected, pd.DataFrame(eligibility)


def find_raw(raw_dir, gsm, suffix, required=True):
    candidates = list(raw_dir.glob(gsm + "*_" + suffix)) + list((raw_dir / gsm).glob(gsm + "*_" + suffix))
    if len(candidates) > 1 or (required and len(candidates) != 1):
        raise ValueError(f"Exactly one raw file required: {gsm} {suffix}: {candidates}")
    return candidates[0] if candidates else None


def fetch(url, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + ".part")
    if partial.exists():
        raise RuntimeError("Inspect the partial download before another bounded attempt: " + str(partial))
    start = time.time()
    request = urllib.request.Request(url, headers={"User-Agent": "IVDD-public-raw-preparer/1"})
    with urllib.request.urlopen(request, timeout=90) as response, partial.open("wb") as handle:
        expected = response.headers.get("Content-Length")
        for block in iter(lambda: response.read(1024 * 1024), b""):
            handle.write(block)
    if expected is not None and partial.stat().st_size != int(expected):
        raise ValueError("Incomplete public transfer: " + url)
    partial.replace(destination)
    return round(time.time() - start, 3)


def obtain_raw(manifest, raw_dir, root, allow_download):
    receipts = []
    for row in manifest.to_dict("records"):
        for suffix in SUFFIXES:
            url = row["matrix_url"].replace("matrix.mtx.gz", suffix)
            existing = find_raw(raw_dir, row["gsm"], suffix, required=False)
            if existing is None:
                if not allow_download:
                    raise FileNotFoundError("Actual public raw count file unavailable: " + url)
                writable(root, raw_dir)
                path = raw_dir / url.rsplit("/", 1)[-1]
                elapsed, status = fetch(url, path), "DOWNLOADED_THIS_RUN"
            else:
                path, elapsed, status = existing, None, "REUSED_ACTUAL_RAW_FILE"
            receipts.append({"gsm": row["gsm"], "donor_id": row["donor_id"], "tissue": row["tissue"],
                "url": url, "path": str(path), "bytes": path.stat().st_size, "sha256": sha(path),
                "status": status, "download_seconds": elapsed})
    return receipts


def unique_symbols(symbols):
    """Deterministic make-unique names, retaining the first submitted occurrence."""
    values = list(map(str, symbols))
    occupied, seen, counter = set(values), set(), {}
    for index, value in enumerate(values):
        if value not in seen:
            seen.add(value)
            continue
        number = counter.get(value, 0) + 1
        while value + "-" + str(number) in occupied:
            number += 1
        counter[value] = number
        replacement = value + "-" + str(number)
        occupied.add(replacement)
        values[index] = replacement
    return pd.Index(values)


def features(path):
    frame = pd.read_csv(path, sep="\t", header=None, dtype=str, keep_default_na=False)
    if frame.shape[1] < 2:
        raise ValueError("Expected a real 10x feature ID and symbol file")
    keep = frame[2].eq("Gene Expression").to_numpy() if frame.shape[1] >= 3 else np.ones(len(frame), bool)
    submitted = frame.loc[keep, 1].to_numpy()
    if (submitted == "").any():
        raise ValueError("Empty published gene symbol; an explicit ID mapping policy is required")
    return frame, keep, submitted, unique_symbols(submitted)


def reconstruct_counts(meta, manifest, raw_dir, output_dir, audit_dir, allow_intersection=False,
                       require_cached_umi=True):
    sample_index = manifest.set_index("gsm")
    per_sample, gene_rows = {}, []
    common, first = None, None
    for gsm in sorted(meta.gsm.unique()):
        frame, keep, submitted, names = features(find_raw(raw_dir, gsm, "features.tsv.gz"))
        per_sample[gsm] = (frame, keep, submitted, names)
        common = set(names) if common is None else common.intersection(names)
        if first is None:
            first = names
        for raw_index, feature_id, symbol, made_unique in zip(np.flatnonzero(keep), frame.loc[keep, 0], submitted, names):
            gene_rows.append({"gsm": gsm, "raw_feature_row_zero_based": int(raw_index), "feature_id": feature_id,
                "submitted_symbol": symbol, "output_unique_symbol": made_unique})
    if not common:
        raise ValueError("No shared Gene Expression feature universe")
    genes = pd.Index([name for name in first if name in common])
    gene_table = pd.DataFrame(gene_rows)
    gene_table["in_common_output"] = gene_table.output_unique_symbol.isin(common)
    inconsistent = gene_table[gene_table.in_common_output].groupby("output_unique_symbol").feature_id.nunique()
    if (inconsistent > 1).any():
        raise ValueError("A shared gene symbol maps to different submitted feature IDs across raw samples")
    gene_table.to_csv(audit_dir / "RAW_GENE_FEATURE_ALIGNMENT.tsv.gz", sep="\t", index=False)
    if not allow_intersection and any(set(item[3]) != common for item in per_sample.values()):
        raise ValueError("Raw GEX feature universes differ. Alignment audit preserved; authorise and disclose a common-gene intersection explicitly if intended.")
    matrices, column_ids, audits, library_rows = [], [], [], []
    for gsm in sorted(per_sample):
        frame, keep, submitted, names = per_sample[gsm]
        matrix_path = find_raw(raw_dir, gsm, "matrix.mtx.gz")
        barcodes = pd.read_csv(find_raw(raw_dir, gsm, "barcodes.tsv.gz"), sep="\t", header=None, dtype=str)[0]
        if barcodes.duplicated().any():
            raise ValueError("Duplicate raw barcodes: " + gsm)
        with gzip.open(matrix_path, "rb") as handle:
            counts = io.mmread(handle).tocsc()
        if counts.shape != (len(frame), len(barcodes)):
            raise ValueError("10x feature/barcode/matrix dimensions disagree: " + gsm)
        if (not np.isfinite(counts.data).all() or (counts.data < 0).any()
                or not np.equal(counts.data, np.floor(counts.data)).all()
                or (counts.data > np.iinfo(np.int64).max).any()):
            raise ValueError("The count input is not finite nonnegative integer raw data: " + gsm)
        counts = counts[keep].tocsc()
        selected = meta[meta.gsm == gsm]
        if set(selected.donor_id) != {sample_index.loc[gsm, "donor_id"]} or set(selected.tissue) != {sample_index.loc[gsm, "tissue"]}:
            raise ValueError("Raw sample identity disagrees with selected metadata")
        raw_ids = pd.Index(gsm + ":" + barcodes)
        positions = raw_ids.get_indexer(selected.index)
        if (positions < 0).any():
            raise ValueError("Selected actual barcode is absent from the public raw matrix: " + gsm)
        selected_all = counts[:, positions].astype(np.int64)
        selected_all.eliminate_zeros()
        raw_umi = np.asarray(selected_all.sum(axis=0)).ravel()
        raw_features = np.asarray((selected_all > 0).sum(axis=0)).ravel()
        canonical = pd.Index([ALIASES_FOR_HISTORICAL_NP_QC.get(name, name) for name in submitted])
        canonical_unique = pd.Index(pd.unique(canonical))
        collapse = sparse.csr_matrix((np.ones(len(canonical), dtype=np.int64),
            (canonical_unique.get_indexer(canonical), np.arange(len(canonical)))),
            shape=(len(canonical_unique), len(canonical)))
        canonical_features = np.asarray(((collapse @ selected_all) > 0).sum(axis=0)).ravel()
        result = selected_all[names.get_indexer(genes), :].tocsc()
        final_umi = np.asarray(result.sum(axis=0)).ravel()
        for index, cell in enumerate(selected.index):
            row = selected.loc[cell]
            cached_umi, cached_features = float(row.n_counts), int(row.n_features)
            historical_features = int(canonical_features[index]) if row.tissue == "NP" else int(raw_features[index])
            audits.append({"cell_id": cell, "gsm": gsm, "donor_id": row.donor_id, "tissue": row.tissue,
                "bridge_group": row.bridge_group, "barcode": cell.split(":", 1)[1], "raw_barcode_found": True,
                "raw_total_UMI": int(raw_umi[index]), "cached_total_UMI": cached_umi,
                "total_UMI_equal": bool(raw_umi[index] == cached_umi),
                "output_common_gene_UMI": int(final_umi[index]), "common_gene_UMI_loss": int(raw_umi[index] - final_umi[index]),
                "raw_detected_genes": int(raw_features[index]), "cached_detected_genes": cached_features,
                "detected_genes_equal": bool(raw_features[index] == cached_features),
                "historical_equivalent_detected_genes": historical_features,
                "historical_gene_count_equal": historical_features == cached_features,
                "historical_feature_policy": "canonical_symbol_duplicate_sum" if row.tissue == "NP" else "submitted_symbols_made_unique",
                "source_count_file": str(matrix_path),
                "GEP2_top25_UCell": float(row.GEP2) if pd.notna(row.get("GEP2")) else np.nan})
        matrices.append(result)
        column_ids.extend(selected.index)
        library_rows.append({"gsm": gsm, "n_raw_features": len(frame), "n_GEX_feature_rows": int(keep.sum()),
            "n_GEX_submitted_unique_symbols": len(set(submitted)), "n_canonical_NP_QC_symbols": len(canonical_unique),
            "n_common_output_features": len(genes), "n_noncommon_GEX_rows": len(names) - len(genes),
            "n_selected_cells": len(selected), "gene_order_equal_to_first_library": list(names) == list(first)})
        print("REAL_RAW_COUNTS_VERIFIED", gsm, list(counts.shape), "selected", list(result.shape), flush=True)
    audits = pd.DataFrame(audits).set_index("cell_id").loc[meta.index].reset_index()
    audits.to_csv(audit_dir / "NEW_INPUT_CELL_DONOR_BARCODE_AUDIT.tsv", sep="\t", index=False)
    pd.DataFrame(library_rows).to_csv(audit_dir / "RAW_GENE_UNIVERSE_BY_LIBRARY.tsv", sep="\t", index=False)
    if require_cached_umi and not audits.total_UMI_equal.all():
        raise ValueError("Raw UMI totals disagree with actual identity metadata; barcode audit preserved")
    combined = sparse.hstack(matrices, format="csc")
    order = pd.Index(column_ids).get_indexer(meta.index)
    if (order < 0).any():
        raise ValueError("Cannot restore independently seeded metadata order")
    combined = combined[:, order]
    with gzip.open(output_dir / "counts.mtx.gz", "wb") as handle:
        io.mmwrite(handle, combined, field="integer")
    pd.Series(genes).to_csv(output_dir / "genes.tsv", sep="\t", index=False, header=False)
    meta.to_csv(output_dir / "meta.tsv", sep="\t", index=False)
    return {"shape_genes_by_cells": list(combined.shape), "nnz": int(combined.nnz),
        "raw_UMI_equal_to_cached": int(audits.total_UMI_equal.sum()),
        "raw_features_equal_to_cached": int(audits.detected_genes_equal.sum()),
        "historical_symbol_policy_features_equal_to_cached": int(audits.historical_gene_count_equal.sum()),
        "output_common_gene_UMI_loss_total": int(audits.common_gene_UMI_loss.sum()),
        "common_gene_policy": "exact shared GEX features, first-library order, duplicate symbols made unique; no score-derived counts",
        "raw_gene_alias_changes": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--project-root")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--from-raw", action="store_true")
    mode.add_argument("--from-cached-inputs", action="store_true")
    parser.add_argument("--no-download", action="store_true", help="Require present real GEO triplets")
    args = parser.parse_args()
    text = args.config.read_text(encoding="utf-8")
    if args.config.suffix.lower() == ".json":
        config = json.loads(text)
    else:
        import yaml
        config = yaml.safe_load(text)
    root = data_path(REPO, args.project_root or config.get("project_root", "."))
    settings = config.get("raw_selection", {})
    output = writable(root, data_path(root, settings.get("output_dir", "02_sigxtalk_official_new/input")))
    audit = writable(root, data_path(root, settings.get("audit_dir", "02_sigxtalk_official_new")))
    raw_dir = data_path(root, settings.get("raw_dir", "data_cache/discovery_original_route/02_geo_downloads/GSE230808"))
    for directory in [output, audit, root / "00_audit", root / "logs"]:
        directory.mkdir(parents=True, exist_ok=True)
    manifest_path = resource_path(settings.get("sample_manifest", "resources/metadata/SIGXTALK_GSE230808_SOURCE_MANIFEST.tsv"))
    manifest = validate_manifest(manifest_path)
    proof = {"started_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "mode": "from-raw" if args.from_raw else "from-cached-inputs", "steps": [], "inputs": [],
        "sample_manifest_sha256": sha(manifest_path), "scientific_parameters": {
            "seed": SEED, "NP_quartiles": [.25, .75], "minimum_cells_per_tail": MIN_TAIL,
            "cap_per_donor_group": CAP, "expected_eligible_pairs": 6,
            "AF_identity_seed": SEED, "counts_source": "genuine GEO 10x matrices"}}
    python = settings.get("python", sys.executable)
    receipts = []
    def input_file(path, role):
        if not path.exists():
            raise FileNotFoundError(f"Actual {role} input missing: {path}")
        proof["inputs"].append({"role": role, "path": str(path), "bytes": path.stat().st_size, "sha256": sha(path)})
    def run(name, command):
        path = root / "logs" / ("raw_selection_" + name + ".log")
        started = time.time()
        with path.open("w") as handle:
            result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT)
        proof["steps"].append({"step": name, "command": command, "exit_code": result.returncode,
            "seconds": time.time() - started, "log_sha256": sha(path)})
        if result.returncode:
            raise RuntimeError("Raw preparation step failed: " + name + "; inspect " + str(path))
    try:
        if args.from_raw:
            receipts.extend(obtain_raw(manifest, raw_dir, root, not args.no_download))
            human_config_path = resource_path(settings.get("human_config", "config/human.example.yaml"))
            input_file(human_config_path, "human_discovery_configuration")
            human_text = human_config_path.read_text(encoding="utf-8")
            if human_config_path.suffix.lower() == ".json":
                human_config = json.loads(human_text)
            else:
                import yaml
                human_config = yaml.safe_load(human_text)
            discovery_config = human_config.get("discovery", {})
            discovery_workspace = writable(root, data_path(root, discovery_config.get("phase2_workspace", "data_cache/discovery_original_route")))
            rebuild_workspace = writable(root, data_path(root, discovery_config.get("rebuild_workspace", "data_cache/discovery_identity_rebuild")))
            flat = writable(root, discovery_workspace / "02_geo_downloads/GSE230808")
            flat.mkdir(parents=True, exist_ok=True)
            for row in manifest[manifest.tissue == "NP"].to_dict("records"):
                for suffix in SUFFIXES:
                    source = find_raw(raw_dir, row["gsm"], suffix)
                    destination = flat / source.name
                    if source.resolve() != destination.resolve():
                        if destination.exists():
                            if sha(source) != sha(destination):
                                raise ValueError("Conflicting discovery raw staging input")
                        else:
                            try:
                                os.link(source, destination)
                            except OSError:
                                shutil.copyfile(source, destination)
            run("discovery_QC_identity_UCell", [python, str(resource_path(settings.get("human_discovery_driver", "scripts/01_discovery_np_cnmf.py"))),
                "--config", str(human_config_path), "--project-root", str(root), "--step", "all"])
            np_meta_path = rebuild_workspace / "tables/discovery_cell_metadata.tsv.gz"
            score_path = rebuild_workspace / "tables/discovery_programs_UCell.tsv.gz"
            af_dir = writable(root, data_path(root, settings.get("af_identity_dir", "data_cache/af_identity")))
            af_meta_path = af_dir / "tables/AF_CELL_METADATA.tsv.gz"
            run("AF_QC_scrublet_fixed_identity", [settings.get("af_python", python),
                str(Path(__file__).with_name("support_af_identity.py")), "--raw-dir", str(raw_dir),
                "--out-dir", str(af_dir / "tables"), "--cache-dir", str(af_dir / "work"),
                "--sample-manifest", str(manifest_path), "--seed", str(SEED)])
            proof["steps"].append({"step": "program_fit", "status": "SKIPPED_FROZEN_K5_REUSED", "refit_performed": False})
        else:
            np_meta_path = data_path(root, settings["cached_discovery_metadata"])
            score_path = data_path(root, settings["cached_discovery_scores"])
            af_meta_path = data_path(root, settings["cached_af_metadata"])
            proof["steps"].append({"step": "QC_identity_UCell", "status": "SKIPPED_REUSED_ACTUAL_CELL_LEVEL_INPUTS"})
        for path, role in [(np_meta_path, "NP_metadata"), (score_path, "NP_UCell_scores"), (af_meta_path, "AF_metadata")]:
            input_file(path, role)
        common_fields = ["cell_id", "gsm", "donor_id", "n_counts", "n_features"]
        np_meta = read_cells(np_meta_path, common_fields + ["high_confidence_NP"])
        af_meta = read_cells(af_meta_path, common_fields + ["broad_stromal"])
        scores = read_cells(score_path, ["cell_id", "GEP2_top25_UCell"])
        validate_cell_sources(np_meta, manifest, "NP")
        validate_cell_sources(af_meta, manifest, "AF")
        try:
            selected, eligibility = select_cells(np_meta, scores, af_meta)
        except EligibilityError as error:
            error.eligibility.to_csv(output / "RECOMPUTED_DONOR_ELIGIBILITY.tsv", sep="\t", index=False)
            raise
        eligibility.to_csv(output / "RECOMPUTED_DONOR_ELIGIBILITY.tsv", sep="\t", index=False)
        donor_groups = selected.groupby(["donor_id", "gsm", "bridge_group"], observed=True).size().rename("n_cells").reset_index()
        donor_groups.to_csv(output / "NEW_DONOR_GROUP_COUNTS.tsv", sep="\t", index=False)
        selection = {"seed": SEED, "n_selected": len(selected), "eligible_donors": sorted(selected.donor_id.unique()),
            "n_selected_libraries": int(selected.gsm.nunique()), "selection_recomputed_independently": True,
            "criteria": "within-donor high-confidence NP GEP2 UCell Q1/Q4; min20/tail and q75>q25; cap150; matched broad-stromal AF; sorted-ID seeded sampling"}
        reference = settings.get("reference_selected_metadata")
        if reference:
            reference_path = data_path(root, reference)
            input_file(reference_path, "optional_private_selection_comparison_only")
            prior = read_cells(reference_path, ["cell_id", "bridge_group"])
            selection.update(reference_cells=len(prior), same_cell_set=set(prior.index) == set(selected.index),
                same_cell_order=list(prior.index) == list(selected.index))
            if selection["same_cell_set"]:
                selection["same_groups"] = bool((prior.loc[selected.index, "bridge_group"] == selected.bridge_group).all())
        dump(selection, output / "BARCODE_SELECTION_PROOF.json")
        selected_manifest = manifest[manifest.gsm.isin(selected.gsm.unique())]
        if not args.from_raw:
            receipts.extend(obtain_raw(selected_manifest, raw_dir, root, not args.no_download))
        summary = reconstruct_counts(selected, selected_manifest, raw_dir, output, audit,
            allow_intersection=bool(settings.get("allow_common_gene_intersection", False)),
            require_cached_umi=bool(settings.get("require_raw_umi_equal_to_metadata", True)))
        summary.update(selection, completed_utc=dt.datetime.now(dt.timezone.utc).isoformat())
        dump(summary, output / "NEW_RAW_INPUT_PROOF.json")
        for name in ["counts.mtx.gz", "genes.tsv", "meta.tsv"]:
            input_file(output / name, "new_selected_raw_contract_" + name)
        proof["status"] = "DONE_FROM_RAW" if args.from_raw else "DONE_RESELECTED_FROM_ACTUAL_CACHES_AND_REAL_RAW_COUNTS"
        proof["result"] = summary
        print("REAL_RAW_INPUT_COMPLETE", json.dumps(summary), flush=True)
    except Exception as error:
        proof["status"], proof["error"] = "PARTIAL_FAILED", str(error)
        raise
    finally:
        pd.DataFrame(receipts).to_csv(root / "00_audit/NEW_RAW_INPUT_FILES.tsv", sep="\t", index=False)
        proof["ended_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        dump(proof, root / "00_audit/RAW_SELECTION_EXECUTION.json")


if __name__ == "__main__":
    main()
