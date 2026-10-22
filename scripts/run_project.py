#!/usr/bin/env python3
"""Dispatch the real study modules with explicit cache reuse and failure receipts."""
from __future__ import annotations
import argparse
import csv
import importlib.util
import json
from pathlib import Path
import sys
from project_common import *

DEPTH_FILES = ["ROI_FIXED_PROVENANCE.tsv", "DISC_LEVEL_PROGRAM_AND_DEPTH.tsv", "DEPTH_SENSITIVITY_VALUES.tsv",
    "DEPTH_DISC_PAIR_COMPARISON.tsv", "DEPTH_RETENTION_BY_CAPTURE_REGION.tsv", "GEP2_SCORE_UMI_ASSOCIATION.tsv",
    "CANDIDATE_GENES_PER_DISC_REGION.tsv", "DEPTH_SENSITIVITY_REPORT.md", "SPATIAL_REVISED_INTERPRETATION.md"]


def import_script(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def spatial_context(cfg, raw=False):
    path = module_config(cfg, "spatial")
    values = load(path)
    for key in list(values):
        if key.endswith("_config"):
            values[key] = data_path(path.parent, values[key])
    root = cfg["_root"]
    for key in ["baseline_cache", "cached_spatial", "depth_output", "figure_output"]:
        values[key] = data_path(root, values[key])
    return values, values["baseline_cache" if raw else "cached_spatial"]


def protected_spatial_path(old_path, baseline, roi):
    """Support moved caches only through unambiguous recorded relative suffixes."""
    old = Path(old_path)
    if old.is_file():
        return old
    text = str(old_path).replace("\\", "/")
    for marker in ["prepared48/", "roi_validation/"]:
        if marker in text:
            candidate = baseline / text[text.index(marker):]
            if candidate.is_file():
                return candidate
    if old.name == roi.name:
        return roi
    candidate = baseline / old.name
    if candidate.is_file():
        return candidate
    raise FileNotFoundError("Cannot resolve a protected input from the completed depth check: " + text)


def verify_completed_depth(values, baseline):
    output = values["depth_output"]
    status_path = output / "DEPTH_ANALYSIS_STATUS.json"
    if not status_path.exists():
        if (output / "DEPTH_SENSITIVITY_VALUES.tsv").exists():
            raise RuntimeError("Existing depth draw lacks completion proof; preserve it for review, do not redraw automatically")
        return None
    status = load(status_path)
    required = {"status": "DONE_SINGLE_PREDECLARED_DEPTH_CHECK", "target_total_UMI": 100,
        "single_draw_only": True, "protected_old_inputs_unchanged": True,
        "config_sha256": sha(values["depth_config"]), "ROI_sha256": sha(values["roi_config"])}
    required["random_seed"] = load(values["depth_config"])["seed"]
    if any(status.get(key) != value for key, value in required.items()):
        raise RuntimeError("Existing fixed-depth completion/config differs; preserve the original draw and investigate")
    files = require_files([output / name for name in DEPTH_FILES])
    audit = output / "PRIVATE_FROZEN_INPUT_HASH_AUDIT.tsv"
    with audit.open(encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if not rows:
        raise ValueError("Completed fixed-depth check has no protected-input hash audit")
    for row in rows:
        current = protected_spatial_path(row["path"], baseline, values["roi_config"])
        if row["sha256_before"] != row["sha256_after"] or sha(current) != row["sha256_after"]:
            raise ValueError("Protected input changed since the single fixed-depth draw: " + str(current))
    return {"status": "SKIPPED_VERIFIED_EXISTING_SINGLE_DEPTH_DRAW", "depth_status": status,
        "outputs": files, "single_draw_repeated": False}


def b_paths(cfg):
    path, values = resolved_b_config(cfg)
    module = data_path(cfg["_root"], values["module_dir"])
    return path, values, module


def verify_b_completion(cfg, mode="formal"):
    config_path, values, module = b_paths(cfg)
    folder = module / ("formal_new" if mode == "formal" else values["smoke"].get("directory_name", "technical_smoke"))
    final = folder / "NEW_OFFICIAL_FINAL_STAGE_STATUS.json"
    if not final.exists():
        return None
    status = load(final)
    if not status.get("complete_R_python_R") or (mode == "formal" and not (
            status.get("formal_design") and status.get("new_result_not_legacy_relabel"))):
        raise ValueError("Existing B status does not prove a new completed official run")
    provenance = load(folder / "audit/NEW_INPUT_PROVENANCE.json")
    if provenance["config_sha256"] != sha(config_path):
        raise ValueError("B configuration differs from completed input provenance")
    hgnn = load(folder / "audit/HGNN_STATUS.json")
    expected_epochs = 100 if mode == "formal" else values["smoke"]["epochs"]
    if hgnn.get("epochs_completed") != expected_epochs or hgnn.get("seed") != 20261017 or not hgnn.get("new_computation"):
        raise ValueError("B epoch/seed/new-computation receipt is incomplete")
    prediction = Path(hgnn["prediction_file"])
    if not prediction.is_absolute():
        prediction = cfg["_root"] / prediction
    if sha(prediction) != hgnn["prediction_sha256"] or not hgnn.get("all_tested_rows", 0):
        raise ValueError("New complete prediction file hash/row proof is missing")
    check_b_input_hashes(cfg, provenance, values)
    official = import_script("pinned_b_verifier", repo_path("scripts/run_sigxtalk_pipeline.py"))
    official.verify_upstream(data_path(cfg["_root"], values["official_source"]), repo_path("scripts/sigxtalk_upstream_lock.json"))
    require_files([folder / "audit/NEW_R_PYTHON_R_COMPLETION.json", prediction])
    return {"status": "SKIPPED_VERIFIED_CURRENT_NEW_OFFICIAL_RUN", "mode": mode,
        "run_class": status.get("run_class"), "completion": status, "config_sha256": sha(config_path)}


def check_b_input_hashes(cfg, provenance, values):
    hashes = provenance.get("new_input_file_hashes", {})
    if provenance.get("input_kind") == "real_GSE230808_raw_counts":
        inputs = data_path(cfg["_root"], values["input_dir"])
        pairs = data_path(cfg["_root"], values["pairing_tsv"])
        for key, path in [("raw_matrix_sha256", inputs / "counts.mtx.gz"), ("raw_genes_sha256", inputs / "genes.tsv"),
                          ("cell_identity_sha256", inputs / "meta.tsv"), ("patient_pairing_sha256", pairs)]:
            if hashes.get(key) != sha(path):
                raise ValueError("Completed B trained input changed: " + str(path))
        if provenance.get("started_from_old_pathway_output"):
            raise ValueError("Old pathway output cannot establish a new B run")
    elif "pbmc_rds_sha256" in hashes:
        if sha(data_path(cfg["_root"], values["pbmc_rds"])) != hashes["pbmc_rds_sha256"]:
            raise ValueError("The official technical-demo RDS changed")
    else:
        raise ValueError("B input hash provenance is unrecognised")


def existing_raw_reuse(cfg, raw_config):
    """Recheck selection/source hashes without rewriting trained gzip inputs."""
    import pandas as pd
    raw = load(raw_config).get("raw_selection", {})
    root = cfg["_root"]
    output = data_path(root, raw.get("output_dir", "02_sigxtalk_official_new/input"))
    present = [output / name for name in ["counts.mtx.gz", "genes.tsv", "meta.tsv"]]
    if not any(path.exists() for path in present):
        return None
    require_files(present + [output / "NEW_RAW_INPUT_PROOF.json"])
    helper = import_script("own_raw_selection", repo_path("scripts/prepare_sigxtalk_raw.py"))
    np_path = data_path(root, raw["cached_discovery_metadata"])
    score_path = data_path(root, raw["cached_discovery_scores"])
    af_path = data_path(root, raw["cached_af_metadata"])
    np_meta = helper.read_cells(np_path, ["cell_id", "donor_id", "gsm", "n_counts", "n_features", "high_confidence_NP"])
    scores = helper.read_cells(score_path, ["cell_id", "GEP2_top25_UCell"])
    af_meta = helper.read_cells(af_path, ["cell_id", "donor_id", "gsm", "n_counts", "n_features", "broad_stromal"])
    selected, _ = helper.select_cells(np_meta, scores, af_meta)
    actual = helper.read_cells(output / "meta.tsv", ["cell_id", "bridge_group"])
    if selected.index.tolist() != actual.index.tolist() or not (selected.bridge_group == actual.bridge_group).all():
        raise ValueError("Existing raw selection differs from independently resampled actual cell-level inputs")
    source_receipts = root / "00_audit/NEW_RAW_INPUT_FILES.tsv"
    sources = pd.read_csv(source_receipts, sep="\t")
    if sources.empty or not {"path", "sha256"}.issubset(sources):
        raise ValueError("Actual raw-source hash receipt is unavailable")
    accession = "gsm" if "gsm" in sources else "accession" if "accession" in sources else None
    if accession is None:
        raise ValueError("Raw source receipts need their actual sample accession")
    sources = sources[sources[accession].isin(selected.gsm.unique())]
    if set(sources[accession]) != set(selected.gsm.unique()):
        raise ValueError("Raw source receipts do not cover all selected libraries")
    for gsm, group in sources.groupby(accession):
        for suffix in ["matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"]:
            if sum(Path(path).name.endswith(suffix) for path in group.path) != 1:
                raise ValueError("Raw source receipt lacks an unambiguous triplet: " + gsm)
    for row in sources.to_dict("records"):
        if sha(Path(row["path"])) != row["sha256"]:
            raise ValueError("Public raw count source changed: " + row["path"])
    private_receipt = root / "00_audit/RAW_SELECTION_EXECUTION.json"
    contract_proven = False
    if private_receipt.exists():
        saved = load(private_receipt)
        by_role = {row["role"]: row["sha256"] for row in saved.get("inputs", [])}
        contract_proven = all(by_role.get("new_selected_raw_contract_" + path.name) == sha(path) for path in present)
    if not contract_proven:
        _, values, module = b_paths(cfg)
        prov_path = module / "formal_new/audit/NEW_INPUT_PROVENANCE.json"
        if not prov_path.exists():
            raise ValueError("Existing selected input lacks its own hash receipt or new B input provenance; preserve for review")
        check_b_input_hashes(cfg, load(prov_path), values)
    return {"status": "SKIPPED_VERIFIED_CURRENT_RAW_INPUTS", "independent_selection_recomputed": True,
        "counts_rewritten": False, "files": require_files(present), "source_hashes_rechecked": len(sources),
        "cached_identity_input_hashes": {name: sha(path) for name, path in [("NP", np_path), ("UCell", score_path), ("AF", af_path)]}}


def materialize_pairing(cfg, raw_config):
    import pandas as pd
    _, b, _ = b_paths(cfg)
    raw = load(raw_config).get("raw_selection", {})
    manifest = pd.read_csv(repo_path(raw.get("sample_manifest", "resources/metadata/SIGXTALK_GSE230808_SOURCE_MANIFEST.tsv")), sep="\t")
    pairs = manifest.pivot(index="donor_id", columns="tissue", values="gsm").reset_index().rename(columns={"AF": "AF_gsm", "NP": "NP_gsm"})
    path = data_path(cfg["_root"], b["pairing_tsv"])
    if path.exists():
        existing = pd.read_csv(path, sep="\t")
        if len(existing) not in [6, 8] or existing.donor_id.duplicated().any():
            raise ValueError("Existing pairing must contain the six selected or eight verified unique donors")
        expected = pairs.set_index("donor_id")
        for row in existing.itertuples():
            if row.donor_id not in expected.index or row.AF_gsm != expected.loc[row.donor_id, "AF_gsm"] or row.NP_gsm != expected.loc[row.donor_id, "NP_gsm"]:
                raise ValueError("Existing B pairing file conflicts with public exact donor pairing")
    else:
        writable(cfg["_root"], path)
        path.parent.mkdir(parents=True, exist_ok=True)
        pairs.to_csv(path, sep="\t", index=False)
    if b.get("state_thresholds_tsv"):
        thresholds = data_path(cfg["_root"], b["state_thresholds_tsv"])
        source = data_path(cfg["_root"], raw.get("output_dir", "02_sigxtalk_official_new/input")) / "RECOMPUTED_DONOR_ELIGIBILITY.tsv"
        actual = pd.read_csv(source, sep="\t")[["donor_id", "q25", "q75"]]
        if not thresholds.exists():
            writable(cfg["_root"], thresholds)
            thresholds.parent.mkdir(parents=True, exist_ok=True)
            actual.to_csv(thresholds, sep="\t", index=False)
        else:
            import numpy as np
            prior = pd.read_csv(thresholds, sep="\t").set_index("donor_id")
            comparison = actual.set_index("donor_id")
            if not set(prior.index).issubset(comparison.index) or not np.allclose(
                    prior[["q25", "q75"]], comparison.loc[prior.index, ["q25", "q75"]], equal_nan=True):
                raise ValueError("Existing B score thresholds differ from actual independently recomputed quartiles")


def stage_human(cfg, mode):
    root = cfg["_root"]
    config = module_config(cfg, "human")
    steps = [execute(cfg, "human", [science_python(cfg), repo_path("scripts/02_crosscohort_np.py"),
        "--config", config, "--project-root", root, mode, "--skip-figures"])]
    values = cfg.get("human", {})
    candidate = data_path(root, values.get("candidate_output", "01_human_programs/TU_CANDIDATE_GENE_BY_DONOR_IDENTITY.tsv"))
    steps.append(execute(cfg, "human_candidate_expression", [science_python(cfg), repo_path("scripts/extract_tu_candidate_expression.py"),
        "--tu-sparse-root", data_path(root, values.get("tu_sparse_root", "data_cache/Tu_sparse")),
        "--program-mapping", data_path(root, values.get("candidate_mapping", "01_human_programs/PROGRAM_GENE_MAPPING.tsv")), "--out", candidate]))
    return {"status": "DONE_REGENERATED_HUMAN_TABLES", "executions": steps, "outputs": require_files([
        root / "01_human_programs/PROGRAM_SCORE_DONOR_SUMMARY.tsv", root / "01_human_programs/GEP2_DONOR_SUPPORT.tsv", candidate])}


def stage_raw(cfg, mode):
    config = module_config(cfg, "raw_selection")
    raw = load(config).get("raw_selection", {})
    _, b, _ = b_paths(cfg)
    if data_path(cfg["_root"], raw.get("output_dir", "02_sigxtalk_official_new/input")) != data_path(cfg["_root"], b["input_dir"]):
        raise ValueError("raw_selection.output_dir must match the frozen B input_dir; no silent input-path substitution")
    if mode == "--from-cached-inputs":
        previous = existing_raw_reuse(cfg, config)
        if previous:
            materialize_pairing(cfg, config)
            return previous
    receipt = execute(cfg, "raw_selection", [science_python(cfg), repo_path("scripts/prepare_sigxtalk_raw.py"),
        "--config", config, "--project-root", cfg["_root"], mode])
    materialize_pairing(cfg, config)
    return {"status": "DONE_REAL_RAW_INPUT_PREPARATION", "execution": receipt}


def stage_spatial(cfg, mode):
    values, baseline = spatial_context(cfg, mode == "--from-raw")
    existing = verify_completed_depth(values, baseline)
    if existing:
        return existing
    command = [science_python(cfg), repo_path("scripts/07_prepare_spatial.py"), "--config", module_config(cfg, "spatial"),
        "--project-root", cfg["_root"], mode]
    if mode == "--from-raw":
        command.append("--download-public")
    steps = [execute(cfg, "spatial_baseline", command)]
    command = [science_python(cfg), repo_path("scripts/08_spatial_depth_check.py"), "--config", module_config(cfg, "spatial"), "--project-root", cfg["_root"]]
    if mode == "--from-raw":
        command.append("--from-raw")
    steps.append(execute(cfg, "spatial_single_depth", command))
    result = verify_completed_depth(values, baseline)
    result.update(status="DONE_ONE_NEW_FIXED_DEPTH_DRAW", executions=steps)
    return result


def stage_sigxtalk(cfg):
    config, values, module = b_paths(cfg)
    steps = []
    for mode in ["smoke", "formal"]:
        completed = verify_b_completion(cfg, mode)
        if completed:
            steps.append(completed)
            continue
        folder = module / ("formal_new" if mode == "formal" else values["smoke"].get("directory_name", "technical_smoke"))
        hgnn = folder / "audit/HGNN_STATUS.json"
        resume = "all"
        if hgnn.exists():
            previous = load(hgnn)
            provenance = load(folder / "audit/NEW_INPUT_PROVENANCE.json")
            if previous.get("status") != "DONE" or provenance.get("config_sha256") != sha(config):
                raise ValueError("Existing partial/changed model requires review; no automatic retraining")
            check_b_input_hashes(cfg, provenance, values)
            resume = "prs"
        steps.append(execute(cfg, "sigxtalk_" + mode, [science_python(cfg), repo_path("scripts/run_sigxtalk_pipeline.py"),
            "--config", config, "--mode", mode, "--stage", resume]))
        verify_b_completion(cfg, mode)
    steps.append(execute(cfg, "sigxtalk_reports", [science_python(cfg), repo_path("scripts/06_sigxtalk_reports.py"),
        "--config", config, "--mode", "formal", "--plot-official-metrics"]))
    finalize = [science_python(cfg), repo_path("scripts/finalize_new_sigxtalk.py"), "--config", config]
    options = cfg.get("sigxtalk_finalize", {})
    if options.get("legacy_root"):
        finalize += ["--legacy-root", data_path(cfg["_root"], options["legacy_root"])]
    for record in options.get("environment_records", []):
        finalize += ["--environment-record", data_path(cfg["_root"], record)]
    steps.append(execute(cfg, "sigxtalk_finalize", finalize))
    return {"status": "DONE_NEW_OFFICIAL_CLOSURE_VERIFIED", "steps": steps}


def stage_integrate(cfg):
    config, _, _ = b_paths(cfg)
    _, baseline = spatial_context(cfg, cfg.get("_raw_mode", False))
    values = cfg.get("human", {})
    command = [science_python(cfg), repo_path("scripts/09_integrate_evidence.py"), "--config", module_config(cfg, "integration"),
        "--project-root", cfg["_root"], "--tu-candidate-table", data_path(cfg["_root"], values.get("candidate_output", "01_human_programs/TU_CANDIDATE_GENE_BY_DONOR_IDENTITY.tsv")),
        "--spatial-input", baseline, "--sigxtalk-config", config]
    options = cfg.get("integration", {})
    if options.get("github_proof"):
        command += ["--github-proof", data_path(cfg["_root"], options["github_proof"])]
    if options.get("allow_incomplete_b", False):
        command.append("--allow-incomplete-b")
    receipt = execute(cfg, "integrate", command)
    status = load(cfg["_root"] / "04_integrated/INTEGRATION_STATUS.json")
    if status.get("status") != "DONE_NEW_FORMAL_INTEGRATION":
        raise RuntimeError("Integration remains partial; draft outputs preserved without full-pipeline success")
    return {"status": "DONE_NEW_FORMAL_INTEGRATION", "execution": receipt}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/pipeline.example.yaml")
    parser.add_argument("--project-root")
    parser.add_argument("--stage", choices=["all", "human", "raw-selection", "sigxtalk", "spatial", "integrate", "figures", "audit"], default="all")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--from-raw", action="store_true")
    mode.add_argument("--from-cached-inputs", action="store_true")
    args = parser.parse_args()
    if args.stage not in ["figures", "audit", "integrate"] and not (args.from_raw or args.from_cached_inputs):
        parser.error("Explicit --from-raw or --from-cached-inputs is required for computation")
    cfg = pipeline(args.config, args.project_root)
    cfg["_raw_mode"] = args.from_raw
    cfg["_root"].mkdir(parents=True, exist_ok=True)
    mode = "--from-raw" if args.from_raw else "--from-cached-inputs"
    proof = {"started_utc": now(), "stage_requested": args.stage, "mode": mode, "config_sha256": sha(cfg["_config"]),
        "status": "RUNNING", "steps": [], "frozen_K5_refit_requested": False}
    plan = ["download", "raw-selection", "human", "spatial", "sigxtalk", "figures", "integrate", "audit"] if args.stage == "all" else [args.stage]
    try:
        for stage in plan:
            print("PROJECT_STAGE", stage, flush=True)
            if stage == "download":
                if mode == "--from-cached-inputs":
                    result = execute(cfg, "download_cache_skip", [science_python(cfg), repo_path("scripts/00_download_public_data.py"), "--config", cfg["_config"], "--project-root", cfg["_root"], "--from-cached-inputs"])
                else:
                    result = execute(cfg, "download", [science_python(cfg), repo_path("scripts/00_download_public_data.py"), "--config", cfg["_config"], "--project-root", cfg["_root"]])
            elif stage == "raw-selection": result = stage_raw(cfg, mode)
            elif stage == "human": result = stage_human(cfg, mode)
            elif stage == "spatial": result = stage_spatial(cfg, mode)
            elif stage == "sigxtalk": result = stage_sigxtalk(cfg)
            elif stage == "integrate": result = stage_integrate(cfg)
            else:
                script = "10_make_figures.py" if stage == "figures" else "11_quality_audit.py"
                command = [science_python(cfg), repo_path("scripts/" + script), "--config", cfg["_config"], "--project-root", cfg["_root"]]
                if args.from_raw:
                    command.append("--from-raw")
                result = execute(cfg, stage, command)
            proof["steps"].append({"stage": stage, "result": result})
            dump(proof, cfg["_root"] / "00_audit/PROJECT_EXECUTION_STATUS.json")
        proof["status"] = "DONE_ALL_REQUIRED_NEW_OUTPUTS_VERIFIED" if args.stage == "all" else "DONE_REQUESTED_STAGE"
    except Exception as error:
        proof["status"], proof["error"] = "PARTIAL_FAILED", str(error)
        raise
    finally:
        proof["ended_utc"] = now()
        dump(proof, cfg["_root"] / "00_audit/PROJECT_EXECUTION_STATUS.json")


if __name__ == "__main__":
    main()
