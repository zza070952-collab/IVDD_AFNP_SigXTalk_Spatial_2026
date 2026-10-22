#!/usr/bin/env python3
"""Reproducible raw/cached human module with explicit reuse and execution receipts."""
import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

import pandas as pd
from human_config import load_config, data_path, resource


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--project-root")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--from-raw", action="store_true")
    mode.add_argument("--from-cached-inputs", action="store_true")
    parser.add_argument("--skip-figures", action="store_true")
    args = parser.parse_args()
    config = load_config(args.config, args.project_root)
    root, repo = config["_root"], config["_repo"]
    out = root / "01_human_programs"
    for path in [root / "00_audit", root / "logs", out]:
        path.mkdir(parents=True, exist_ok=True)
    started = dt.datetime.now(dt.timezone.utc).isoformat()
    proof = {"started_utc": started, "mode": "from-raw" if args.from_raw else "from-cached-inputs", "steps": [], "inputs": []}
    python = config.get("executables", {}).get("human_python", sys.executable)
    environment = dict(__import__("os").environ)
    addons = config.get("executables", {}).get("python_addons")
    if addons:
        environment["PYTHONPATH"] = str(data_path(config, addons)) + __import__("os").pathsep + environment.get("PYTHONPATH", "")
    def run(name, command):
        t = time.time()
        log_path = root / "logs" / ("human_" + name + ".log")
        with log_path.open("w") as log:
            result = subprocess.run(command, env=environment, stdout=log, stderr=subprocess.STDOUT)
        proof["steps"].append({"step": name, "status": "DONE" if result.returncode == 0 else "FAILED", "exit_code": result.returncode, "seconds": time.time()-t, "log_sha256": sha(log_path), "command": command})
        if result.returncode:
            raise RuntimeError(f"Human step {name} failed; inspect {log_path}")
    def cache_input(path, role):
        if not path.exists():
            raise FileNotFoundError(f"Required {role} input unavailable: {path}; do not reconstruct it from result summaries")
        proof["inputs"].append({"role": role, "path": str(path), "bytes": path.stat().st_size, "sha256": sha(path)})
    try:
        common = ["--config", str(Path(args.config).resolve()), "--project-root", str(root)]
        if args.from_raw:
            run("download", [python, str(repo / "scripts/00_download_human_public_data.py"), *common])
            run("discovery_raw_qc_identity_ucell", [python, str(repo / "scripts/01_discovery_np_cnmf.py"), *common, "--step", "all"])
            rebuilt = data_path(config, config["discovery"]["rebuild_workspace"])
            discovery_metadata = rebuilt / "tables/discovery_cell_metadata.tsv.gz"
            discovery_scores = rebuilt / "tables/discovery_programs_UCell.tsv.gz"
            discovery_origin = "reconstructed_official_UCell"
        else:
            discovery_metadata = data_path(config, config["discovery"]["cached_metadata"])
            discovery_scores = data_path(config, config["discovery"]["cached_scores"])
            discovery_origin = "verified_legacy_UCell"
            run("frozen_program_audit", [python, str(repo / "scripts/01_discovery_np_cnmf.py"), *common, "--step", "audit-frozen"])
        cache_input(discovery_metadata, "discovery_cell_metadata")
        cache_input(discovery_scores, "discovery_cell_program_scores")
        # Cache mode may reuse actual cell-level Tu inputs, never summary-only scores.
        reused = config.get("tu", {}).get("cached_tables_dir")
        if args.from_cached_inputs and reused:
            source = data_path(config, reused)
            required = ["TU_CELL_METADATA.tsv.gz", "TU_PROGRAM_AUCELL_CELL_SCORES.tsv.gz", "PROGRAM_GENE_MAPPING.tsv",
                "PROGRAM_GENE_COVERAGE.tsv", "TU_MARKER_DOTPLOT_SOURCE.tsv", "TU_PROGRAM_GENE_DOTPLOT_SOURCE.tsv",
                "TU_MATRIX_AND_QC_MANIFEST.tsv", "COHORT_AND_DONOR_TABLE.tsv", "AUCELL_RUN_METADATA.json",
                "AUCELL_EFFECTIVE_PROGRAM_TARGETS.tsv", "AUCELL_EFFECTIVE_PROGRAM_COVERAGE.tsv"]
            for name in required:
                path = source / name
                cache_input(path, "tu_actual_cache_" + name)
                if path.resolve() != (out / name).resolve():
                    shutil.copyfile(path, out / name)
            shutil.copyfile(resource(config, "programs", "resources/frozen_k5/FROZEN_K5_PROGRAMS.tsv"), out / "FROZEN_K5_PROGRAMS.tsv")
            proof["steps"].append({"step": "Tu_QC_identifiers_AUCell", "status": "SKIPPED_REUSED_ACTUAL_CELL_LEVEL_INPUTS", "recalculated": False})
        else:
            tu = config["tu"]
            run("tu_parse_qc_identifiers", [python, str(repo / "scripts/human/parse_tu.py"), "--work", str(root),
                "--input-dir", str(data_path(config, tu["raw_dir"])),
                "--frozen-signatures", str(resource(config, "programs", "resources/frozen_k5/FROZEN_K5_PROGRAMS.tsv")),
                "--crosswalk", str(resource(config, "tu_crosswalk", "resources/metadata/TU_GEO_PAPER_GRADE_CROSSWALK.tsv")),
                "--alias-map", str(resource(config, "hgnc_alias_candidates", "resources/hgnc/ALIAS_CANDIDATES.tsv")),
                "--min-features", str(tu.get("min_features", 200)), "--max-mito-fraction", str(tu.get("max_mito_fraction", .2))])
            run("tu_official_aucell", [python, str(repo / "scripts/human/score_tu_aucell.py"), "--work", str(root),
                "--n-up", str(tu.get("n_up", 1500)), "--tmin", str(tu.get("tmin", 5)),
                "--bsize", str(tu.get("bsize", 256)), "--threads", str(tu.get("threads", 2))])
        run("donor_summary", [python, str(repo / "scripts/human/summarize.py"), "--work", str(root),
            "--tu-score-method", "AUCell", "--discovery-metadata", str(discovery_metadata), "--discovery-scores", str(discovery_scores),
            "--discovery-score-origin", discovery_origin])
        support = [python, str(repo / "scripts/human/support_tables.py"), "--work", str(root),
            "--frozen-dir", str(resource(config, "frozen_dir", "resources/frozen_k5")),
            "--sample-manifest", str(resource(config, "discovery_samples", "resources/metadata/DISCOVERY_SAMPLE_MANIFEST.tsv"))]
        if config.get("interpretation", {}).get("review_decisions"):
            support += ["--review-decisions", str(resource(config, "review_decisions", config["interpretation"]["review_decisions"]))]
        run("support_tables_and_methods", support)
        if not args.skip_figures:
            command = [python, str(repo / "scripts/human/plot_figures12.py"), "--work", str(root)]
            af = config.get("discovery", {}).get("af_donor_counts")
            if af:
                command += ["--af-donor-counts", str(data_path(config, af))]
            run("figures1_2", command)
        proof["status"] = "DONE_FROM_RAW_WITH_FROZEN_K5" if args.from_raw else "DONE_FROM_CACHED_CELL_LEVEL_INPUTS"
    except Exception as error:
        proof["status"] = "PARTIAL_FAILED"
        proof["error"] = str(error)
        raise
    finally:
        proof["ended_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        (root / "00_audit/HUMAN_MODULE_EXECUTION.json").write_text(json.dumps(proof, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
