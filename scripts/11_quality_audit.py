#!/usr/bin/env python3
"""Light source/actual-output audit; CI never downloads data or trains models."""
from __future__ import annotations
import argparse
import ast
import csv
import json
from pathlib import Path
import re
from project_common import *

REQUIRED_SCRIPTS = ["run_project.py", "00_download_public_data.py", "01_discovery_np_cnmf.py", "02_crosscohort_np.py",
    "extract_tu_candidate_expression.py", "10_make_human_figures.py", "prepare_sigxtalk_raw.py", "support_af_identity.py", "run_sigxtalk_pipeline.py", "03_prepare_sigxtalk.R",
    "05_sigxtalk_prs.R", "06_sigxtalk_reports.py", "07_prepare_spatial.py", "08_spatial_depth_check.py",
    "09_integrate_evidence.py", "10_make_figures.py", "10_plot_new_sigxtalk_figure3.py", "plot_figure4.py", "11_quality_audit.py"]


def source_files(root):
    folders = ["scripts", "config", "configs", "resources", "envs", "tests", "docs", ".github"]
    files = []
    for folder in folders:
        files.extend(path for path in (root / folder).rglob("*") if path.is_file() and not any(part in {"__pycache__", "node_modules", ".smoke", "analysis", "outputs"} for part in path.relative_to(root).parts) and path.suffix != ".pyc")
    files.extend(path for path in root.iterdir() if path.is_file() and path.suffix.lower() in [".md", ".cff", ".sh", ".toml"])
    return sorted(set(files))


def audit_source(root):
    required = [root / "scripts" / name for name in REQUIRED_SCRIPTS]
    require_files(required + [root / "run_all.sh", root / "resources/frozen_k5/PROGRAM_LOCK.json", root / "scripts/sigxtalk_upstream_lock.json"])
    files, parsed = source_files(root), 0
    forbidden_suffixes = {".h5ad", ".rds", ".hdf5", ".pt", ".pth", ".pkl", ".rdata", ".zip"}
    findings = []
    for path in files:
        relative = path.relative_to(root).as_posix()
        if path.suffix.lower() in forbidden_suffixes or "private" in path.name.lower():
            findings.append(relative + ": runtime/private object in source tree")
        if path.suffix.lower() in [".tsv", ".csv", ".gz"] and re.search(r"gene_spectra|consensus_spectra|all_gene_ranking|cell_usage", path.name, re.I):
            findings.append(relative + ": fitted model weights in public tree")
        if path.suffix.lower() in [".tsv", ".csv"]:
            with path.open(encoding="utf-8") as stream:
                if "gene_spectra_score" in stream.readline().split():
                    findings.append(relative + ": fitted coefficient column in public table")
        if path.suffix == ".py":
            ast.parse(path.read_text(encoding="utf-8"), filename=relative)
            parsed += 1
        if path.suffix.lower() in [".py", ".r", ".md", ".yaml", ".yml", ".json", ".sh", ".cff"]:
            text = path.read_text(encoding="utf-8")
            if re.search(r"[A-Za-z]:[/\\]Users[/\\]|/home/asset(?:/|\\)|ghp_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,}|BEGIN [A-Z ]*PRIVATE KEY", text):
                findings.append(relative + ": personal path or credential pattern")
    if findings:
        raise ValueError("Source audit failed: " + "; ".join(findings))
    frozen = json.loads((root / "resources/frozen_k5/PROGRAM_LOCK.json").read_text())
    programs = frozen["programs"]
    if len(programs) != 5 or any(len(genes) != 25 or len(set(genes)) != 25 for genes in programs.values()):
        raise ValueError("Frozen program contract is not five exact 25-gene signatures")
    if sha(root / "resources/frozen_k5/PROGRAM_LOCK.json") != "12905c6475cd4f3dd07c660efeee15e60d59f592c8f74c6cc60202f549110da2":
        raise ValueError("The original retained K5 lock bytes changed")
    return {"status": "PASS_SOURCE_STRUCTURE", "python_files_parsed": parsed,
        "source_files": len(files), "frozen_lock_sha256": sha(root / "resources/frozen_k5/PROGRAM_LOCK.json"),
        "large_raw_downloads_or_training_executed": False,
        "third_party_source_vendoring_audit": "explicit allowlists must also be reviewed by maintainer"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/pipeline.example.yaml")
    parser.add_argument("--project-root")
    parser.add_argument("--repo-root", type=Path, help="Source tree to audit; useful in lightweight CI")
    parser.add_argument("--source-only", action="store_true")
    parser.add_argument("--from-raw", action="store_true")
    args = parser.parse_args()
    source = audit_source((args.repo_root or REPO).resolve())
    if args.source_only:
        print(json.dumps(source, indent=2))
        return
    cfg = pipeline(args.config, args.project_root)
    root = cfg["_root"]
    from run_project import verify_b_completion, spatial_context, verify_completed_depth
    formal = verify_b_completion(cfg)
    spatial, baseline = spatial_context(cfg, args.from_raw)
    depth = verify_completed_depth(spatial, baseline)
    if formal is None or depth is None:
        raise RuntimeError("New B/C completion is missing; no all-stage quality pass is generated")
    integration = load(root / "04_integrated/INTEGRATION_STATUS.json")
    if integration.get("status") != "DONE_NEW_FORMAL_INTEGRATION" or integration.get("old_SigXTalk_rows_used"):
        raise ValueError("Integration is incomplete or uses historical SigXTalk results")
    figures = load(root / "05_figures/FIGURE_EXPORT_STATUS.json")
    if figures.get("main_figure_count") != 4:
        raise ValueError("Expected four main computational figures")
    require_files([row["path"] for row in figures["outputs"]])
    for row in figures["outputs"]:
        if sha(row["path"]) != row["sha256"]:
            raise ValueError("An exported figure changed after its receipt")
    required = [root / "01_human_programs/PROGRAM_SCORE_DONOR_SUMMARY.tsv",
        root / "01_human_programs/GEP2_DONOR_SUPPORT.tsv", root / "04_integrated/MAIN_CANDIDATE_EVIDENCE_MATRIX.tsv",
        root / "04_integrated/FINAL_CANDIDATE_AND_WETLAB_DECISION.md", root / "06_manuscript/RESULTS_COMPUTATION_ACTUALLY_EXECUTED.md"]
    require_files(required)
    records = []
    for folder in ["01_human_programs", "02_sigxtalk_official_new/formal_new/tables", "03_spatial_gse287499", "04_integrated", "05_figures", "06_manuscript"]:
        for path in sorted((root / folder).rglob("*")):
            if path.is_file() and path.suffix.lower() in [".tsv", ".json", ".md", ".png", ".pdf", ".svg"]:
                records.append({"relative_path": path.relative_to(root).as_posix(), "bytes": path.stat().st_size, "sha256": sha(path)})
    manifest = root / "00_audit/PIPELINE_OUTPUT_MANIFEST.tsv"
    manifest.parent.mkdir(parents=True, exist_ok=True)
    with manifest.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["relative_path", "bytes", "sha256"], delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    status = {"status": "PASS_COMPUTATIONAL_SOURCE_AND_OUTPUT_AUDIT", "finished_utc": now(), "source": source,
        "output_manifest_sha256": sha(manifest), "manifest_rows": len(records), "B_run_class": formal.get("run_class"),
        "single_spatial_depth_draw_repeated": False, "clinical_grade_P_values_created": False,
        "visual_QA": "requires separately recorded actual-image inspection", "wetlab_results_generated": False,
        "github_release": integration.get("github_code_release", "UNVERIFIED"),
        "full_clean_machine_from_raw_validated": False}
    dump(status, root / "00_audit/QUALITY_AUDIT_STATUS.json")
    print(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
