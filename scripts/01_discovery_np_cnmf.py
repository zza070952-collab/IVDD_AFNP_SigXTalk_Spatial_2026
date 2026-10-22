#!/usr/bin/env python3
"""Executable historical discovery route; cNMF rerun is explicit, never default."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import pandas as pd
from human_config import load_config, data_path, resource
from audit_frozen_models import audit_private_models


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--project-root")
    parser.add_argument("--step", choices=["audit-frozen", "raw-qc", "identity", "ucell", "cnmf", "all"], default="audit-frozen")
    parser.add_argument("--execute-historical-cnmf", action="store_true")
    args = parser.parse_args()
    config = load_config(args.config, args.project_root)
    root = config["_root"]
    private_models = audit_private_models(config, root)
    discovery_config = config.get("discovery", {})
    discovery = data_path(config, discovery_config.get("phase2_workspace", "data_cache/discovery_original_route"))
    rebuilt = data_path(config, discovery_config.get("rebuild_workspace", "data_cache/discovery_identity_rebuild"))
    (root / "00_audit").mkdir(parents=True, exist_ok=True)
    if args.step != "audit-frozen":
        for workspace in [discovery, rebuilt]:
            if not workspace.resolve().is_relative_to(root.resolve()):
                raise ValueError("Writable reconstruction workspaces must remain inside the new project root; cached source inputs can remain outside")
        for path in [discovery / "01_metadata_audit", discovery / "03_discovery_NP",
                     discovery / "04_cNMF", discovery / "logs", rebuilt / "work", rebuilt / "tables", rebuilt / "envs"]:
            path.mkdir(parents=True, exist_ok=True)
        manifest = resource(config, "discovery_samples", "resources/metadata/DISCOVERY_SAMPLE_MANIFEST.tsv")
        shutil.copyfile(manifest, discovery / "01_metadata_audit/sample_manifest.tsv")
    lock_file = resource(config, "lock", "resources/frozen_k5/PROGRAM_LOCK.json")
    lock = json.loads(lock_file.read_text())
    frozen = pd.read_csv(resource(config, "programs", "resources/frozen_k5/FROZEN_K5_PROGRAMS.tsv"), sep="\t")
    for name, genes in lock["programs"].items():
        if frozen[frozen.program == name].sort_values("rank").gene.tolist() != genes:
            raise ValueError("Frozen signature order drift: " + name)
    environment = os.environ.copy()
    environment["IVDD_DISCOVERY_ROOT"] = str(discovery)
    environment["IVDD_DISCOVERY_REBUILD_ROOT"] = str(rebuilt)
    scripts = Path(__file__).resolve().parent / "human"
    python = config.get("executables", {}).get("discovery_python", sys.executable)
    rscript = config.get("executables", {}).get("rscript", "Rscript")
    records = []
    def run(name, command):
        start = time.time()
        with (discovery / "logs" / (name + ".log")).open("w") as log:
            result = subprocess.run(command, env=environment, stdout=log, stderr=subprocess.STDOUT)
        records.append({"step": name, "command": command, "exit_code": result.returncode, "seconds": time.time()-start})
        if result.returncode:
            raise RuntimeError(f"{name} failed; inspect its preserved log")
    try:
        if args.step in ["raw-qc", "all"]:
            run("raw_qc", [python, str(scripts / "discovery_qc.py"), "--study", "GSE230808"])
        if args.step in ["identity", "all"]:
            run("identity", [python, str(scripts / "discovery_identity.py"), "discovery"])
            doublets = pd.read_csv(rebuilt / "tables/discovery_doublet_summary.tsv", sep="\t")
            if not doublets.status.eq("completed").all():
                raise RuntimeError("Historical code recorded a Scrublet failure; do not call the rebuilt identity stage complete")
        if args.step in ["ucell", "all"]:
            run("ucell", [rscript, "--vanilla", str(scripts / "discovery_ucell.R"), str(rebuilt),
                           str(resource(config, "programs", "resources/frozen_k5/FROZEN_K5_PROGRAMS.tsv"))])
        if args.step == "cnmf":
            if not args.execute_historical_cnmf:
                raise RuntimeError("cNMF recomputation requires --execute-historical-cnmf; normal V3 reuses the frozen programs")
            run("cnmf_historical_candidates", [python, str(scripts / "discovery_cnmf.py")])
        elif args.step in ["audit-frozen", "all"]:
            records.append({"step": "cnmf", "status": "SKIPPED_FROZEN_K5_REUSED", "refit_performed": False})
    finally:
        proof = {"frozen_lock_sha256": hashlib.sha256(lock_file.read_bytes()).hexdigest(),
                 "step_requested": args.step, "executions": records, "private_frozen_model_artifacts": private_models,
                 "historical_fit_seed": 20261012, "rebuild_identity_seed": 20261016,
                 "public_code_scope": "Original own-code preprocessing and explicit cNMF invocation retained; frozen main K is not reselected."}
        (root / "00_audit/DISCOVERY_EXECUTION_RECORD.json").write_text(json.dumps(proof, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
