#!/usr/bin/env python3
"""Figure 1-2 adapter for a project-wide figure driver."""
import argparse
from pathlib import Path
import subprocess
import sys
from human_config import load_config, data_path

parser = argparse.ArgumentParser()
parser.add_argument("--config", required=True)
parser.add_argument("--project-root")
args = parser.parse_args()
config = load_config(args.config, args.project_root)
command = [config.get("executables", {}).get("human_python", sys.executable),
           str(Path(__file__).resolve().parent / "human/plot_figures12.py"), "--work", str(config["_root"])]
if config.get("discovery", {}).get("af_donor_counts"):
    command += ["--af-donor-counts", str(data_path(config, config["discovery"]["af_donor_counts"]))]
raise SystemExit(subprocess.call(command))
