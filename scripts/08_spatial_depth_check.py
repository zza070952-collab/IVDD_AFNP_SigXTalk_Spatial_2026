#!/usr/bin/env python3
"""Run the single frozen depth check after raw or cached spatial preparation."""
from __future__ import annotations
import argparse
import importlib.util
from pathlib import Path
import subprocess
import sys

def main():
    here = Path(__file__).resolve().parent
    spec = importlib.util.spec_from_file_location('spatial_entry', here/'07_prepare_spatial.py')
    entry = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(entry)
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config', required=True, type=Path)
    ap.add_argument('--project-root', required=True, type=Path)
    ap.add_argument('--cached-spatial', type=Path, help='Read-only original baseline override')
    ap.add_argument('--from-raw', action='store_true', help='Use baseline_cache created by 07 raw mode')
    a = ap.parse_args()
    cfg = entry.load_config(a.config, a.project_root)
    root = a.cached_spatial.resolve() if a.cached_spatial else cfg['baseline_cache' if a.from_raw else 'cached_spatial']
    entry.validate_cached(root)
    result = subprocess.run([sys.executable, str(here/'spatial_depth_check.py'),
                             '--spatial-input', str(root), '--roi-json', str(cfg['roi_config']),
                             '--config', str(cfg['depth_config']), '--out', str(cfg['depth_output'])], check=False)
    raise SystemExit(result.returncode)

if __name__ == '__main__':
    main()
