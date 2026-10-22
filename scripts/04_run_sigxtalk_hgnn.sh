#!/usr/bin/env bash
# Own orchestration; upstream SigXTalk is installed from its pinned repository.
set -euo pipefail
if [[ $# -lt 1 ]]; then
  echo 'Usage: bash 04_run_sigxtalk_hgnn.sh configs/sigxtalk.json [smoke|formal|all] [prepare|hgnn|prs|all]' >&2
  exit 2
fi
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
exec "${ORCHESTRATOR_PYTHON:-python3}" "$script_dir/run_sigxtalk_pipeline.py" \
  --config "$1" --mode "${2:-all}" --stage "${3:-all}"
