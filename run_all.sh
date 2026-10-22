#!/usr/bin/env bash
set -euo pipefail
SCRIPT_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
exec "${PYTHON_SCIENCE:-python}" "$SCRIPT_ROOT/scripts/run_project.py" "$@"
