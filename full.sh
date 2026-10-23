#!/usr/bin/env bash
# Compatibility entry point for the original full.sh.
set -euo pipefail
base_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
exec python "$base_dir/methmap.py" build "$@"
