#!/usr/bin/env bash
# sbatch [cluster options] scripts/submit_slurm.sh build [MethMap options]
# Activate the methmap environment before submitting. No cluster-specific paths.
set -euo pipefail
# SLURM copies this script into its spool; BASH_SOURCE is no longer the repo path.
repo_dir="${METHMAP_HOME:-${SLURM_SUBMIT_DIR:-$PWD}}"
if [[ ! -f "$repo_dir/methmap.py" ]]; then
  echo "Submit from the MethMap repository root, or set METHMAP_HOME to that directory." >&2
  exit 1
fi
exec python "$repo_dir/methmap.py" "$@"
