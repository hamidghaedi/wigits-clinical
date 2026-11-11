#!/usr/bin/env bash
# Arguments: explicit config JSON, NEW output directory, frozen source directory.
# Small annotation/test job, not a caller or native workflow submission.
set -eo pipefail
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
module load StdEnv/2023 apptainer/1.3.5
exec python3 "$3/run_deepsomatic_adapter.py" --config "$1" --out "$2" --source "$3"
