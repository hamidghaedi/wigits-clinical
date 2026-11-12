#!/usr/bin/env bash
# Arguments are passed directly to the frozen test_evidence_consumers.py.
set -eo pipefail
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
module load StdEnv/2023 apptainer/1.3.5
exec python3 "$@"
