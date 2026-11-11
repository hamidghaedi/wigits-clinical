#!/usr/bin/env bash
# Arguments: frozen script, config JSON, annotated pilot VCF, NEW output directory.
set -eo pipefail
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
module load StdEnv/2023 apptainer/1.3.5
exec python3 "$1" --config "$2" --adapted "$3" --out "$4"
