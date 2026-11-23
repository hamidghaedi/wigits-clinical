#!/usr/bin/env bash
# Arguments: frozen script, config JSON, annotated pilot VCF, NEW output directory.
set -eo pipefail
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
module load StdEnv/2023 apptainer/1.3.5
script=$1; config=$2; adapted=$3; out=$4
shift 4
exec python3 "$script" --config "$config" --adapted "$adapted" --out "$out" "$@"
