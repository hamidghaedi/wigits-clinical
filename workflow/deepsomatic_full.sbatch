#!/usr/bin/env bash
# Resources/account/partition/dependency/output are explicit sbatch arguments.
# Arguments: frozen Python, config JSON, output directory, stage.
set -eo pipefail
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
module load StdEnv/2023 apptainer/1.3.5
exec python3 "$1" --config "$2" --out "$3" --stage "$4"
