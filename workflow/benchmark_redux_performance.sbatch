#!/usr/bin/env bash
#SBATCH --job-name=redux-performance-pilot
#SBATCH --cpus-per-task=64
#SBATCH --mem=192G
#SBATCH --time=02:00:00
# Submit with --account and an external --output; arguments: script native-work out image.
set -eo pipefail
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
module load StdEnv/2023 apptainer/1.3.5
exec python3 "$1" --native-work "$2" --out "$3" --image "$4"
