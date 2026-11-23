#!/usr/bin/env bash
#SBATCH --job-name=google-deepsomatic-seam-test
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --gres=gpu:1
#SBATCH --time=04:00:00
# Supply account/output on submission. Arguments: frozen Python script, JSON, new output root.
# Small bounded test only; full-WGS CPU example preparation must be sized separately.
set -eo pipefail
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
module load StdEnv/2023 apptainer/1.3.5
exec python3 "$1" --config "$2" --out "$3"
