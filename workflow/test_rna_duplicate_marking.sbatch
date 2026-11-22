#!/usr/bin/env bash
#SBATCH --job-name=rna-duplicate-policy-test
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
set -eo pipefail
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
module load apptainer/1.3.5
unset PYTHONPATH
exec python3 "$1" --container "$2" --out "$3"
