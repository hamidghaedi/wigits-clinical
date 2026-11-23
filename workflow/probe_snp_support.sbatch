#!/usr/bin/env bash
#SBATCH --job-name=targeted-snp-support
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
set -eo pipefail
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
module load apptainer/1.3.5
python_bin=$1
script=$2
shift 2
env -u PYTHONPATH "$python_bin" -S "$script" "$@"
