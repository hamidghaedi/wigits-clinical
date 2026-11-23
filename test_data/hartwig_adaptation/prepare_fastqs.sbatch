#!/usr/bin/env bash
#SBATCH --job-name=chr17-19-fastqs
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
set -eo pipefail
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
module load StdEnv/2023 samtools/1.22.1
set -u
python3 "$1" --bam "$2" --sample "$3" --outdir "$4" \
    --threads "${SLURM_CPUS_PER_TASK:-8}"
