#!/usr/bin/env bash
# Resources are issued from measured BAM size in the execution package.
#SBATCH --job-name=full-bam-to-fastq
set -eo pipefail
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
module load StdEnv/2023 samtools/1.22.1
set -u
python3 "$1" --bam "$2" --sample "$3" --outdir "$4" \
    --threads "${SLURM_CPUS_PER_TASK:?Slurm allocation required}"
