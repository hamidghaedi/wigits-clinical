#!/usr/bin/env bash
# Resource requests are issued by the input-sized intake dispatcher.
#SBATCH --job-name=onco-bam-intake
set -eo pipefail
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
module load StdEnv/2023 samtools/1.22.1 apptainer/1.3.5
set -u
exec python3 "$1" prepare --bam "$2" --sample "$3" --outdir "$4" \
  --gatk-command "$5" --threads "${SLURM_CPUS_PER_TASK:?}"
