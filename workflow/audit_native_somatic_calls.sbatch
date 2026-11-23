#!/usr/bin/env bash
#SBATCH --job-name=native-call-qc
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
# Read-only VCF summaries; no BAM/FASTQ scan, caller or report modification.
set -euo pipefail
python_bin=$1
script=$2
manifest=$3
output=$4
env -u PYTHONPATH "$python_bin" -S "$script" --manifest "$manifest" \
    --index "$SLURM_ARRAY_TASK_ID" --out "$output"
