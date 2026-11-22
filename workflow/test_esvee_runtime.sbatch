#!/usr/bin/env bash
#SBATCH --job-name=esvee-jvm-smoke
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
# No sample data: verify the existing container accepts the exact JVM exclusion.
set -eo pipefail
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
module load apptainer/1.3.5
apptainer exec --cleanenv --env-file "$2" "$1" esvee -version 2>&1 | tee "$3"
grep -F 'CompileCommand: exclude htsjdk/samtools/AbstractBAMFileIndex.query' "$3"
grep -F 'Esvee version 2.0.1' "$3"
