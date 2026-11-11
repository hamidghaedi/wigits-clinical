#!/usr/bin/env bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=13-00:00:00
# Replay an explicitly recorded native command with one final recovery config.
# Does not prepare BAMs, regenerate resources or replace sample sheets.
set -eo pipefail
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
module load nextflow/26.04.4 apptainer/1.3.5
package=$1
cd "$(cat "$package/run_directory.txt")"
export NXF_SINGULARITY_CACHEDIR=/global/project/hpcg6049/somatic_pipeline/containers/NXF_SINGULARITY_CACHEDIR
export NXF_APPTAINER_CACHEDIR="$NXF_SINGULARITY_CACHEDIR"
export APPTAINER_CACHEDIR=/global/scratch/hpc6049/tmp
export APPTAINER_TMPDIR="${SLURM_TMPDIR:-/global/scratch/hpc6049/tmp}"
export NXF_OPTS='-Xms1g -Xmx4g'
export NXF_ANSI_LOG=false
exec bash "$package/resume_command.sh"
