#!/usr/bin/env bash
#SBATCH --account=def-hpcg6049_cpu
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --job-name=oncoanalyser_test_head
#
# run_oncoanalyser_test_frontenac.sh — smoke-test nf-core/oncoanalyser on Frontenac.
#
# The Nextflow HEAD runs inside this Slurm job (not on the login node); it submits every
# task as its own Slurm job via conf/frontenac.config. Safe to re-run: -resume is on.
#
#   sbatch workflow/run_oncoanalyser_test_frontenac.sh
#
# Override with env vars: ONCO_REV (default 3.0.0), RUN_ROOT, SITE_CONFIG.
set -eo pipefail   # -u is enabled AFTER the cvmfs profile is sourced: that script reads $SKIP_CC_CVMFS unset

ONCO_REV="${ONCO_REV:-3.0.0}"
RUN_ROOT="${RUN_ROOT:-/global/scratch/hpc6049/oncoanalyser/test_${ONCO_REV}}"
SITE_CONFIG="${SITE_CONFIG:-/global/project/hpcg6049/somatic_pipeline/wigits-clinical/workflow/conf/frontenac.config}"

source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
module load nextflow/26.04.4 apptainer/1.3.5          # oncoanalyser 3.0.0 needs Nextflow >= 25.10.4
set -u

export NXF_SINGULARITY_CACHEDIR=/global/project/hpcg6049/somatic_pipeline/containers/NXF_SINGULARITY_CACHEDIR
export NXF_APPTAINER_CACHEDIR="$NXF_SINGULARITY_CACHEDIR"
export APPTAINER_CACHEDIR=/global/scratch/hpc6049/tmp
# Image pulls run mksquashfs on the HEAD node. With build-temp on GPFS a 2-CPU head needed >35 min for
# the GRIDSS image (2026-10-14); node-local $SLURM_TMPDIR is much faster. First run of a new
# pipeline version pulls ~20 WiGiTS images — give the head 4-8 CPUs for that run.
export APPTAINER_TMPDIR="${SLURM_TMPDIR:-/global/scratch/hpc6049/tmp}"
export NXF_OPTS='-Xms1g -Xmx4g'                        # login shell default is -Xmx2g with 1 CPU; the head gets 8 GB here
export NXF_ANSI_LOG=false
mkdir -p "$APPTAINER_CACHEDIR" "$RUN_ROOT"
cd "$RUN_ROOT"

echo "node=$(hostname) rev=$ONCO_REV run_root=$RUN_ROOT"
echo "nextflow=$(nextflow -version 2>/dev/null | grep -i version) apptainer=$(apptainer --version)"
echo "site_config=$SITE_CONFIG"; sha256sum "$SITE_CONFIG"

nextflow run nf-core/oncoanalyser \
    -r "$ONCO_REV" \
    -profile test,apptainer \
    -c "$SITE_CONFIG" \
    -work-dir "$RUN_ROOT/work" \
    --outdir "$RUN_ROOT/results" \
    -resume

status=$?
echo "nextflow exit=$status"
# Harvest the head log + pipeline_info alongside the sbatch log for whoever launched us.
if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -d "$SLURM_SUBMIT_DIR" ]; then
  cp -f .nextflow.log "$SLURM_SUBMIT_DIR/nextflow_head.log" 2>/dev/null || true
  cp -rf results/pipeline_info "$SLURM_SUBMIT_DIR/pipeline_info" 2>/dev/null || true
fi
exit $status
