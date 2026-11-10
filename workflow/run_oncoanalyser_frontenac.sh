#!/usr/bin/env bash
#SBATCH --account=def-hpcg6049_cpu
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --job-name=oncoanalyser_head
#
# run_oncoanalyser_frontenac.sh — run nf-core/oncoanalyser (as designed) on Frontenac for one run_id.
#
#   sbatch workflow/run_oncoanalyser_frontenac.sh <run_id> <samplesheet.csv> [extra nextflow args...]
#
# The Nextflow HEAD runs inside this Slurm job; every pipeline task is its own Slurm job
# (conf/frontenac.config). Paths come from conf/frontenac_paths.env (docs/directory_layout.md).
# Run layout: $RUN_ROOT/<run_id>/{samplesheet.csv,launch.sh,work,results,logs}. -resume is on.
set -eo pipefail   # -u only after the cvmfs profile is sourced (it reads unset SKIP_CC_CVMFS)

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ ! -f "$HERE/conf/frontenac_paths.env" ]]; then
    # sbatch copies this script to its spool; configuration stays in the checkout.
    HERE="${WIGITS_WORKFLOW_DIR:-${SLURM_SUBMIT_DIR:?submit from the repository root}/workflow}"
fi
[[ -f "$HERE/conf/frontenac_paths.env" ]] || {
    echo "workflow configuration not found in $HERE; set WIGITS_WORKFLOW_DIR" >&2
    exit 2
}
source "$HERE/conf/frontenac_paths.env"
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
module load "$NEXTFLOW_MODULE" "$APPTAINER_MODULE"
set -u

RUN_ID="${1:?usage: $0 <run_id> <samplesheet.csv> [extra args]}"
SAMPLESHEET="$(readlink -f "${2:?usage: $0 <run_id> <samplesheet.csv> [extra args]}")"
shift 2
ONCO_REV="${ONCO_REV:-$ONCO_REFERENCE_VERSION}"
ONCO_SOURCE="${ONCO_SOURCE:-nf-core/oncoanalyser}"
ONCO_GENOME="${ONCO_GENOME:-GRCh38_hmf}"
REF_CONFIG="${REF_CONFIG:-$ONCO_REF_ROOT/$ONCO_REFERENCE_VERSION/reference_data.config}"
RUN_DIR="$RUN_ROOT/$RUN_ID"

export NXF_SINGULARITY_CACHEDIR="$CONTAINER_DIR"
export NXF_APPTAINER_CACHEDIR="$CONTAINER_DIR"
export APPTAINER_CACHEDIR=/global/scratch/hpc6049/tmp
export APPTAINER_TMPDIR="${SLURM_TMPDIR:-/global/scratch/hpc6049/tmp}"
export NXF_OPTS='-Xms1g -Xmx4g'
export NXF_ANSI_LOG=false

mkdir -p "$RUN_DIR/logs" "$APPTAINER_CACHEDIR"
cd "$RUN_DIR"
cp -f "$SAMPLESHEET" samplesheet.csv
[ -f "$REF_CONFIG" ] || { echo "reference config not found: $REF_CONFIG" >&2; exit 2; }

if [[ "$ONCO_SOURCE" == nf-core/oncoanalyser ]]; then
    CMD=(nextflow run "$ONCO_SOURCE" -r "$ONCO_REV")
else
    [[ -f "$ONCO_SOURCE/main.nf" ]] || { echo "pipeline source not found: $ONCO_SOURCE" >&2; exit 2; }
    CMD=(nextflow run "$ONCO_SOURCE")
fi
CMD+=(-profile apptainer
     -c "$HERE/conf/frontenac.config"
     -c "$HERE/conf/frontenac_limits.config"
     -c "$REF_CONFIG"
     --mode wgts --genome "$ONCO_GENOME"
     --input "$RUN_DIR/samplesheet.csv"
     --outdir "$RUN_DIR/results"
     -work-dir "$RUN_DIR/work"
     -resume "$@")

{ echo "# $(date -u +%FT%TZ) Slurm job ${SLURM_JOB_ID:-none} on $(hostname)";
  echo "# nextflow $(nextflow -version 2>/dev/null | grep -i version) | $(apptainer --version)";
  echo "# samplesheet sha256: $(sha256sum samplesheet.csv | cut -d' ' -f1)";
  echo "# reference config: $REF_CONFIG";
  echo "# source: $ONCO_SOURCE";
  if [[ -d "$ONCO_SOURCE/.git" || -f "$ONCO_SOURCE/.git" ]]; then
      echo "# source commit: $(git -C "$ONCO_SOURCE" rev-parse HEAD)";
  fi;
  printf '%q ' "${CMD[@]}"; echo; } | tee -a launch.sh

"${CMD[@]}"
status=$?
echo "nextflow exit=$status"
cp -f .nextflow.log "logs/nextflow_head_$(date -u +%Y%m%dT%H%M%SZ).log" 2>/dev/null || true
exit $status
