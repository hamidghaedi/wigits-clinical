#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=lab_run_qc
set -eo pipefail
HERE="${WIGITS_WORKFLOW_DIR:?Set WIGITS_WORKFLOW_DIR to the frozen release workflow directory}"
source "$HERE/conf/frontenac_paths.env"
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
module load "$NEXTFLOW_MODULE" "$APPTAINER_MODULE"
set -u
MANIFEST="$(readlink -f "${1:?manifest.json}")"
TASK_ROOT="$(readlink -m "${2:?execution directory}")"
shift 2
mkdir -p "$TASK_ROOT/logs"
cd "$TASK_ROOT"
RELEASE_COMMIT="$(git -C "$HERE" rev-parse HEAD)"
export NXF_APPTAINER_CACHEDIR="$CONTAINER_DIR" NXF_SINGULARITY_CACHEDIR="$CONTAINER_DIR"
export APPTAINER_CACHEDIR="$CONTAINER_DIR" NXF_ANSI_LOG=false
export NXF_SYNTAX_PARSER=v1
nextflow run "$HERE/laboratory_qc.nf" -entry LAB_RUN_QC \
  -c "$HERE/conf/frontenac.config" --manifest "$MANIFEST" --outdir "$TASK_ROOT/results" \
  --picard_image "$CONTAINER_DIR/picard-3.1.0.sif" \
  --fastp_image "$CONTAINER_DIR/fastp-1.0.1--heae3180_0.sif" \
  --multiqc_image "$CONTAINER_DIR/multiqc-1.32--pyhdfd78af_0.sif" \
  --release_tag lab-qc/v0.1.0 --release_commit "$RELEASE_COMMIT" \
  -work-dir "$TASK_ROOT/work" -with-trace "$TASK_ROOT/execution_trace.txt" "$@"
