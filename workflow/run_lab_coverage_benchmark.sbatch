#!/usr/bin/env bash
#SBATCH --job-name=lab-coverage-compare
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=1-00:00:00
set -eo pipefail
HERE="${WIGITS_WORKFLOW_DIR:?Frozen workflow source}"
MANIFEST="$(readlink -f "${1:?benchmark manifest}")"
RUN="$(readlink -m "${2:?new benchmark execution}")"
source "$HERE/conf/frontenac_paths.env"
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
module load "$NEXTFLOW_MODULE" "$APPTAINER_MODULE"
set -u
mkdir -p "$RUN"
cd "$RUN"
export NXF_ANSI_LOG=false NXF_SYNTAX_PARSER=v1
export NXF_APPTAINER_CACHEDIR="$CONTAINER_DIR" NXF_SINGULARITY_CACHEDIR="$CONTAINER_DIR"
nextflow run "$HERE/benchmark_lab_coverage.nf" -entry COVERAGE_BENCHMARK \
  -c "$HERE/conf/frontenac.config" --manifest "$MANIFEST" --outdir "$RUN/results" \
  --parabricks_image "$CONTAINER_DIR/clara-parabricks-4.6.0-1.sif" \
  -work-dir "$RUN/work" -with-trace "$RUN/execution_trace.txt" -resume
