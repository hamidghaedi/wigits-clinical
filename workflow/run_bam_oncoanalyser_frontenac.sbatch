#!/usr/bin/env bash
#SBATCH --job-name=onco-native-bam-head
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=13-00:00:00
set -eo pipefail
HERE="${WIGITS_WORKFLOW_DIR:-${SLURM_SUBMIT_DIR:?submit from repository root}/workflow}"
source "$HERE/conf/frontenac_paths.env"
RUN_ID="${1:?run ID required}"
SHEET="$(readlink -f "${2:?BAM samplesheet required}")"
shift 2
mkdir -p "$RUN_ROOT/$RUN_ID/state"
python3 "$HERE/build_input_resources.py" --samplesheet "$SHEET" \
    --out "$RUN_ROOT/$RUN_ID/state/input_resources.config" \
    --evidence "$RUN_ROOT/$RUN_ID/state/input_resources.json"
export ONCO_SOURCE=nf-core/oncoanalyser ONCO_REV=3.0.0 ONCO_GENOME=GRCh38_d1_vd1
export WIGITS_WORKFLOW_DIR="$HERE"
exec bash "$HERE/run_oncoanalyser_frontenac.sh" "$RUN_ID" "$SHEET" \
    -c "$HERE/conf/d1_vd1.config" --genome_version 38 --genome_type no_alt --force_genome \
    -c "$RUN_ROOT/$RUN_ID/state/input_resources.config" "$@"
