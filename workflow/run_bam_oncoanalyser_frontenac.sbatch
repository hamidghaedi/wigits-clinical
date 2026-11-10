#!/usr/bin/env bash
#SBATCH --job-name=onco-native-bam-head
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=13-00:00:00
set -eo pipefail
HERE="${WIGITS_WORKFLOW_DIR:-${SLURM_SUBMIT_DIR:?submit from repository root}/workflow}"
source "$HERE/conf/frontenac_paths.env"
source "$HERE/conf/bam_intake.env"
RUN_ID="${1:?run ID required}"
CONFIG="$(readlink -f "${2:?BAM YAML configuration required}")"
shift 2
mkdir -p "$RUN_ROOT/$RUN_ID/state"
PREPARED_SHEET="$RUN_ROOT/$RUN_ID/state/bam_samplesheet_v0.3.0.csv"
python3 "$HERE/bam_intake.py" batch --config "$CONFIG" \
    --outdir "$RUN_ROOT/$RUN_ID/inputs/prepared_bams_v0.3.0" --prepared-samplesheet "$PREPARED_SHEET" \
    --gatk-command "$INTAKE_GATK_COMMAND" \
    --account "$SLURM_ACCOUNT_CPU"
SHEET="$PREPARED_SHEET"
python3 "$HERE/build_input_resources.py" --samplesheet "$SHEET" \
    --out "$RUN_ROOT/$RUN_ID/state/input_resources.config" \
    --evidence "$RUN_ROOT/$RUN_ID/state/input_resources.json"
export ONCO_SOURCE=nf-core/oncoanalyser ONCO_REV=3.0.0 ONCO_GENOME=GRCh38_d1_vd1
export WIGITS_WORKFLOW_DIR="$HERE"
# Existing Nextflow cache may still be held while a failed attempt finishes
# reference preparation. Wait without cancelling or replacing that work.
if [[ -n "${ONCO_PRIOR_HEAD_JOB:-}" ]]; then
    while true; do
        prior_state="$(sacct -n -X -j "$ONCO_PRIOR_HEAD_JOB" --format=State -P)"
        case "$prior_state" in
            COMPLETED*|FAILED*|CANCELLED*|TIMEOUT*|OUT_OF_MEMORY*|NODE_FAIL*) break ;;
            '') echo 'Unable to establish prior controller terminal state' >&2; exit 2 ;;
            *) sleep 30 ;;
        esac
    done
fi
exec bash "$HERE/run_oncoanalyser_frontenac.sh" "$RUN_ID" "$SHEET" \
    -c "$HERE/conf/d1_vd1.config" --genome_version 38 --genome_type no_alt --force_genome \
    -c "$RUN_ROOT/$RUN_ID/state/input_resources.config" "$@"
