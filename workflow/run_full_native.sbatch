#!/usr/bin/env bash
#SBATCH --job-name=onco-full-native-head
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=13-00:00:00
set -euo pipefail
PACKAGE="$(readlink -f "${1:?execution package required}")"
RUN_ID="${2:?run ID required}"
RUN_DIR="$(dirname "$PACKAGE")"

# Original FASTQs are already staged. Keep one native row per lane and size
# the libraries from the existing inputs; no BAM recovery or recovery receipts.
mkdir -p "$RUN_DIR/state"
python3 "$PACKAGE/workflow/build_input_resources.py" \
    --samplesheet "$PACKAGE/samplesheet.csv" \
    --out "$RUN_DIR/state/input_resources.config" \
    --evidence "$RUN_DIR/state/input_resources.json"

export ONCO_SOURCE=nf-core/oncoanalyser ONCO_REV=3.0.0 ONCO_GENOME=GRCh38_hmf
export WIGITS_WORKFLOW_DIR="$PACKAGE/workflow"
bash "$PACKAGE/workflow/run_oncoanalyser_frontenac.sh" \
    "$RUN_ID" "$PACKAGE/samplesheet.csv" -c "$RUN_DIR/state/input_resources.config"
