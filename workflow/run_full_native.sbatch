#!/usr/bin/env bash
#SBATCH --job-name=onco-full-native-head
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=4-00:00:00
set -euo pipefail
PACKAGE="$(readlink -f "${1:?execution package required}")"
RUN_ID="${2:?run ID required}"
shift 2
RUN_DIR="$(dirname "$PACKAGE")"

# Original FASTQs are already staged. Keep one native row per lane and size
# the libraries from the existing inputs; no BAM recovery or recovery receipts.
mkdir -p "$RUN_DIR/state"
resource_resume=()
if [[ -f "$PACKAGE/completed_alignment_cpus.json" ]]; then
    resource_resume=(--completed-alignment-cpus "$PACKAGE/completed_alignment_cpus.json")
fi
python3 "$PACKAGE/workflow/build_input_resources.py" \
    --samplesheet "$PACKAGE/samplesheet.csv" \
    --out "$RUN_DIR/state/input_resources.config" \
    --evidence "$RUN_DIR/state/input_resources.json" "${resource_resume[@]}"

if [[ -z "${WIGITS_SITE_ENV:-}" ]]; then
    export WIGITS_SITE_ENV="$PACKAGE/workflow/conf/frontenac_paths.env"
    [[ ! -f "$PACKAGE/site.env" ]] || export WIGITS_SITE_ENV="$PACKAGE/site.env"
fi
export WIGITS_WORKFLOW_DIR="$PACKAGE/workflow"
source "$WIGITS_SITE_ENV"
# A direct sbatch/site dependency belongs to this head, not its children.
unset SBATCH_DEPENDENCY
export ONCO_SOURCE="${ONCO_SOURCE:-$ONCO_SOURCE_DIR}" ONCO_REV=3.0.0 ONCO_GENOME=GRCh38_hmf
site_final=()
if [[ -n "${WIGITS_FINAL_CONFIG:-}" ]]; then
    site_final=(-c "$WIGITS_FINAL_CONFIG")
fi
bash "$PACKAGE/workflow/run_oncoanalyser_frontenac.sh" \
    "$RUN_ID" "$PACKAGE/samplesheet.csv" -c "$RUN_DIR/state/input_resources.config" \
    -c "$PACKAGE/workflow/conf/rna_duplicate_marking.config" "$@" "${site_final[@]}"

# The native receipt is independent of the selected downstream research work.
# Absent configuration preserves the historical native execution exactly.
if [[ -f "$PACKAGE/optional_extensions.json" ]]; then
    python3 - "$RUN_DIR/state/native_workflow_completed.json" "$RUN_ID" <<'PY'
import json, sys
from datetime import datetime, timezone
with open(sys.argv[1], 'w') as stream:
    json.dump(dict(run_id=sys.argv[2], status='COMPLETED', native_exit=0,
                   recorded_at=datetime.now(timezone.utc).isoformat()), stream, indent=2)
PY
    # Modules loaded in the native child shell do not persist in this shell.
    set +u
    source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
    module load "$NEXTFLOW_MODULE" "$APPTAINER_MODULE"
    set -u
    # A downstream failure has its own nonzero exit and status; it never erases
    # the native completion receipt or causes alignment to be repeated.
    bash "$PACKAGE/workflow/run_optional_extensions.sh" \
        "$PACKAGE/optional_extensions.json" "$RUN_DIR/extensions"
fi
