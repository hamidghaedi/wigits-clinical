#!/usr/bin/env bash
# Internal head for explicit optional entries selected by the common launcher.
#SBATCH --job-name=optional-neo-recovery
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
set -euo pipefail
unset SBATCH_DEPENDENCY
HERE="${WIGITS_RECOVERY_WORKFLOW:?Export the new pinned workflow directory with sbatch}"
RUN="$(readlink -f "${1:?original native run root required}")"
ENTRY="${2:?neo-prepare, neo-predict or neo-review required}"
PREDICTIONS="${3:-}"
EXPORT="${4:-}"
case "$ENTRY" in neo-prepare|neo-predict|neo-review) ;; *) echo 'Unsupported optional entry' >&2; exit 2 ;; esac
RECOVERY="$RUN/extensions/recovery/native-intake-v1.3.0/$ENTRY"
source "$RUN/package/site.env"
unset SBATCH_DEPENDENCY
set +u
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
module load "$NEXTFLOW_MODULE" "$APPTAINER_MODULE"
set -u
case_id="$(python3 - "$RUN" <<'PY'
import json, sys
from pathlib import Path
root = Path(sys.argv[1])
native = json.loads((root/'state/native_workflow_completed.json').read_text())
if native.get('status') != 'COMPLETED':
    raise SystemExit('Native completion is required; this entry never reruns native work')
config = json.loads((root/'package/optional_extensions.json').read_text())
if len(config['cases']) != 1:
    raise SystemExit('Recovery expects one native case per run')
print(config['cases'][0]['case_id'])
PY
)"
extra=()
case "$ENTRY" in
    neo-prepare) export WIGITS_EXTENSIONS_ENTRY=PREPARE_EXISTING_NEO ;;
    neo-predict) export WIGITS_EXTENSIONS_ENTRY=PREDICT_EXISTING_NEO ;;
    neo-review)
        export WIGITS_EXTENSIONS_ENTRY=REVIEW_EXISTING_NEO
        extra=(--existing_predictions "${PREDICTIONS:?Completed prediction directory required}") ;;
esac
if [[ "$ENTRY" != neo-prepare ]]; then
    extra+=(--existing_export "${EXPORT:-$RUN/extensions/results/$case_id/neo_export}")
fi
export WIGITS_EXTENSIONS_PRIOR_RESULTS="$RUN/extensions/results"
mkdir -p "$RECOVERY/state"
printf '%s\n' "$HERE" > "$RECOVERY/state/recovery_workflow.txt"
bash "$HERE/run_optional_extensions.sh" "$RUN/package/optional_extensions.json" "$RECOVERY" \
    --case_id "$case_id" "${extra[@]}"
