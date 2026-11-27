#!/usr/bin/env bash
# Reuse a completed native case and export; run failed neo prediction only.
# Submit this script from the NEW pinned release, passing the original run root.
#SBATCH --job-name=optional-neo-recovery
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
set -euo pipefail
unset SBATCH_DEPENDENCY
HERE="${WIGITS_RECOVERY_WORKFLOW:?Export the new pinned workflow directory with sbatch}"
RUN="$(readlink -f "${1:?original native run root required}")"
RECOVERY="$RUN/extensions/recovery/native-intake-v1.1.5"
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
export WIGITS_EXTENSIONS_ENTRY=PREDICT_EXISTING_NEO
export WIGITS_EXTENSIONS_PRIOR_RESULTS="$RUN/extensions/results"
mkdir -p "$RECOVERY/state"
printf '%s\n' "$HERE" > "$RECOVERY/state/recovery_workflow.txt"
bash "$HERE/run_optional_extensions.sh" "$RUN/package/optional_extensions.json" "$RECOVERY" \
    --case_id "$case_id" --existing_export "$RUN/extensions/results/$case_id/neo_export"
