#!/usr/bin/env bash
# Run explicitly selected post-native extensions; reusable for completed runs.
# Source/site configuration stays in the frozen native package. No data copying.
set -euo pipefail
# Standalone or post-native execution must not inherit head dependencies.
unset SBATCH_DEPENDENCY
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="$(readlink -f "${1:?optional_extensions.json required}")"
RUN_DIR="$(readlink -m "${2:?extension execution root required}")"
shift 2
mkdir -p "$RUN_DIR/state" "$RUN_DIR/logs"
cd "$RUN_DIR"
python3 "$HERE/optional_extensions.py" plan --config "$CONFIG" --out "$RUN_DIR/state/plan.json"
# Nextflow and Apptainer must already be loaded by the caller. A site config is
# explicit for standalone execution, and inherited from the native launcher.
site_config="${WIGITS_EXTENSIONS_SITE_CONFIG:-$HERE/conf/frontenac.config}"
status=0
NXF_SYNTAX_PARSER=v1 nextflow run "$HERE/optional_extensions.nf" -entry "${WIGITS_EXTENSIONS_ENTRY:-OPTIONAL_EXTENSIONS}" \
  -c "$site_config" -c "$HERE/conf/optional_extensions.config" --plan "$RUN_DIR/state/plan.json" --outdir "$RUN_DIR/results" \
  -work-dir "$RUN_DIR/work" -with-trace "$RUN_DIR/logs/trace.tsv" -resume "$@" || status=$?
# Export-only is an entry result, not a claim that prediction is completed.
if [[ "${WIGITS_EXTENSIONS_ENTRY:-}" == PREPARE_EXISTING_NEO ]]; then
    python3 - "$RUN_DIR" <<'PY_ENTRY' || status=1
import json, sys
from pathlib import Path
root = Path(sys.argv[1])
products = list((root/'results').glob('*/neo_export/status.json'))
records = [json.loads(p.read_text()) for p in products]
ok = len(records) == 1 and records[0]['status'] in {'PREPARED','NO_CANDIDATES','NOT_ASSESSED'}
(root/'state/entry_status.json').write_text(json.dumps(dict(entry='neo-prepare',
    status='COMPLETED' if ok else 'INCOMPLETE', prediction='NOT_RUN', products=records), indent=2)+'\n')
sys.exit(not ok)
PY_ENTRY
    exit "$status"
fi
prior=()
[[ -z "${WIGITS_EXTENSIONS_PRIOR_RESULTS:-}" ]] || prior=(--prior-results "$WIGITS_EXTENSIONS_PRIOR_RESULTS")
python3 "$HERE/optional_extensions.py" status --plan "$RUN_DIR/state/plan.json" \
  --results "$RUN_DIR/results" --trace "$RUN_DIR/logs/trace.tsv" --out "$RUN_DIR/state/status.json" "${prior[@]}"
# Missing required inputs for a selected eligible extension cannot be success.
python3 - "$RUN_DIR/state/status.json" <<'PY' || status=1
import json, sys
record=json.load(open(sys.argv[1]))
bad={"FAILED", "INCOMPLETE", "MISSING_INPUTS"}
sys.exit(any(c[f]["status"] in bad for c in record["cases"] for f in ("blood_ch", "rna_neoantigens")))
PY
exit "$status"
