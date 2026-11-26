#!/usr/bin/env bash
# Run explicitly selected post-native extensions; reusable for completed runs.
# Source/site configuration stays in the frozen native package. No data copying.
set -euo pipefail
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
NXF_SYNTAX_PARSER=v1 nextflow run "$HERE/optional_extensions.nf" -entry OPTIONAL_EXTENSIONS \
  -c "$site_config" -c "$HERE/conf/optional_extensions.config" --plan "$RUN_DIR/state/plan.json" --outdir "$RUN_DIR/results" \
  -work-dir "$RUN_DIR/work" -with-trace "$RUN_DIR/logs/trace.tsv" -resume "$@" || status=$?
python3 "$HERE/optional_extensions.py" status --plan "$RUN_DIR/state/plan.json" \
  --results "$RUN_DIR/results" --trace "$RUN_DIR/logs/trace.tsv" --out "$RUN_DIR/state/status.json"
# Missing required inputs for a selected eligible extension cannot be success.
python3 - "$RUN_DIR/state/status.json" <<'PY' || status=1
import json, sys
record=json.load(open(sys.argv[1]))
bad={"FAILED", "INCOMPLETE", "MISSING_INPUTS"}
sys.exit(any(c[f]["status"] in bad for c in record["cases"] for f in ("blood_ch", "rna_neoantigens")))
PY
exit "$status"
