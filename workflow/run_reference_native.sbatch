#!/usr/bin/env bash
# Arguments: frozen package, run ID. Allocation supplied by sbatch.
set -euo pipefail
package=$(readlink -f "$1")
run_id=$2
run_dir=$(dirname "$package")
mkdir -p "$run_dir/state"
python3 "$package/workflow/build_input_resources.py" \
    --samplesheet "$package/samplesheet.csv" \
    --out "$run_dir/state/input_resources.config" \
    --evidence "$run_dir/state/input_resources.json"
export WIGITS_WORKFLOW_DIR="$package/workflow"
export ONCO_SOURCE=nf-core/oncoanalyser ONCO_REV=3.0.0 ONCO_GENOME=GRCh38_hmf
bash "$package/workflow/run_oncoanalyser_frontenac.sh" \
    "$run_id" "$package/samplesheet.csv" \
    -c "$run_dir/state/input_resources.config" -c "$package/runtime.config"
