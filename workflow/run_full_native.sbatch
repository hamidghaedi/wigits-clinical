#!/usr/bin/env bash
#SBATCH --job-name=onco-full-native-head
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=14-00:00:00
set -euo pipefail
PACKAGE="$(readlink -f "${1:?execution package required}")"
RUN_ID="${2:?run ID required}"
RUN_DIR="$(dirname "$PACKAGE")"

# Dependencies must finish before this script starts. Completion records account
# for every source primary read and all paired FASTQ names, counts and CRCs.
python3 - "$PACKAGE/samplesheet.csv" "$RUN_DIR" <<'PY'
import csv, json, subprocess, sys
from pathlib import Path
sheet, root = Path(sys.argv[1]), Path(sys.argv[2])
rows = list(csv.DictReader(sheet.open()))
assert len(rows) in (2, 3) and all(r['filetype'] == 'fastq' for r in rows)
command = [sys.executable, str(sheet.parent/'workflow/build_input_resources.py'),
           '--out', str(root/'state/input_resources.config'),
           '--evidence', str(root/'state/input_resources.json')]
for row in rows:
    completion_path = root/'inputs'/row['sample_id']/'completion.json'
    completion = json.loads(completion_path.read_text())
    assert completion['sample_id'] == row['sample_id']
    assert row['filepath'].split(';') == [completion['outputs'][x] for x in ['R1', 'R2']]
    role = 'rna' if row['sequence_type'] == 'rna' else row['sample_type']
    command += ['--library', role, str(completion_path)]
subprocess.run(command, check=True)
print('Full-input samplesheet matches every verified recovery record', flush=True)
PY

export ONCO_SOURCE=nf-core/oncoanalyser ONCO_REV=3.0.0
export WIGITS_WORKFLOW_DIR="$PACKAGE/workflow"
bash "$PACKAGE/workflow/run_oncoanalyser_frontenac.sh" \
    "$RUN_ID" "$PACKAGE/samplesheet.csv" -c "$RUN_DIR/state/input_resources.config"
