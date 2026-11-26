#!/usr/bin/env bash
#SBATCH --job-name=lab-qc-event-test
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:05:00
set -eo pipefail
HERE="${WIGITS_WORKFLOW_DIR:?Frozen workflow directory}"
NATIVE="${1:?Prepared isolated native source}"
TEST_ROOT="$(readlink -m "${2:?New synthetic test directory}")"
source "$HERE/conf/frontenac_paths.env"
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
module load "$NEXTFLOW_MODULE"
set -u
mkdir -p "$TEST_ROOT/lib"
cd "$TEST_ROOT"
export NXF_OPTS='-Xms256m -Xmx2g' NXF_ANSI_LOG=false
nextflow -log native_help.log run "$NATIVE" \
  -c "$ONCO_REF_ROOT/$ONCO_REFERENCE_VERSION/reference_data.config" \
  --genome GRCh38_hmf --mode wgts --help > native_help.txt
cp "$HERE/../tests/fixtures/lab_qc_events/notification.nf" test.nf
cp "$HERE/overlays/lab-qc-events/lib/LabQcEvents.groovy" lib/
python3 - "$HERE" "$TEST_ROOT" <<'PY'
import json,sys
from pathlib import Path
here,root=map(Path,sys.argv[1:])
(root/'event_config.json').write_text(json.dumps(dict(python=sys.executable,
    dispatcher=str(here.parent/'tests/fixtures/lab_qc_events/capture_event.py'),
    receipt=str(root/'notification.json')))+'\n')
PY
export WIGITS_LAB_QC_CONFIG="$TEST_ROOT/event_config.json"
nextflow -log canary.log run test.nf -process.executor local -process.cpus 1 \
  -process.memory '512 MB' -with-trace execution_trace.txt
python3 - "$TEST_ROOT" <<'PY'
from pathlib import Path
from datetime import datetime
import json,sys
root=Path(sys.argv[1]);event=json.loads((root/'notification.json').read_text())
completed=[int(p.read_text().strip()) for p in (root/'work').glob('*/*/caller_completed.txt')]
assert len(completed)==1 and event['completed_output_present']
assert datetime.fromisoformat(event['received_at']).timestamp()<completed[0]
(root/'verification.json').write_text(json.dumps(dict(event=event,caller_completed_epoch=completed[0],
    notified_before_caller_finished=True,native_help_exit=0),indent=2)+'\n')
print('BAM event observed before caller completion; no periodic scan or analytical execution')
PY
