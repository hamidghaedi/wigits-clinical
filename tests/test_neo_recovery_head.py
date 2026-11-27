"""Execute the reusable recovery head with synthetic receipts and mocked modules."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class NeoRecoveryHeadTests(unittest.TestCase):
    def run_head(self, native_status, entry='neo-predict'):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        base = Path(temporary.name)
        run, workflow = base/'original', base/'new_release/workflow'
        workflow.mkdir(parents=True)
        (run/'package').mkdir(parents=True)
        (run/'state').mkdir()
        native = run/'state/native_workflow_completed.json'
        native.write_text(json.dumps({'status':native_status}))
        config = run/'package/optional_extensions.json'
        config.write_text(json.dumps({'cases':[{'case_id':'SYNTHETIC'}]}))
        (run/'package/site.env').write_text('NEXTFLOW_MODULE=synthetic-nextflow\nAPPTAINER_MODULE=synthetic-apptainer\nexport SBATCH_DEPENDENCY=afterany:123\n')
        ch = run/'extensions/results/SYNTHETIC/blood_ch/status.json'
        ch.parent.mkdir(parents=True)
        ch.write_text('{"status":"COMPLETED"}\n')
        modules = base/'modules.sh'
        modules.write_text('''source() {
if [[ "$1" == /cvmfs/soft.computecanada.ca/config/profile/bash.sh ]]; then return 0; fi
builtin source "$@"
}
module() { :; }
''')
        recorder = base/'record.py'
        recorder.write_text('''import json,os,sys
from pathlib import Path
Path(os.environ['TEST_RECORD']).write_text(json.dumps({'args':sys.argv[1:],
'entry':os.environ.get('WIGITS_EXTENSIONS_ENTRY'),
'prior':os.environ.get('WIGITS_EXTENSIONS_PRIOR_RESULTS'),
'dependency':os.environ.get('SBATCH_DEPENDENCY')}))
''')
        (workflow/'run_optional_extensions.sh').write_text('python3 "$TEST_RECORDER" "$@"\n')
        record = base/'record.json'
        env = dict(os.environ, BASH_ENV=str(modules), WIGITS_RECOVERY_WORKFLOW=str(workflow),
                   TEST_RECORDER=str(recorder), TEST_RECORD=str(record), SBATCH_DEPENDENCY='afterany:456')
        result = subprocess.run(['bash',str(ROOT/'workflow/run_optional_recovery.sbatch'),str(run),entry,str(base/'predictions')],
                                env=env,text=True,capture_output=True)
        self.assertEqual(native.read_text(), json.dumps({'status':native_status}))
        self.assertEqual(ch.read_text(), '{"status":"COMPLETED"}\n')
        return result, record, run, config

    def test_head_selects_prediction_only_and_preserves_prior_ch(self):
        result, record, run, config = self.run_head('COMPLETED')
        self.assertEqual(result.returncode, 0, result.stderr)
        invocation = json.loads(record.read_text())
        self.assertEqual(invocation['entry'], 'PREDICT_EXISTING_NEO')
        self.assertEqual(invocation['prior'], str(run/'extensions/results'))
        self.assertIsNone(invocation['dependency'])
        self.assertEqual(invocation['args'], [str(config),str(run/'extensions/recovery/native-intake-v1.3.0/neo-predict'),
                         '--case_id','SYNTHETIC','--existing_export',str(run/'extensions/results/SYNTHETIC/neo_export')])

    def test_head_cannot_start_when_native_did_not_complete(self):
        result, record, _, _ = self.run_head('FAILED')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Native completion is required', result.stderr)
        self.assertFalse(record.exists())

    def test_prepare_selects_export_only_without_existing_export(self):
        result, record, run, config = self.run_head('COMPLETED','neo-prepare')
        self.assertEqual(result.returncode,0,result.stderr)
        invocation=json.loads(record.read_text())
        self.assertEqual(invocation['entry'],'PREPARE_EXISTING_NEO')
        self.assertNotIn('--existing_export',invocation['args'])
        self.assertEqual(invocation['args'],[str(config),
            str(run/'extensions/recovery/native-intake-v1.3.0/neo-prepare'),'--case_id','SYNTHETIC'])

    def test_review_entry_passes_existing_predictions_without_prediction(self):
        result, record, run, config = self.run_head('COMPLETED', 'neo-review')
        self.assertEqual(result.returncode, 0, result.stderr)
        invocation = json.loads(record.read_text())
        self.assertEqual(invocation['entry'], 'REVIEW_EXISTING_NEO')
        self.assertIn('--existing_predictions', invocation['args'])
        self.assertEqual(invocation['args'][1],str(run/'extensions/recovery/native-intake-v1.3.0/neo-review'))


if __name__ == '__main__':
    unittest.main()
