"""The existing owner launcher selects reusable post-native entry points."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class OptionalEntryLauncherTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.run = self.root/'original'
        (self.run/'package').mkdir(parents=True)
        (self.run/'state').mkdir()
        (self.run/'state/head_job.txt').write_text('777\n')
        (self.run/'package/site.env').write_text('SLURM_ACCOUNT_CPU=synthetic\nexport SBATCH_DEPENDENCY=afterany:999\n')
        self.workflow = self.root/'new/workflow'
        self.workflow.mkdir(parents=True)
        (self.workflow/'run_optional_recovery.sbatch').touch()
        commands = self.root/'bin'
        commands.mkdir()
        capture = self.root/'capture.json'
        (commands/'sbatch').write_text('#!/usr/bin/env python3\nimport json,os,sys\n'
            'json.dump({"args":sys.argv[1:],"dependency":os.environ.get("SBATCH_DEPENDENCY")},open(os.environ["CAPTURE"],"w"))\n'
            'print("12345")\n')
        (commands/'squeue').write_text('#!/bin/sh\nexit 0\n')
        for p in commands.iterdir():p.chmod(0o755)
        self.env = {k:v for k,v in os.environ.items() if k!='BASH_ENV'}
        self.env.update(PATH=str(commands)+os.pathsep+self.env['PATH'],CAPTURE=str(capture))
        self.capture = capture

    def launch(self, entry='neo-predict', extra=()):
        return subprocess.run(['bash',str(ROOT/'workflow/run_native_fir.sh'),'--entry',entry,
            '--run-dir',str(self.run),'--workflow-dir',str(self.workflow),*extra],
            env=self.env,text=True,capture_output=True)

    def test_prediction_uses_pinned_recovery_source_and_preserves_native_receipt(self):
        result = self.launch()
        self.assertEqual(result.returncode,0,result.stderr)
        capture=json.loads(self.capture.read_text())
        args=capture['args']
        self.assertIn(str(self.workflow/'run_optional_recovery.sbatch'),args)
        self.assertNotIn('run_full_native.sbatch',' '.join(args))
        self.assertIn('--time=1-00:00:00',args)
        self.assertIsNone(capture['dependency'])
        self.assertEqual((self.run/'state/head_job.txt').read_text(),'777\n')
        recovery=self.run/'extensions/recovery/native-intake-v1.2.0/neo-predict/state'
        self.assertEqual((recovery/'head_job.txt').read_text(),'12345\n')
        self.assertEqual((recovery/'entry.txt').read_text(),'neo-predict\n')

    def test_optional_entry_retains_dependency_and_walltime_overrides(self):
        result=self.launch(extra=['--dependency','afterok:42','--controller-time','08:00:00'])
        self.assertEqual(result.returncode,0,result.stderr)
        args=json.loads(self.capture.read_text())['args']
        self.assertIn('--dependency=afterok:42',args)
        self.assertIn('--time=08:00:00',args)

    def test_review_accepts_explicit_completed_predictions(self):
        predictions=self.root/'completed_predictions';predictions.mkdir()
        result=self.launch('neo-review',['--existing-predictions',str(predictions)])
        self.assertEqual(result.returncode,0,result.stderr)
        args=json.loads(self.capture.read_text())['args']
        self.assertEqual(args[-3:],[str(self.run),'neo-review',str(predictions)])

    def test_review_without_predictions_never_submits(self):
        result=self.launch('neo-review')
        self.assertNotEqual(result.returncode,0)
        self.assertFalse(self.capture.exists())

    def test_repeated_recovery_requires_resume(self):
        self.assertEqual(self.launch().returncode,0)
        self.assertNotEqual(self.launch().returncode,0)
        self.assertEqual(self.launch(extra=['--resume']).returncode,0)

    def test_unknown_entry_never_submits(self):
        self.assertNotEqual(self.launch('invented').returncode,0)
        self.assertFalse(self.capture.exists())


if __name__ == '__main__':
    unittest.main()
