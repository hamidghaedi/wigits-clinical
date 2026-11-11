import importlib.util
from pathlib import Path
import shlex
import unittest

spec=importlib.util.spec_from_file_location('recovery',Path(__file__).resolve().parents[1]/'workflow/prepare_run_recovery.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class RecoveryTests(unittest.TestCase):
    def test_preserves_original_command_except_explicit_resume_and_last_config(self):
        original='nextflow run nf-core/oncoanalyser -r 3.0.0 -c original.config --input /run/samples.csv -work-dir /run/work -resume -c resources.config'
        new=shlex.split(m.resume_command(original,'original-session',Path('/recovery/fix.config')))
        expected=shlex.split(original);expected.insert(expected.index('-resume')+1,'original-session')
        self.assertEqual(new,expected+['-c','/recovery/fix.config'])

    def test_replaces_existing_session_without_new_pipeline(self):
        cmd=m.resume_command('nextflow run nf-core/oncoanalyser -resume old -r 3.0.0','saved',Path('/fix'))
        self.assertIn('-resume saved',cmd)
        self.assertNotIn(' old ',cmd)

    def test_rejects_unrelated_launch(self):
        with self.assertRaises(ValueError):m.resume_command('nextflow run other -resume','id',Path('/fix'))
