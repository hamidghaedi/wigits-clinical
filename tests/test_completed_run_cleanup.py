import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

SOURCE = Path(__file__).resolve().parents[1] / 'workflow/cleanup_completed_run.py'
spec = importlib.util.spec_from_file_location('completed_cleanup', SOURCE)
cleanup = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cleanup)


class CompletedCleanupTest(unittest.TestCase):
    def fixture(self, root):
        run = root / 'run'
        inputs = run / 'inputs/R/T'
        inputs.mkdir(parents=True)
        (run / 'work/task').mkdir(parents=True)
        fastq = inputs / 'T_R1.fastq.gz'; fastq.write_bytes(b'restored')
        archive = inputs / 'T.tar.gz'; archive.write_bytes(b'archive')
        receipt = root / 'restore.json'
        receipt.write_text(json.dumps(dict(status='complete', files=[dict(path=str(fastq), bytes=8)])))
        libraries = [dict(sample_id='T', archive=str(archive), receipt=str(receipt))]
        plan = root / 'input.json'; plan.write_text(json.dumps(dict(run_directory=str(run), libraries=libraries)))
        results = run / 'results/C'
        (results / 'orange').mkdir(parents=True)
        (results / 'orange/T.orange.pdf').write_text('report prose')
        (results / 'alignments/T').mkdir(parents=True)
        (results / 'alignments/T/T.bam').write_bytes(b'final alignment')
        (run / 'results/multiqc_report.html').write_text('QC')
        return run, plan, libraries

    def test_scheduler_must_report_success(self):
        for text in ['RUNNING|0:0\n', 'FAILED|1:0\n', '']:
            with patch.object(cleanup.subprocess, 'run', return_value=subprocess.CompletedProcess([], 0, text, '')):
                with self.assertRaises(ValueError): cleanup.completed('1')
        with patch.object(cleanup.subprocess, 'run', return_value=subprocess.CompletedProcess([], 0, 'COMPLETED|0:0\n', '')):
            cleanup.completed('1')

    def test_large_candidates_skip_symlink_and_final_results(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); run, plan, libraries = self.fixture(root)
            bam = run / 'work/task/lane.bam'
            with bam.open('wb') as handle: handle.truncate(1024**3)
            (run / 'work/task/foreign.bam').symlink_to(run / 'results/C/alignments/T/T.bam')
            rows = cleanup.candidates(run, dict(libraries=libraries))
            self.assertEqual({Path(r['path']).name for r in rows}, {'T.tar.gz', 'T_R1.fastq.gz', 'lane.bam'})
            self.assertTrue((run / 'results/C/alignments/T/T.bam').is_file())

    def test_final_alignment_required_for_each_library(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); run, plan, libraries = self.fixture(root)
            with self.assertRaises(ValueError):
                cleanup.preserve_results(run, root / 'release', libraries + [dict(sample_id='N')])

    def test_apply_is_repeatable_and_preserves_results_and_evidence(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); run, inputs, libraries = self.fixture(root)
            plan = root / 'cleanup.json'
            args = ['cleanup', '--run-dir', str(run), '--input-plan', str(inputs), '--head-job', '1',
                    '--release-dir', str(root / 'release'), '--plan', str(plan), '--apply']
            with patch('sys.argv', args), patch.object(cleanup, 'completed'):
                cleanup.main()
                self.assertFalse((run / 'inputs/R/T/T_R1.fastq.gz').exists())
                self.assertTrue((run / 'results/C/alignments/T/T.bam').exists())
                self.assertEqual((root / 'release/results/C/orange/T.orange.pdf').read_text(), 'report prose')
                self.assertTrue(Path(libraries[0]['receipt']).exists())
                cleanup.main()
            result = json.loads(plan.with_suffix('.complete.json').read_text())
            self.assertEqual(result['removed_files'], 0)
            self.assertEqual(len(result['already_absent']), 2)

    def test_changed_file_stops_apply(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); run, inputs, libraries = self.fixture(root)
            plan = root / 'cleanup.json'
            args = ['cleanup', '--run-dir', str(run), '--input-plan', str(inputs), '--head-job', '1',
                    '--release-dir', str(root / 'release'), '--plan', str(plan)]
            with patch('sys.argv', args), patch.object(cleanup, 'completed'): cleanup.main()
            (run / 'inputs/R/T/T.tar.gz').write_bytes(b'changed size')
            with patch('sys.argv', args + ['--apply']), patch.object(cleanup, 'completed'):
                with self.assertRaises(ValueError): cleanup.main()
            self.assertTrue((run / 'results/C/alignments/T/T.bam').exists())


if __name__ == '__main__':
    unittest.main()
