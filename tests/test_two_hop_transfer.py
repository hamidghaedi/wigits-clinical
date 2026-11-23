import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch


class TwoHopTests(unittest.TestCase):
    def test_resumes_existing_task_and_completes_hops_in_order(self):
        script = Path(__file__).resolve().parents[1] / 'workflow/transfer_fastqs_two_hop.py'
        spec = importlib.util.spec_from_file_location('transfer', script)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            dest = root/'reads.fastq.gz'
            dest.write_bytes(b'test')
            manifest = root/'manifest.json'
            manifest.write_text(json.dumps(dict(case='synthetic', source_endpoint='nas',
                relay_endpoint='relay', destination_endpoint='local', objects=[
                    dict(source='/source/reads.fastq.gz', relay='/relay/reads.fastq.gz',
                         destination=str(dest), expected_bytes=4)])))
            state = root/'state'
            state.mkdir()
            (state/'state.json').write_text(json.dumps({'files': {
                '/source/reads.fastq.gz': {'nas_to_fir_task': 'existing'}}}))
            calls = []

            def fake(*args):
                calls.append(args)
                if args[0] == 'transfer':
                    self.assertIn('--verify-checksum', args)
                    self.assertTrue(args[1].startswith('relay:'))
                    return {'task_id': 'second'}
                src, dst = ('nas', 'relay') if args[2] == 'existing' else ('relay', 'local')
                return dict(status='SUCCEEDED', source_endpoint_id=src,
                            destination_endpoint_id=dst, files_transferred=1, files_skipped=0)

            argv = [str(script), '--manifest', str(manifest), '--state-dir', str(state)]
            with patch.object(sys, 'argv', argv), patch.object(module, 'globus', fake):
                module.main()
            self.assertEqual([x[0] for x in calls], ['task', 'transfer', 'task'])
            result = json.loads((state/'state.json').read_text())
            self.assertEqual(result['status'], 'ALL_FASTQS_VERIFIED_LOCAL')
            self.assertTrue(dest.exists())


if __name__ == '__main__':
    unittest.main()
