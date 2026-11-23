import gzip
import importlib.util
from pathlib import Path
import shutil
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('existing_fastqs', Path(__file__).resolve().parents[1]/'workflow/check_existing_fastqs.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


@unittest.skipUnless(shutil.which('pigz'), 'pigz required')
class ExistingFastqTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.paths = [Path(self.tmp.name)/f'R{i}.fastq.gz' for i in (1, 2)]
        for i, p in enumerate(self.paths, 1):
            p.write_bytes(gzip.compress(f'@read/{i}\nAC\n+\nII\n'.encode()))

    def test_complete_pair(self):
        self.assertEqual(module.check_pair('synthetic', self.paths)['read_pairs'], 1)

    def test_mismatch(self):
        self.paths[1].write_bytes(gzip.compress(b'@other/2\nAC\n+\nII\n'))
        with self.assertRaisesRegex(ValueError, 'different mate names'):
            module.check_pair('synthetic', self.paths)

    def test_corruption(self):
        self.paths[1].write_bytes(self.paths[1].read_bytes()[:-4])
        with self.assertRaisesRegex(ValueError, 'gzip integrity'):
            module.check_pair('synthetic', self.paths)

    def test_unequal_counts(self):
        self.paths[1].write_bytes(gzip.compress(b''))
        with self.assertRaisesRegex(ValueError, 'unequal mate counts'):
            module.check_pair('synthetic', self.paths)

    def test_malformed_record(self):
        self.paths[1].write_bytes(gzip.compress(b'@read/2\nAC\n+\nI\n'))
        with self.assertRaisesRegex(ValueError, 'malformed'):
            module.check_pair('synthetic', self.paths)
