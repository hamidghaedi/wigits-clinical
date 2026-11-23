import importlib.util
from pathlib import Path
import tempfile
import gzip
import unittest

spec=importlib.util.spec_from_file_location('distribution',Path(__file__).resolve().parents[1]/'workflow/benchmark_deepsomatic.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class DistributionTests(unittest.TestCase):
    def test_partition_complete_and_disjoint(self):
        a=m.shard_ids(0,2,64);b=m.shard_ids(1,2,64)
        self.assertEqual(len(a),32);self.assertFalse(set(a)&set(b))
        self.assertEqual(sorted(a+b),list(range(64)))

    def test_scientific_options_preserved(self):
        cmd='time seq 0 63 | parallel -q --halt 2 --line-buffer /opt/deepvariant/bin/make_examples_somatic --mode calling --ref r.fa --reads_tumor t.bam --reads_normal n.bam --examples /old/intermediate/make_examples_somatic.tfrecord@64.gz --checkpoint /model --call_small_model_examples --task {}'
        argv=m.native_command(cmd,Path('/old'),Path('/new'))
        self.assertIn('--call_small_model_examples',argv)
        self.assertEqual(argv[argv.index('--checkpoint')+1],'/model')
        self.assertIn('/new/intermediate/make_examples_somatic.tfrecord@64.gz',argv)
        self.assertNotIn('/old', ' '.join(argv))
        self.assertEqual(argv[argv.index('--task')+1],'{}')
        with self.assertRaises(ValueError):m.native_command(cmd.replace('--mode calling','--regions chr17'),Path('/old'),Path('/new'))
        with self.assertRaises(ValueError):m.native_command(cmd.replace('seq 0 63','seq 0 31'),Path('/old'),Path('/new'))

    def test_digest_ignores_gzip_metadata(self):
        with tempfile.TemporaryDirectory() as d:
            a=Path(d)/'a.gz';b=Path(d)/'b.gz'
            a.write_bytes(gzip.compress(b'payload',mtime=1));b.write_bytes(gzip.compress(b'payload',mtime=2))
            self.assertNotEqual(a.read_bytes(),b.read_bytes())
            self.assertEqual(m.digest(a),m.digest(b))


if __name__=='__main__':unittest.main()
