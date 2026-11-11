import sys
from pathlib import Path
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflow'))
from probe_sage_context import stage_redux_sidecars


class ProbeStagingTests(unittest.TestCase):
    def test_separate_bam_dirs_share_native_sidecar_directory(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);config={}
            for role,sample in [('normal','N'),('tumor','T')]:
                parent=root/sample;parent.mkdir()
                config[role+'_id']=sample;config[role+'_bam']=str(parent/(sample+'.redux.bam'))
                for suffix in ('bqr.tsv','jitter_params.tsv','ms_table.tsv.gz'):
                    (parent/f'{sample}.redux.{suffix}').write_bytes(b'native measurement')
            out=root/'probe';out.mkdir();target,sources=stage_redux_sidecars(config,out)
            self.assertEqual(len(sources),6)
            for name,source in sources.items():
                self.assertTrue((target/name).is_symlink())
                self.assertEqual((target/name).resolve(),Path(source))
                self.assertEqual(Path(source).read_bytes(),b'native measurement')

    def test_missing_sidecar_not_replaced_by_skipping_bqr(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(FileNotFoundError):
                stage_redux_sidecars(dict(normal_id='N',normal_bam=d+'/N.bam'),Path(d))


if __name__=='__main__':unittest.main()
