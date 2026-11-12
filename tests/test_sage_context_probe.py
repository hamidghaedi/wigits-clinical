import sys
from pathlib import Path
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflow'))
from probe_sage_context import stage_redux_sidecars, evidence_command


class ProbeStagingTests(unittest.TestCase):
    def test_evidence_retains_native_measurement_and_separates_tinc(self):
        cmd=['apptainer','exec','image','sage','-Xmx24g','-run_tinc',
             '-write_fit_variants','-jitter_bqr_dir','/both/',
             '-specific_var_file','requested.csv','-output_vcf','evidence.vcf.gz']
        out=evidence_command(cmd,Path('/isolated/sage.jar'))
        self.assertIn('-Ddeepsomatic.evidence_only=true',out)
        self.assertNotIn('-run_tinc',out)
        self.assertIn('-jitter_bqr_dir',out)
        self.assertIn('requested.csv',out)
        self.assertEqual(cmd[3],'sage')
        self.assertFalse(any('disable' in x for x in out))

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
