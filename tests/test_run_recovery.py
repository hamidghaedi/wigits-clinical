import importlib.util
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest

spec=importlib.util.spec_from_file_location('recovery',Path(__file__).resolve().parents[1]/'workflow/prepare_run_recovery.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class RecoveryTests(unittest.TestCase):
    def test_rna_resource_recovery_preserves_launch_and_only_overlays_rna(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); run=root/'run'; run.mkdir(); out=root/'recovery'
            original='nextflow run /synthetic/oncoanalyser -c original.config --input sample.csv -work-dir work -resume'
            (run/'launch.sh').write_text(original+'\n')
            (run/'samplesheet.csv').write_text('sample\nsynthetic\n')
            subprocess.run([sys.executable,str(Path(m.__file__)), '--run',str(run),
                            '--out',str(out),'--route','rna','--rna-resources','--session','saved',
                            '--native-source','/synthetic/oncoanalyser'],check=True)
            config=(out/'recovery.config').read_text()
            self.assertIn('--READ_NAME_REGEX null',config)
            self.assertIn('RNA input-v5',config)
            self.assertNotIn('BWAMEM2',config)
            self.assertNotIn('REDUX',config)
            args=shlex.split((out/'resume_command.sh').read_text().split('exec ',1)[1])
            expected=shlex.split(original);expected.insert(expected.index('-resume')+1,'saved')
            self.assertEqual(args,expected+['-c',str(out/'recovery.config')])

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

    def test_local_source_requires_explicit_original_source_binding(self):
        original='nextflow run /synthetic/oncoanalyser/3.0.0 -c original -resume'
        with self.assertRaises(ValueError):
            m.resume_command(original,'saved',Path('/fix'))
        args=shlex.split(m.resume_command(original,'saved',Path('/fix'),
                                         Path('/synthetic/oncoanalyser/3.0.0')))
        self.assertEqual(args[:3],shlex.split(original)[:3])
        self.assertEqual(args[-2:],['-c','/fix'])
        self.assertEqual(args[args.index('-resume')+1],'saved')

    def test_explicit_local_source_does_not_accept_a_different_checkout(self):
        with self.assertRaises(ValueError):
            m.resume_command('nextflow run /other/main.nf -resume','id',Path('/fix'),
                             Path('/synthetic/oncoanalyser/3.0.0'))

    def test_esvee_route_changes_only_esvee_runtime_and_keeps_recorded_launch(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);run=root/'run';run.mkdir();out=root/'recovery'
            original='nextflow run /synthetic/oncoanalyser -c resources.config --input sample.csv -work-dir work -resume'
            (run/'launch.sh').write_text(original+'\n')
            (run/'samplesheet.csv').write_text('patient,sample\nSYNTHETIC,T\n')
            subprocess.run([sys.executable,str(Path(m.__file__)), '--run',str(run),
                            '--out',str(out),'--route','esvee','--session','saved',
                            '--native-source','/synthetic/oncoanalyser'],check=True)
            config=(out/'recovery.config').read_text()
            self.assertIn("withName: 'ESVEE'",config)
            self.assertIn(str(out/'esvee_jvm_recovery.env'),config)
            for other in ['CIDER','REDUX','BWAMEM2_ALIGN','SAGE','reference']:
                self.assertNotIn(other,config)
            expected=shlex.split(original);expected.insert(expected.index('-resume')+1,'saved')
            actual=shlex.split((out/'resume_command.sh').read_text().split('exec ',1)[1])
            self.assertEqual(actual,expected+['-c',str(out/'recovery.config')])
            self.assertEqual((out/'esvee_jvm_recovery.env').read_text(),
                             (Path(m.__file__).parent/'conf/esvee_jvm_recovery.env').read_text())

    def test_kraken_reuse_keeps_jar_setting_and_requires_all_completed_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            with self.assertRaises(FileNotFoundError):m.kraken_config(root,'SYNTHETIC')
            for suffix in ['readnames.txt','kraken2.report.all.txt','viralreferencescontigs.txt']:
                (root/f'SYNTHETIC.virusbreakend.vcf.{suffix}').write_text('synthetic\n')
            config=m.kraken_config(root,'SYNTHETIC')
            self.assertEqual(config.count('ln -sfn'),3)
            self.assertNotIn('ext.args',config)
            for p in root.iterdir():self.assertEqual(p.read_text(),'synthetic\n')
            self.assertIn("--jar /synthetic/runtime.jar",m.kraken_config(root,'SYNTHETIC',Path('/synthetic/runtime.jar')))
