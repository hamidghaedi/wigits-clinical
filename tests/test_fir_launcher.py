"""Execute the actual head/launcher with synthetic FASTQs and stub executables."""
import csv
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class FirLauncherTests(unittest.TestCase):
    def test_fir_site_reaches_native_command_without_frontenac_bindings(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run = root / 'analysis/WiGiTS/vTEST/synthetic'
            package = run / 'package'
            package.mkdir(parents=True)
            (package / 'workflow').symlink_to(ROOT / 'workflow', target_is_directory=True)
            project = root / 'project'
            source = project / 'releases/oncoanalyser/3.0.0-qlmp.2.2'
            source.mkdir(parents=True)
            (source / 'main.nf').touch()
            ref = project / 'resources/oncoanalyser/3.0.0'
            ref.mkdir(parents=True)
            (ref / 'reference_data.config').touch()
            with (package / 'samplesheet.csv').open('w') as stream:
                writer = csv.writer(stream)
                writer.writerow('group_id,subject_id,sample_id,sample_type,sequence_type,filetype,info,filepath'.split(','))
                for role in ('normal', 'tumor'):
                    mates = [root / f'{role}_R{i}.fastq.gz' for i in (1, 2)]
                    for path in mates:
                        path.write_bytes(b'synthetic size proxy, not parsed')
                    writer.writerow(['synthetic', 'synthetic', role, role, 'dna', 'fastq',
                                     f'library_id:{role};lane:001;flowcell:SYNTHETIC', ';'.join(map(str, mates))])
            (package / 'site.env').write_text(
                f'WIGITS_PROJECT_ROOT="{project}"\nWIGITS_SCRATCH_ROOT="{root}/scratch"\n'
                f'WIGITS_RUN_ROOT="{run.parent}"\n'
                'source "$WIGITS_WORKFLOW_DIR/conf/fir_paths.env"\n')
            commands = root / 'bin'
            commands.mkdir()
            fake = commands / 'nextflow'
            fake.write_text('#!/usr/bin/env python3\nimport json,os,sys\n'
                            'if "-version" in sys.argv: print("nextflow version synthetic")\n'
                            'else:\n open(os.environ["CAPTURE"],"w").write(json.dumps(dict(argv=sys.argv[1:],'
                            'cache=os.environ["APPTAINER_CACHEDIR"],plugins=os.environ["NXF_PLUGINS_DIR"],'
                            'site=os.environ["WIGITS_EXTENSIONS_SITE_CONFIG"])))\n'
                            ' open(".nextflow.log","w").write("synthetic")\n')
            fake.chmod(0o755)
            apptainer = commands / 'apptainer'
            apptainer.write_text('#!/bin/sh\necho "apptainer synthetic"\n')
            apptainer.chmod(0o755)
            # Replace only the external module profile, retaining normal source semantics.
            wrapper = 'source() { if [[ "$1" != /cvmfs/soft.computecanada.ca/config/profile/bash.sh ]]; then builtin source "$@"; fi; }; module() { :; }; export -f source module; exec bash "$@"'
            env = {k: v for k, v in os.environ.items() if k != 'BASH_ENV' and not k.startswith(('WIGITS_', 'SLURM_', 'ONCO_', 'REF_CONFIG'))}
            env.update(PATH=str(commands)+os.pathsep+env['PATH'], CAPTURE=str(root/'capture.json'))
            reply = subprocess.run(['bash', '-c', wrapper, 'test', str(ROOT/'workflow/run_full_native.sbatch'),
                                    str(package), 'synthetic'], env=env, text=True, capture_output=True)
            self.assertEqual(reply.returncode, 0, reply.stdout+reply.stderr)
            data = json.loads((root/'capture.json').read_text())
            args = data['argv']
            self.assertIn(str(source), args)
            self.assertIn(str(package/'workflow/conf/fir.config'), args)
            self.assertNotIn(str(package/'workflow/conf/frontenac.config'), args)
            self.assertEqual(args[-2:], ['-c', str(package/'workflow/conf/fir_limits.config')])
            self.assertEqual(data['cache'], str(root/'scratch/.cache/apptainer'))
            self.assertEqual(data['plugins'], str(project/'resources/nextflow/plugins'))
            self.assertEqual(data['site'], str(package/'workflow/conf/fir.config'))
            resources = json.loads((run/'state/input_resources.json').read_text())
            self.assertEqual(len(resources['alignment_lanes']), 2)
            self.assertTrue(all(lane['cpus']==64 for lane in resources['alignment_lanes']))

    def test_owner_launcher_submits_once_with_ch_and_reuses_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);run=root/'runs/SYNTHETIC';workflow=root/'release/workflow'
            workflow.mkdir(parents=True);(workflow/'run_full_native.sbatch').touch()
            manifest=root/'samplesheet.csv';manifest.write_text('synthetic supplied manifest\n')
            ch=root/'blood_ch.json';ch.write_text(json.dumps(dict(cases=[dict(options=dict(blood_ch=True),laboratory=dict(normal_specimen_type='Blood'))])))
            config=root/'site.env';config.write_text(f'RUN_ROOT="{run.parent}"\nWIGITS_RUN_ID=SYNTHETIC\nWIGITS_WORKFLOW_DIR="{workflow}"\nSLURM_ACCOUNT_CPU=synthetic-account\nWIGITS_BLOOD_CH_CONFIG="{ch}"\n')
            commands=root/'bin';commands.mkdir()
            sbatch=commands/'sbatch';sbatch.write_text('#!/bin/sh\nprintf "%s\\n" "$@" > "$CAPTURE"\necho 12345\n');sbatch.chmod(0o755)
            env={k:v for k,v in os.environ.items() if k!='BASH_ENV'}
            env.update(PATH=str(commands)+os.pathsep+env['PATH'],CAPTURE=str(root/'submit.txt'))
            command=['bash',str(ROOT/'workflow/run_native_fir.sh'),'--manifest',str(manifest),'--config',str(config),'--blood-ch']
            first=subprocess.run(command,env=env,text=True,capture_output=True)
            self.assertEqual(first.returncode,0,first.stderr)
            self.assertEqual((run/'state/head_job.txt').read_text(),'12345\n')
            self.assertEqual((run/'package/samplesheet.csv').read_bytes(),manifest.read_bytes())
            self.assertEqual((run/'package/optional_extensions.json').read_bytes(),ch.read_bytes())
            self.assertIn('--account=synthetic-account',(root/'submit.txt').read_text())
            again=subprocess.run(command,env=env,text=True,capture_output=True)
            self.assertNotEqual(again.returncode,0)
            self.assertIn('Already submitted',again.stderr)
            manifest.write_text('changed manifest\n')
            changed=subprocess.run(command,env=env,text=True,capture_output=True)
            self.assertNotEqual(changed.returncode,0)
            self.assertIn('Existing package differs',changed.stderr)

    def test_controller_walltime_default_site_and_cli_override_reach_sbatch(self):
        for site_time,extra,expected in [('',[], '4-00:00:00'),
                ('3-00:00:00',[], '3-00:00:00'),
                ('3-00:00:00',['--controller-time','2-12:00:00'],'2-12:00:00')]:
            with self.subTest(site_time=site_time,extra=extra), tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);workflow=root/'release/workflow';workflow.mkdir(parents=True)
                (workflow/'run_full_native.sbatch').touch()
                manifest=root/'samplesheet.csv';manifest.write_text('prepared synthetic input\n')
                site=root/'site.env';site.write_text(f'RUN_ROOT="{root}/analysis"\nWIGITS_RUN_ID=SYNTHETIC\n'
                    f'WIGITS_WORKFLOW_DIR="{workflow}"\nSLURM_ACCOUNT_CPU=synthetic\n'
                    +(f'WIGITS_CONTROLLER_TIME={site_time}\n' if site_time else ''))
                commands=root/'bin';commands.mkdir()
                stub=commands/'sbatch';stub.write_text('#!/bin/sh\nprintf "%s\\n" "$@" > "$CAPTURE"\necho 12345\n');stub.chmod(0o755)
                env={k:v for k,v in os.environ.items() if k!='BASH_ENV' and not k.startswith('WIGITS_')}
                env.update(PATH=str(commands)+os.pathsep+env['PATH'],CAPTURE=str(root/'sbatch.args'))
                result=subprocess.run(['bash',str(ROOT/'workflow/run_native_fir.sh'),
                    '--manifest',str(manifest),'--config',str(site),*extra],env=env,text=True,capture_output=True)
                self.assertEqual(result.returncode,0,result.stderr)
                submitted=(root/'sbatch.args').read_text().splitlines()
                self.assertIn('--time='+expected,submitted)
                self.assertEqual(sum(arg.startswith('--time=') for arg in submitted),1)
                self.assertEqual((root/'analysis/SYNTHETIC/state/head_walltime.txt').read_text(),expected+'\n')
                self.assertEqual((root/'analysis/SYNTHETIC/package/site.env').read_bytes(),site.read_bytes())
                self.assertEqual((root/'analysis/SYNTHETIC/state/head_job.txt').read_text(),'12345\n')

    def test_missing_explicit_site_never_falls_back_silently(self):
        source = (ROOT/'workflow/run_full_native.sbatch').read_text()
        self.assertIn('if [[ -z "${WIGITS_SITE_ENV:-}" ]]', source)
        self.assertIn('source "$WIGITS_SITE_ENV"', source)
        self.assertNotIn('if [[ ! -f "$WIGITS_SITE_ENV" ]]', source)


if __name__ == '__main__':
    unittest.main()
