import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

SCRIPT = Path(__file__).resolve().parents[1]/'workflow/report_resources.py'
spec=importlib.util.spec_from_file_location('report_resources',SCRIPT)
report=importlib.util.module_from_spec(spec);spec.loader.exec_module(report)


def parent(job='101',state='COMPLETED',cpu='00:32:00',elapsed='60'):
    return dict(Cluster='synthetic',JobIDRaw=job,JobName='example',State=state,
                ExitCode='0:0',AllocCPUS='64',AllocTRES='cpu=64,mem=24G,node=1',
                ElapsedRaw=elapsed,TotalCPU=cpu,MaxRSS='',Eligible='2026-01-01T00:00:00',
                Start='2026-01-01T00:00:20')


class ResourceReportTests(unittest.TestCase):
    def test_completion_hook_preserves_native_exit_even_if_reporting_fails(self):
        launcher=SCRIPT.with_name('run_oncoanalyser_frontenac.sh').read_text()
        tail='status=0\n'+launcher.split('status=0\n',1)[1]
        for native_exit in [0,7]:
            for report_exit in [0,1]:
                with self.subTest(native_exit=native_exit,report_exit=report_exit), tempfile.TemporaryDirectory() as tmp:
                    root=Path(tmp);(root/'logs').mkdir()
                    marker=root/'called.json'
                    (root/'report_resources.py').write_text('import json,pathlib,sys\n'
                        f'pathlib.Path({str(marker)!r}).write_text(json.dumps(sys.argv[1:]))\n'
                        f'raise SystemExit({report_exit})\n')
                    prefix=f'set -e\nHERE={str(root)!r}\nRUN_DIR={str(root)!r}\nSLURM_JOB_ID=101\nCMD=(bash -c "exit {native_exit}")\n'
                    result=subprocess.run(['bash','-c',prefix+tail],cwd=root,capture_output=True,text=True)
                    self.assertEqual(result.returncode,native_exit,result.stderr)
                    self.assertEqual(json.loads(marker.read_text()),['--run-dir',str(root),'--job-id','101'])

    def test_slurm_times_and_memory(self):
        self.assertEqual(report.seconds('1-12:46:17'),132377)
        self.assertEqual(report.seconds('00:00.250'),.25)
        self.assertEqual(report.gib('1048576K'),1)
        self.assertIsNone(report.gib(''))

    def test_parent_not_added_to_steps_and_memory_uses_batch(self):
        p=parent()
        batch={**p,'JobIDRaw':'101.batch','MaxRSS':'12G'}
        rows=report.usage_rows({'101':{}},[p,batch],[])
        self.assertEqual(len(rows),1)
        self.assertEqual(rows[0]['cpu_seconds'],1920)
        self.assertEqual(rows[0]['cpu_use_percent'],50)
        self.assertEqual(rows[0]['average_busy_cpus'],32)
        self.assertEqual(rows[0]['rss_fraction_percent'],50)
        self.assertEqual(rows[0]['queue_seconds'],20)

    def test_live_zero_cpu_is_not_zero_efficiency(self):
        p=parent(state='RUNNING',cpu='00:00:00')
        row=report.usage_rows({'101':{}},[p],[])[0]
        self.assertEqual(row['cpu_use_percent'],'')
        live=[dict(JobID='101.batch',NTasks='1',AveCPU='00:16:00',MaxRSS='6G')]
        row=report.usage_rows({'101':{}},[p],live)[0]
        self.assertEqual(row['cpu_use_percent'],25)
        self.assertEqual(row['cpu_measurement'],'live_batch_estimate')

    def test_unsuccessful_costs_and_gpu_counts(self):
        p=parent(state='CANCELLED by 1',cpu='00:00:00')
        p['AllocTRES']='cpu=64,mem=24G,node=1,gres/gpu=2,gres/gpu:a100=2'
        row=report.usage_rows({'101':{}},[p],[])[0]
        self.assertEqual(row['cpu_use_percent'],'')
        self.assertAlmostEqual(row['allocated_core_hours'],64/60,places=4)
        self.assertEqual(row['allocated_gpus'],2)
        self.assertAlmostEqual(row['gpu_card_hours'],2/60,places=4)

    def test_repeated_submission_logs_and_trace_count_each_job_once(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            (root/'logs').mkdir();(root/'results/pipeline_info').mkdir(parents=True)
            line='[SLURM] submitted process NFCORE_TEST:ALIGN (example) > jobId: 101; workDir: /none\n'
            (root/'.nextflow.log').write_text(line+line)
            (root/'logs/nextflow_head_old.log').write_text(line)
            (root/'results/pipeline_info/execution_trace.txt').write_text('native_id\tname\t%cpu\tpeak_rss\n101\tNFCORE_TEST:ALIGN (example)\t3200%\t4 GB\n')
            jobs=report.discover(root)
            self.assertEqual(list(jobs),['101'])
            self.assertEqual(jobs['101']['nextflow_cpu_percent'],'3200%')

    def test_html_escapes_labels_and_uses_weighted_cpu_use(self):
        rows=report.usage_rows({'101':{'task':'<script>'},'102':{}},
                               [parent(),parent(job='102',cpu='00:00:00',elapsed='120')],[])
        rendered=report.render(Path('/synthetic'),rows,[],'2026-01-01')
        self.assertIn('&lt;script&gt;',rendered)
        self.assertNotIn('<script>',rendered)
        self.assertIn('16.67',rendered)  # 1920/(64*(60+120))*100, not mean percentages.

    def test_offline_cli_writes_reports_without_scheduler(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            raw=root/'input.psv'
            p=parent()
            raw.write_text('|'.join(p)+'\n'+'|'.join(p.values())+'\n')
            result=subprocess.run(['python3',str(SCRIPT),'--run-dir',str(root),'--job-id','101',
                                   '--sacct-file',str(raw)],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            out=root/'results/pipeline_info/resource_usage'
            self.assertTrue((out/'resource_usage.csv').exists())
            self.assertTrue((out/'resource_usage.html').exists())
            self.assertIn('50.0',(out/'resource_usage.csv').read_text())


if __name__=='__main__':unittest.main()
