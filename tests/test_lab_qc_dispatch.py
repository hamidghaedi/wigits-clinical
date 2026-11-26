import csv
import importlib.util
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

SPEC=importlib.util.spec_from_file_location('lab_dispatch',Path(__file__).parents[1]/'workflow/dispatch_laboratory_qc.py')
tool=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(tool)


class LabDispatchTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.root=Path(self.tmp.name)
        self.ref=self.root/'ref.fa';self.ref.write_text('synthetic-reference')
        for suffix in ['.fai','.dict']:Path(str(self.ref)+suffix).write_text('synthetic-metadata')
        self.queue=self.root/'queue.json';self.tracker=self.root/'tracker.json';self.lab=self.root/'laboratory.csv'
        self.config=dict(state_directory=str(self.root/'qc_state'),run_root=str(self.root/'runs'),
            native_tracker_config=str(self.tracker),queue_roots=[],laboratories=[str(self.lab)],baseline_job_ids=[],
            workflow_version='lab-qc/v0.2.0',workflow_commit='a'*40,workflow_directory='/synthetic/workflow',
            account='synthetic_account',reference_profiles=[dict(reference=str(self.ref),fai=str(self.ref)+'.fai',dict=str(self.ref)+'.dict')],
            logo='/synthetic/logo.jpg',css='/synthetic/report.css')
        self.cases=[];self.labrows=[]
        for i in range(2):self.add_case(i)
        tool.save(self.queue,dict(cases=self.cases));tool.save(self.tracker,dict(queues=[dict(path=str(self.queue),state_directory=str(self.root/'state'))]))
        self.write_lab()

    def add_case(self,i):
        case='CASE'+str(i);sid='LIB'+str(i);run=self.root/case;pkg=run/'package';pkg.mkdir(parents=True)
        self.cases.append(dict(case_id=case,package=str(pkg)))
        self.labrows.append(dict(sample_submitter_id=sid,fastq_parent_dir='/nas/SEQ_RUN/fastq',analyte_type='DNA',lanes='1',phase='Phase 2'))
        tool.save(pkg/'input_plan.json',dict(case_id=case,run_directory=str(run),
            libraries={'tumor':dict(sample_submitter_id=sid,analyte_type='DNA',fastq_parent_dir='/nas/SEQ_RUN/fastq')},
            objects=[dict(sample_id=sid,run_id='SEQ_RUN')]))
        tool.save(pkg/'release.json',dict(native_source_version='synthetic-native'))
        tool.save(self.root/'state'/(case+'.json'),dict(analysis_job=str(100+i)))
        paths=[]
        for mate in [1,2]:
            p=run/f'{sid}_S1_L001_R{mate}_001.fastq.gz';p.write_bytes(b'synthetic FASTQ placeholder');paths.append(str(p))
        with (pkg/'samplesheet.csv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=['sample_id','sequence_type','filetype','filepath','info']);w.writeheader()
            w.writerow(dict(sample_id=sid,sequence_type='dna',filetype='fastq',filepath=';'.join(paths),info='lane:001'))
        traces=[]
        for process,head in [('REDUX','aa/111111'),('BAMTOOLS','bb/222222')]:
            work=run/'work'/head;out=work/'outputs';out.mkdir(parents=True)
            (work/'versions.yml').write_text(json.dumps({f'PIPELINE:{process}':{process.lower():'1.0'}}))
            if process=='REDUX':
                for ending in ['.redux.bam','.redux.bam.bai']:
                    name=sid+ending;(out/name).write_bytes(b'complete producer output')
                    dest=run/'results'/case/'alignments'/sid/name;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes((out/name).read_bytes())
            else:
                (work/'ref.fa').symlink_to(self.ref);(work/'.command.sh').write_text('bamtools -ref_genome ref.fa\n')
                name=sid+'.bam_metric.flag_counts.tsv';(out/name).write_text('native flags')
                dest=run/'results'/case/'bamtools'/sid/name;dest.parent.mkdir(parents=True);dest.write_text('native flags')
            traces.append(dict(name=f'PIPELINE:{process} ({case}_{sid})',hash=head,status='COMPLETED',exit='0'))
        # A running caller must not prevent early QC. There is no final combined versions file.
        traces.append(dict(name='PIPELINE:SAGE_SOMATIC (paired)',hash='cc/333333',status='RUNNING',exit='-'))
        trace=run/'results/pipeline_info/execution_trace_test.txt';trace.parent.mkdir(parents=True)
        with trace.open('w') as f:
            w=csv.DictWriter(f,fieldnames=['name','hash','status','exit'],delimiter='\t');w.writeheader();w.writerows(traces)

    def write_lab(self):
        with self.lab.open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(self.labrows[0]));w.writeheader();w.writerows(self.labrows)

    def test_complete_run_submits_once_while_variant_calling_runs(self):
        calls=[]
        def submit(cmd,**kwargs):calls.append(cmd);return SimpleNamespace(stdout='999;cluster\n')
        a=tool.tick(self.config,submit);b=tool.tick(self.config,submit)
        self.assertEqual(len(calls),1);self.assertEqual(a[0]['job_id'],'999');self.assertEqual(b[0]['status'],'submitted')
        state=tool.read(Path(self.config['state_directory'])/'SEQ_RUN.run.json')
        manifest=tool.read(Path(state['execution_root'])/'manifest.json')
        self.assertEqual(len(manifest['samples']),2)
        versions=tool.read(manifest['samples'][0]['source_versions']);self.assertIn('BAMTOOLS',versions)
        self.assertNotIn('SAGE',json.dumps(manifest))

    def test_missing_library_prevents_partial_run_report(self):
        tool.save(self.queue,dict(cases=self.cases[:1]))
        result=tool.tick(self.config,lambda *a,**k:self.fail('must not submit'))
        self.assertEqual(result[0]['status'],'waiting_for_bams');self.assertIn('LIB1',result[0]['reason'])

    def test_existing_coverage_summary_is_bound_without_scanning_bam(self):
        parent=self.root/'CASE0/results/CASE0/bamtools/LIB0'
        source=parent/'LIB0.bam_metric.summary.tsv'
        source.write_text('MeanCoverage\n141\n')
        tool.tick(self.config,lambda *a,**k:SimpleNamespace(stdout='999\n'))
        state=tool.read(Path(self.config['state_directory'])/'SEQ_RUN.run.json')
        manifest=tool.read(Path(state['execution_root'])/'manifest.json')
        rows={r['sample_id']:r for r in manifest['samples']}
        self.assertEqual(rows['LIB0']['bamtools_summary'],str(source))
        self.assertIsNone(rows['LIB1']['bamtools_summary'])

    def test_incomplete_published_bam_copy_waits(self):
        p=self.root/'CASE1/results/CASE1/alignments/LIB1/LIB1.redux.bam';p.write_bytes(b'partial')
        result=tool.tick(self.config,lambda *a,**k:self.fail('must not submit'))
        self.assertIn('copy incomplete',result[0]['reason'])

    def test_failed_bam_task_never_uses_an_old_completed_attempt(self):
        p=self.root/'CASE1/results/pipeline_info/execution_trace_test.txt'
        p.write_text(p.read_text().replace('COMPLETED','FAILED',1))
        result=tool.tick(self.config,lambda *a,**k:self.fail('must not submit'))
        self.assertIn('REDUX not complete',result[0]['reason'])

    def test_old_jobs_are_not_retroactively_launched(self):
        self.config['baseline_job_ids']=['100','101']
        self.assertEqual(tool.tick(self.config,lambda *a,**k:self.fail('must not submit')),[])

    def test_active_external_qc_job_is_not_changed_or_duplicated(self):
        tool.save(Path(self.config['state_directory'])/'SEQ_RUN.run.json',dict(run_id='SEQ_RUN',external_managed=True,job_id='777'))
        result=tool.tick(self.config,lambda *a,**k:self.fail('must not submit'))
        self.assertEqual(result[0],dict(run_id='SEQ_RUN',status='external_managed',job_id='777'))

    def test_unknown_submission_response_is_not_retried(self):
        calls=[]
        def submit(*a,**k):calls.append(a);return SimpleNamespace(stdout='ambiguous response')
        with self.assertRaises(ValueError):tool.tick(self.config,submit)
        result=tool.tick(self.config,submit)
        self.assertEqual(len(calls),1);self.assertEqual(result[0]['status'],'submission_unconfirmed')

    def test_pointer_enrolls_a_new_queue_without_tracker_changes(self):
        self.config['queue_roots']=[str(self.root/'queues')]
        tool.save(self.root/'queues/new/lab_qc_source.json',dict(path=str(self.queue),state_directory=str(self.root/'state'),laboratory=str(self.lab)))
        tool.save(self.tracker,dict(queues=[]))
        self.assertEqual(len(tool.candidates(self.config)),2)

    def test_missing_lane_is_not_relabelled_as_complete(self):
        self.labrows[0]['lanes']='1;2';self.write_lab()
        result=tool.tick(self.config,lambda *a,**k:self.fail('must not submit'))
        self.assertIn('not all laboratory lanes',result[0]['reason'])

    def test_delivery_happens_once_without_caller_or_head_success(self):
        tool.tick(self.config,lambda *a,**k:SimpleNamespace(stdout='999'))
        with patch.object(tool,'deliver',return_value=True) as deliver:
            result=tool.tick(self.config);tool.tick(self.config)
            self.assertEqual(deliver.call_count,1);self.assertEqual(result[0]['status'],'delivered')

    def test_failed_qc_head_is_visible_and_not_resubmitted(self):
        tool.tick(self.config,lambda *a,**k:SimpleNamespace(stdout='999'))
        result=tool.tick(self.config,lambda *a,**k:self.fail('must not resubmit'),
                         status_reader=lambda ids:{'999':['TIMEOUT','0:0']})
        self.assertEqual(result[0]['status'],'execution_failed')

    def test_delivery_archives_only_completed_qc_and_uses_existing_stager(self):
        tool.tick(self.config,lambda *a,**k:SimpleNamespace(stdout='999'))
        state=tool.read(Path(self.config['state_directory'])/'SEQ_RUN.run.json')
        self.config.update(archive_root=str(self.root/'archive'),python='/synthetic/python',stager='/synthetic/stager.py')
        result=Path(state['execution_root'])/'results/SEQ_RUN';result.mkdir(parents=True)
        for name in ['qc','provenance']:(result/name).mkdir()
        report=result/'SEQ_RUN.multiqc_report.html';report.write_text('synthetic report')
        tool.save(result/'_SUCCESS',dict(run_id='SEQ_RUN',report_sha256=tool.hashlib.sha256(report.read_bytes()).hexdigest()))
        (result/'SEQ_RUN.multiqc_report_data.zip').write_bytes(b'synthetic archive')
        with patch.object(tool.subprocess,'run',return_value=SimpleNamespace(stdout='{"reports":[]}')) as run:
            self.assertTrue(tool.deliver(self.config,state))
            self.assertIn('--manifest',run.call_args.args[0])
        self.assertTrue(list((self.root/'archive/Phase_2/qc_reports/SEQ_RUN/evidence').glob('*/qc_sources.tar.gz')))


if __name__=='__main__':unittest.main()


class BamEventTests(unittest.TestCase):
    setUp=LabDispatchTests.setUp
    add_case=LabDispatchTests.add_case
    write_lab=LabDispatchTests.write_lab
    # Reuse the synthetic native output tree; no BAM content is inspected.
    def event(self,i):
        run=self.root/('CASE'+str(i));sid='LIB'+str(i)
        work=run/'work/bb/222222';metrics=work/'outputs'
        (work/'.exitcode').write_text('0')
        (work/(sid+'.redux.bam')).symlink_to(run/'work/aa/111111/outputs'/(sid+'.redux.bam'))
        (work/'.command.sh').write_text('bamtools -ref_genome ref.fa -bam_file '+sid+'.redux.bam\n')
        return metrics

    def test_last_bam_event_triggers_without_published_outputs_or_traces(self):
        import shutil
        e0=self.event(0);e1=self.event(1)
        for case in ['CASE0','CASE1']:shutil.rmtree(self.root/case/'results')
        calls=[]
        def submit(*a,**k):calls.append(a);return SimpleNamespace(stdout='999')
        first=tool.bam_ready(self.config,e0,submit)
        self.assertEqual(first[0]['status'],'waiting_for_bams');self.assertEqual(calls,[])
        second=tool.bam_ready(self.config,e1,submit)
        self.assertEqual(second[0]['status'],'submitted');self.assertEqual(len(calls),1)
        tool.bam_ready(self.config,e1,submit);self.assertEqual(len(calls),1)

    def test_failed_bamtools_never_emits_readiness(self):
        e=self.event(0);(e.parent/'.exitcode').write_text('1')
        with self.assertRaisesRegex(ValueError,'successful completed'):tool.bam_ready(self.config,e)
        self.assertFalse((Path(self.config['state_directory'])/'bam_events').exists())

    def test_report_event_delivers_its_own_run_before_other_run_finishes(self):
        tool.bam_ready(self.config,self.event(0),lambda *a,**k:SimpleNamespace(stdout='999'))
        tool.bam_ready(self.config,self.event(1),lambda *a,**k:SimpleNamespace(stdout='999'))
        state=tool.read(Path(self.config['state_directory'])/'SEQ_RUN.run.json')
        task=Path(state['execution_root'])/'work/report';task.mkdir(parents=True)
        with patch.object(tool,'deliver',return_value=True) as delivery:
            tool.finish(self.config,'SEQ_RUN',task)
            self.assertEqual(delivery.call_args.args[2],task)

    def test_successful_qc_completion_delivers_once_without_polling(self):
        tool.bam_ready(self.config,self.event(0),lambda *a,**k:SimpleNamespace(stdout='999'))
        tool.bam_ready(self.config,self.event(1),lambda *a,**k:SimpleNamespace(stdout='999'))
        with patch.object(tool,'deliver',return_value=True) as delivery:
            tool.finish(self.config,'SEQ_RUN');tool.finish(self.config,'SEQ_RUN')
            self.assertEqual(delivery.call_count,1)
