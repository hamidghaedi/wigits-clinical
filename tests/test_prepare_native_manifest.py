import csv
import gzip
import importlib.util
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('prepare_native',ROOT/'workflow/prepare_native_manifest.py')
tool=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(tool)
SPEC=importlib.util.spec_from_file_location('intake_lab_qc',ROOT/'workflow/dispatch_laboratory_qc.py')
qc=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(qc)


class PrepareNativeTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name);self.project=self.root/'project';self.scratch=self.root/'scratch'
        self.rows=[]
        for case in ('CASE_A','CASE_B'):
            for kind,run,analyte,descriptor,specimen in [('DNA','RUN_T','DNA','Primary','Tissue'),
                                                      ('gDNA','RUN_N','DNA','Normal','Blood')]:
                self.add_library(case,case+'-'+kind,run,analyte,descriptor,specimen)
        self.lab=self.root/'laboratory.csv';self.write_lab()
        self.config=dict(project_root=str(self.project),scratch_root=str(self.scratch),workflow_version='vTEST',
            workflow_directory=str(ROOT/'workflow'),workflow_commit='a'*40,native_source='/synthetic/native',
            native_source_version='synthetic',native_source_commit='b'*40,reference_root='/synthetic/reference',
            resource_version='synthetic',driver_panel_sha256='c'*64,
            extension_resources=dict(bind_roots=[str(self.root)]),
            lab_qc=dict(state_directory=str(self.project/'qc_state'),run_root=str(self.scratch/'qc'),
                baseline_job_ids=[],workflow_version='lab-qc/v0.4.1',workflow_commit='a'*40,
                workflow_directory=str(ROOT/'workflow'),account='synthetic_account',
                logo='/synthetic/logo',css='/synthetic/css',site_environment='/synthetic/shared-site.env'))
        self.site=self.root/'native.json';self.write_config()

    def write_config(self): self.site.write_text(json.dumps(self.config))

    def write_lab(self): tool.write_csv(self.lab,self.rows)

    def add_library(self,case,sample,run,analyte,descriptor,specimen):
        self.rows.append(dict(case_submitter_id=case,sample_submitter_id=sample,analyte_type=analyte,
            tumor_descriptor=descriptor,specimen_type=specimen,tissue_preservation_method='Frozen',
            lanes='1;2',fastq_parent_dir='/nas/'+run+'/fastq'))
        for lane in (1,2):
            self.pair(sample,run,lane)

    def pair(self,sample,run,lane,chunk=1):
        folder=self.scratch/'fastq'/run/sample;folder.mkdir(parents=True,exist_ok=True)
        for mate in (1,2):
            path=folder/f'{sample}_S1_L{lane:03d}_R{mate}_{chunk:03d}.fastq.gz'
            with gzip.open(path,'wt') as stream:
                stream.write(f'@TEST:1:FC_{run}:{lane}:1:1:1 {mate}:N:0:AAAA\nACGT\n+\nIIII\n')

    def args(self,case='CASE_A',extra=()):
        return tool.parser().parse_args(['--config',str(self.site),'--lab-manifest',str(self.lab),
            '--phase','Phase_1','--case-id',case,'--analysis-id',case+'-test',*extra])

    def prepared(self,case='CASE_A',extra=()):
        result=tool.prepare(self.args(case,extra));return Path(result['prepared'])

    def test_case_selection_preserves_complete_run_qc_membership(self):
        package=self.prepared(extra=['--blood-ch'])
        with (package/'samplesheet.csv').open() as stream: rows=list(csv.DictReader(stream))
        self.assertEqual(len(rows),4)
        self.assertEqual({r['sample_id'] for r in rows},{'CASE_A-DNA','CASE_A-gDNA'})
        config=json.loads((self.project/'file_manifests/lab_qc/config.json').read_text())
        for path in config['laboratories']:
            with Path(path).open() as stream: libraries=list(csv.DictReader(stream))
            self.assertEqual(len(libraries),2)
            self.assertTrue(any(r['sample_submitter_id'].startswith('CASE_B') for r in libraries))
        options=json.loads((package/'optional_extensions.json').read_text())['cases'][0]['options']
        self.assertEqual(options,dict(blood_ch=True,rna_neoantigens=False))
        self.assertFalse(list(self.project.rglob('*.fastq.gz')))
        self.assertFalse((self.scratch/'analysis').exists())

    def test_same_tool_prepares_another_case_with_arguments_only(self):
        a=self.prepared();b=self.prepared('CASE_B')
        for p,case in [(a,'CASE_A'),(b,'CASE_B')]:
            plan=json.loads((p/'input_plan.json').read_text())
            self.assertEqual(plan['case_id'],case)
            self.assertEqual(len(plan['objects']),8)
            self.assertEqual(len(list(p.glob('*.sh'))),0)
            self.assertEqual(len(list(p.glob('*.md'))),0)

    def test_missing_lane_or_mate_stops_preparation(self):
        folder=self.scratch/'fastq/RUN_T/CASE_A-DNA'
        next(folder.glob('*L002_R2*')).unlink()
        with self.assertRaisesRegex(ValueError,'missing mate'):self.prepared()
        next(folder.glob('*L002_R1*')).unlink()
        with self.assertRaisesRegex(ValueError,'expected lanes'):self.prepared()

    def test_all_chunks_are_included(self):
        self.pair('CASE_A-DNA','RUN_T',1,chunk=2)
        package=self.prepared()
        plan=json.loads((package/'input_plan.json').read_text())
        self.assertEqual(len(plan['objects']),10)
        self.assertEqual({o['chunk'] for o in plan['objects']},{'001','002'})

    def test_rna_and_neo_selection_are_explicit_and_preserved(self):
        self.add_library('CASE_A','CASE_A-RNA','RUN_R','RNA','Primary','Tissue');self.write_lab()
        package=self.prepared(extra=['--rna-neoantigens'])
        config=json.loads((package/'optional_extensions.json').read_text())['cases'][0]
        self.assertTrue(config['laboratory']['tumour_rna_available'])
        self.assertTrue(config['options']['rna_neoantigens'])
        self.assertEqual(config['rna_sample'],'CASE_A-RNA')

    def test_dna_only_does_not_infer_rna_or_neo(self):
        self.add_library('CASE_A','CASE_A-RNA','RUN_R','RNA','Primary','Tissue');self.write_lab()
        package=self.prepared(extra=['--dna-only'])
        config=json.loads((package/'optional_extensions.json').read_text())['cases'][0]
        self.assertFalse(config['laboratory']['tumour_rna_available'])
        self.assertEqual(config['rna_sample'],'')

    def test_ffpe_and_nonblood_are_not_silently_routed(self):
        self.rows[0]['tissue_preservation_method']='FFPE';self.write_lab()
        with self.assertRaisesRegex(ValueError,'fresh-frozen'):self.prepared()
        self.rows[0]['tissue_preservation_method']='Frozen';self.rows[1]['specimen_type']='Saliva';self.write_lab()
        with self.assertRaisesRegex(ValueError,'blood normal'):self.prepared(extra=['--blood-ch'])

    def test_ambiguous_pair_and_conflicting_laboratory_are_rejected(self):
        self.add_library('CASE_A','CASE_A-OTHER-DNA','RUN_T','DNA','Primary','Tissue');self.write_lab()
        with self.assertRaisesRegex(ValueError,'one tumor library'):self.prepared()
        package=self.prepared(extra=['--tumour-sample','CASE_A-DNA'])
        self.assertTrue(package.is_dir())

    def test_registered_full_run_cannot_be_shrunk_by_case_only_csv(self):
        self.prepared()
        self.rows=[r for r in self.rows if r['case_submitter_id']=='CASE_B'];self.write_lab()
        with self.assertRaisesRegex(ValueError,'shrink'):self.prepared('CASE_B')

    def test_runinfo_avoids_fastq_header_reads(self):
        infos=[]
        for run in ('RUN_T','RUN_N'):
            path=self.root/(run+'.xml');path.write_text(f'<RunInfo><Run Id="{run}"><Flowcell>EXACT_{run}</Flowcell></Run></RunInfo>')
            infos.extend(['--run-info',str(path)])
        with patch.object(tool,'header_flowcell',side_effect=AssertionError('must use supplied RunInfo')):
            package=self.prepared(extra=infos)
        self.assertIn('flowcell:EXACT_RUN_T',(package/'samplesheet.csv').read_text())

    def test_qc_waits_for_other_case_and_submits_once_per_complete_run(self):
        first=self.prepared()
        def launched(package,job):
            plan=json.loads((package/'input_plan.json').read_text());run=Path(plan['run_directory'])
            (run/'state').mkdir(parents=True);(run/'state/head_job.txt').write_text(str(job)+'\n')
        launched(first,101)
        config=json.loads((self.project/'file_manifests/lab_qc/config.json').read_text())
        def sample(item,lab,config,folder):
            return dict(sample_id=lab['sample_submitter_id'],run_id=lab['run_id']),dict(reference='same.fa',fai='same.fa.fai',dict='same.dict')
        calls=[]
        with patch.object(qc,'sample_manifest',side_effect=sample):
            result=qc.tick(config,lambda *a,**k:self.fail('cannot submit a partial run'),only_run='RUN_T')
            self.assertEqual(result[0]['status'],'waiting_for_bams')
            self.assertIn('CASE_B',result[0]['reason'])
            second=self.prepared('CASE_B');launched(second,102)
            config=json.loads((self.project/'file_manifests/lab_qc/config.json').read_text())
            def submit(command,**kwargs):calls.append(command);return SimpleNamespace(stdout='999\n')
            qc.tick(config,submit,only_run='RUN_T');qc.tick(config,submit,only_run='RUN_T')
        self.assertEqual(len(calls),1)
        self.assertIn('WIGITS_SITE_ENV=/synthetic/shared-site.env',str(calls[0]))
        state=qc.read(Path(config['state_directory'])/'RUN_T.run.json')
        manifest=qc.read(Path(state['execution_root'])/'manifest.json')
        self.assertEqual(len(manifest['samples']),2)


if __name__=='__main__':unittest.main()
