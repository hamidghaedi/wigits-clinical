import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

SPEC=importlib.util.spec_from_file_location('labqc',Path(__file__).parents[1]/'workflow/lab_qc_report.py')
lab=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(lab)


class LaboratoryQcTests(unittest.TestCase):
    def test_q30_is_base_weighted_and_mates_are_individual_reads(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths=[]
            for i,(bases,q30,reads) in enumerate([(100,90,2),(900,450,6)]):
                p=Path(tmp)/f'{i}.json';row=dict(total_bases=bases,q30_bases=q30,total_reads=reads)
                p.write_text(json.dumps({'summary':{'before_filtering':row,'after_filtering':row}}));paths.append(p)
            self.assertEqual(lab.raw_summary(paths),dict(total_reads=8,total_bases=1000,q30_bases=540,q30_pct=54))

    def test_qc_must_not_change_yield(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'x.json';p.write_text(json.dumps({'summary':{
                'before_filtering':dict(total_reads=2,total_bases=100,q30_bases=90),
                'after_filtering':dict(total_reads=1,total_bases=50,q30_bases=45)}}))
            with self.assertRaisesRegex(ValueError,'changed'):lab.raw_summary([p])

    def test_duplicate_primary_counts_exclude_supplementary(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'flags';p.write_text('120 + 0 in total (QC-passed reads + QC-failed reads)\n100 + 0 primary\n20 + 0 supplementary\n15 + 0 duplicates\n10 + 0 primary duplicates\n')
            self.assertEqual(lab.duplicate_summary(p)['bam_duplicate_pct'],10)

    def test_picard_table_stops_before_histogram(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'metric';p.write_text('## METRICS CLASS\tpicard.analysis.CollectWgsMetrics$WgsMetrics\nMEAN_COVERAGE\tPCT_30X\n40.5\t0.75\n\n## HISTOGRAM\n0\t1\n')
            self.assertEqual(lab.metric_rows(p),[dict(MEAN_COVERAGE='40.5',PCT_30X='0.75')])

    def test_membership_is_by_sequencing_run_not_pair(self):
        samples=[dict(run_id='T_RUN',sample_id='TEST-DNA',fastqs=[dict(label='S1_L001')]),
                 dict(run_id='N_RUN',sample_id='TEST-gDNA',fastqs=[dict(label='S2_L001')])]
        self.assertEqual([s['sample_id'] for s in lab.load_samples(dict(samples=samples),'T_RUN')],['TEST-DNA'])
        with self.assertRaises(ValueError):lab.load_samples(dict(samples=samples+samples),'T_RUN')


if __name__=='__main__':unittest.main()
