"""Regressions for nonempty proximal records and mounted Python user sites."""
import collections
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow'))
import extension_tasks as tasks
import optional_extensions as planner


class OptionalPortabilityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_sort_preserves_annotations_genotypes_duplicates_and_header(self):
        path = self.root / 'proximal.vep.vcf'
        header = '##fileformat=VCFv4.2\n##contig=<ID=chr3,length=1000>\n##contig=<ID=chr6,length=1000>\n##contig=<ID=chr10,length=1000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOUR\n'
        def row(chrom, pos, alt='T'):
            return f'{chrom}\t{pos}\t.\tA\t{alt}\t52.7\tPASS\tCSQ=retained;LPS=7\tGT:AD\t0/1:20,8\n'
        rows = [row('chr3', 20), row('chr6', 80), row('chr10', 30), row('chr6', 15), row('chr6', 15, 'C,G'), row('chr6', 15)]
        path.write_text(header + ''.join(rows))
        tasks.sort_annotation_vcf(path)
        output = path.read_text().splitlines(keepends=True)
        self.assertEqual(''.join(l for l in output if l.startswith('#')), header)
        records = [l for l in output if not l.startswith('#')]
        self.assertEqual(collections.Counter(records), collections.Counter(rows))
        self.assertEqual([(l.split('\t')[0], int(l.split('\t')[1])) for l in records],
                         [('chr3',20),('chr6',15),('chr6',15),('chr6',15),('chr6',80),('chr10',30)])

    def test_empty_annotation_and_undeclared_contigs_remain_indexable_order(self):
        path = self.root / 'empty.vcf'
        header = '##fileformat=VCFv4.2\n#CHROM\tPOS\n'
        path.write_text(header)
        tasks.sort_annotation_vcf(path)
        self.assertEqual(path.read_text(), header)
        path.write_text(header + 'chrZ\t8\nchrA\t9\nchrZ\t2\n')
        tasks.sort_annotation_vcf(path)
        self.assertEqual(path.read_text(), header + 'chrA\t9\nchrZ\t2\nchrZ\t8\n')

    def test_container_explicit_isolation_overrides_injected_python_settings(self):
        out = self.root / 'output'
        out.mkdir()
        commands = tasks.Commands({'bind_roots':[str(self.root)], 'image':'/pinned/image.sif'}, out)
        with patch.dict(os.environ, {'APPTAINERENV_PYTHONNOUSERSITE':'0', 'APPTAINERENV_PYTHONPATH':'/host/packages'}), \
                patch.object(commands, 'run') as run:
            commands.container('image', 'pvacseq', '--help')
        argv = run.call_args.args[0]
        self.assertIn('PYTHONNOUSERSITE=1', argv)
        self.assertIn('PYTHONPATH=', argv)
        self.assertEqual(argv[-3:], ['/pinned/image.sif','pvacseq','--help'])

    def test_recovery_status_reads_completed_ch_without_copying_it(self):
        prior, current = self.root/'original', self.root/'recovery'
        ch = prior/'CASE/blood_ch/status.json'
        neo = current/'CASE/rna_neoantigens/status.json'
        planner.write_json(ch, {'status':'COMPLETED','reason':'existing CH'})
        planner.write_json(neo, {'status':'COMPLETED','reason':'recovered prediction'})
        plan = {'cases':[{'case_id':'CASE','blood_ch':{'status':'READY'},'rna_neoantigens':{'status':'READY'}}]}
        before = ch.read_bytes()
        result = planner.collect_status(plan, current, prior_results=prior)['cases'][0]
        self.assertEqual(result['blood_ch']['status'], 'COMPLETED')
        self.assertEqual(result['blood_ch']['status_file'], str(ch))
        self.assertEqual(result['rna_neoantigens']['status_file'], str(neo))
        self.assertEqual(ch.read_bytes(), before)
        self.assertFalse((current/'CASE/blood_ch').exists())

    def test_recovery_failure_is_not_hidden_by_prior_success(self):
        prior, current = self.root/'original', self.root/'recovery'
        planner.write_json(prior/'CASE/rna_neoantigens/status.json', {'status':'COMPLETED'})
        planner.write_json(current/'CASE/rna_neoantigens/status.json', {'status':'FAILED','reason':'new task failed'})
        plan = {'cases':[{'case_id':'CASE','blood_ch':{'status':'DISABLED'},'rna_neoantigens':{'status':'READY'}}]}
        result = planner.collect_status(plan, current, prior_results=prior)['cases'][0]
        self.assertEqual(result['rna_neoantigens']['status'], 'FAILED')

    def test_prediction_cannot_overwrite_completed_export(self):
        prepared = self.root/'export'
        prepared.mkdir()
        plan = {'cases':[{'case_id':'CASE'}]}
        with self.assertRaisesRegex(ValueError, 'separate from its preserved export'):
            tasks.execute(plan, 'CASE', 'neo_predict', prepared, 8, prepared)

    def test_wrong_export_identity_stops_prediction(self):
        prepared = self.root/'export'
        case = {'case_id':'CASE','tumour_sample':'T','normal_sample':'N','rna_sample':'R',
                'rna_neoantigens':{'status':'READY'}}
        planner.write_json(prepared/'task_metadata.json', {'status':'PREPARED'})
        planner.write_json(prepared/'CASE.interface_manifest.json',
                           {'case':'CASE','samples':{'tumour_dna':'OTHER','normal_dna':'N','tumour_rna':'R'}})
        with patch.object(tasks, 'neo_predict') as predict:
            with self.assertRaisesRegex(ValueError, 'Export source disagrees'):
                tasks.execute({'cases':[case],'resources':{}}, 'CASE', 'neo_predict', self.root/'out', 8, prepared)
        predict.assert_not_called()

    def test_prediction_reuses_valid_export_without_preparation_or_ch(self):
        prepared = self.root/'export'
        case = {'case_id':'CASE','tumour_sample':'T','normal_sample':'N','rna_sample':'R',
                'rna_neoantigens':{'status':'READY'}}
        planner.write_json(prepared/'task_metadata.json', {'status':'PREPARED'})
        manifest = prepared/'CASE.interface_manifest.json'
        planner.write_json(manifest, {'case':'CASE','samples':{'tumour_dna':'T','normal_dna':'N','tumour_rna':'R'}})
        before = manifest.read_bytes()
        with patch.object(tasks, 'neo_predict', return_value={'status':'COMPLETED'}) as predict, \
                patch.object(tasks, 'chip') as chip, patch.object(tasks, 'neo_prepare') as prepare:
            tasks.execute({'cases':[case],'resources':{}}, 'CASE', 'neo_predict', self.root/'out', 8, prepared)
        predict.assert_called_once()
        chip.assert_not_called()
        prepare.assert_not_called()
        self.assertEqual(manifest.read_bytes(), before)


if __name__ == '__main__':
    unittest.main()
