"""Development tests only; these are not stages in patient execution."""
import argparse
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import yaml

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('preparation', ROOT/'workflow/bam_intake.py')
preparation = importlib.util.module_from_spec(spec)
spec.loader.exec_module(preparation)


class ConfigurationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.config = self.root/'case.yaml'
        self.args = argparse.Namespace(config=self.config, outdir=self.root/'prepared',
                                       prepared_samplesheet=self.root/'samples.csv',
                                       account='example', gatk_command='gatk')

    def tearDown(self):
        self.tmp.cleanup()

    def sample(self, repair=False):
        return dict(group_id='pair', subject_id='example', sample_id='example-DNA',
                    sample_type='tumor', sequence_type='dna', filetype='bam',
                    info='cancer_type:Carcinoid', filepath=str(self.root/'source.bam'), prepare_mates=repair)

    def test_as_is_input_is_written_without_receipt_or_audit(self):
        sample = self.sample()
        self.config.write_text(yaml.safe_dump(dict(samples=[sample])))
        with patch.object(preparation.subprocess, 'check_output') as commands:
            preparation.batch(self.args)
        commands.assert_not_called()
        with self.args.prepared_samplesheet.open() as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(rows[0]['filepath'], sample['filepath'])
        self.assertEqual(rows[0]['info'], sample['info'])
        self.assertNotIn('prepare_mates', rows[0])
        self.assertFalse((self.args.outdir/'example-DNA').exists())

    def test_repair_job_and_completed_normal_are_mapped_to_native_csv(self):
        tumour = self.sample(True)
        Path(tumour['filepath']).write_bytes(b'BAM')
        normal = {**self.sample(), 'sample_id':'example-gDNA', 'sample_type':'normal',
                  'filepath':str(self.root/'already-prepared.bam'), 'info':''}
        self.config.write_text(yaml.safe_dump(dict(samples=[tumour, normal])))
        with patch.object(preparation.subprocess, 'check_output', side_effect=['123\n', 'COMPLETED|0:0\n']) as commands:
            preparation.batch(self.args)
        submit = commands.call_args_list[0].args[0]
        self.assertIn('--cpus-per-task=64', submit)
        self.assertIn('--mem=192G', submit)
        with self.args.prepared_samplesheet.open() as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(rows[0]['filepath'], str(self.args.outdir/'example-DNA/example-DNA.oncoanalyser.bam'))
        self.assertEqual(rows[1]['filepath'], normal['filepath'])

    def test_failed_tool_job_does_not_launch_analysis(self):
        sample = self.sample(True)
        Path(sample['filepath']).write_bytes(b'BAM')
        self.config.write_text(yaml.safe_dump(dict(samples=[sample])))
        with patch.object(preparation.subprocess, 'check_output', side_effect=['123\n', 'FAILED|1:0\n']):
            with self.assertRaises(RuntimeError):
                preparation.batch(self.args)
        self.assertFalse(self.args.prepared_samplesheet.exists())


@unittest.skipUnless(os.environ.get('INTAKE_TEST_GATK_COMMAND'), 'Real tools run under Slurm')
class StandardToolTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def source(self):
        sam = self.root/'source.sam'
        lines = ['@HD\tVN:1.6\tSO:unsorted', '@SQ\tSN:chr1\tLN:2000',
                 '@RG\tID:rg\tSM:example\tPL:ILLUMINA']
        for q,flag,pos,cigar,mate in [('pair',99,101,'20M',151), ('pair',147,151,'10M10S',101),
                                      ('pair',99+2048,301,'10S10M',151),
                                      ('duplicate',99+1024,601,'20M',651),
                                      ('duplicate',147+1024,651,'20M',601)]:
            lines.append('\t'.join([q,str(flag),'chr1',str(pos),'60',cigar,'=',str(mate),'0',
                                    'A'*20,'I'*20,'RG:Z:rg','OQ:Z:'+'G'*20]))
        sam.write_text('\n'.join(lines)+'\n')
        bam = self.root/'source.bam'
        subprocess.run(['samtools','view','-b','-o',str(bam),str(sam)],check=True)
        return bam

    def args(self, source):
        return argparse.Namespace(bam=source,sample='example',outdir=self.root/'out',threads=8,
                                  gatk_command=os.environ['INTAKE_TEST_GATK_COMMAND'])

    def test_picard_supplies_mate_tags_and_preserves_reads_and_qualities(self):
        source = self.source()
        original = hashlib.sha256(source.read_bytes()).hexdigest()
        preparation.prepare(self.args(source))
        output = self.root/'out/example.oncoanalyser.bam'
        rows = [r.split('\t') for r in subprocess.check_output(['samtools','view',str(output)],text=True).splitlines()]
        self.assertEqual(len(rows),5)
        self.assertTrue(all(any(v.startswith('MC:Z:') for v in row[11:]) for row in rows))
        before = [r.split('\t') for r in subprocess.check_output(['samtools','view',str(source)],text=True).splitlines()]
        keep = lambda r:(r[0],r[2],r[3],r[4],r[5],r[9],r[10],*[x for x in r[11:] if x.startswith(('RG:','OQ:'))])
        self.assertEqual(sorted(map(keep,before)),sorted(map(keep,rows)))
        self.assertEqual(hashlib.sha256(source.read_bytes()).hexdigest(), original)
        self.assertTrue(Path(str(output)+'.bai').exists())
        self.assertEqual(json.loads((self.root/'out/preparation.json').read_text())['status'],'PREPARED')
        self.assertFalse((self.root/'out/receipt.json').exists())

    def test_picard_failure_stops_preparation(self):
        source = self.source()
        args = self.args(source)
        args.gatk_command = 'false'
        with self.assertRaises(RuntimeError):
            preparation.prepare(args)
        self.assertEqual(json.loads((args.outdir/'preparation.json').read_text())['status'],'FAILED')
        self.assertFalse((args.outdir/'example.oncoanalyser.bam').exists())


if __name__ == '__main__':
    unittest.main()
