import argparse
import importlib.util
from pathlib import Path
import tempfile
import unittest

spec=importlib.util.spec_from_file_location('probe',Path(__file__).resolve().parents[1]/'workflow/probe_snp_support.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class SnpSupportTests(unittest.TestCase):
    def test_alt_absent_in_pileup_is_measured_zero(self):
        row=m.counts({'AD':'21,0','DP':'21'},'A','T',['A','<*>'])
        self.assertEqual((row['ref'],row['alt'],row['allele_fraction']),(21,0,0))
        self.assertIsNone(m.counts({'AD':'.','DP':'.'},'A','T',['A','T'])['alt'])

    def test_sample_order_and_missing_site_retained(self):
        targets=[dict(chrom='chr1',pos=1,ref='A',alt='T'),dict(chrom='chr1',pos=2,ref='G',alt='C')]
        text='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tN\tT\nchr1\t1\t.\tA\tT,<*>\t.\t.\t.\tAD:DP\t20,0,0:20\t1,19,0:20\n'
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'a.vcf';p.write_text(text);rows=m.parse_pileup(p,targets,'T','N')
        self.assertEqual(rows[0]['tumor']['alt'],19)
        self.assertEqual(rows[1]['status'],'NO_PILEUP_RECORD')
        self.assertIsNone(rows[1]['tumor']['dp'])
        summary=m.summarize(rows)
        self.assertEqual(summary['opposing_homozygous_like'],1)
        self.assertEqual(summary['informative'],1)

    def test_wrong_reference_or_sample_is_rejected(self):
        targets=[dict(chrom='chr1',pos=1,ref='G',alt='T')]
        text='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tN\tT\nchr1\t1\t.\tA\tT\t.\t.\t.\tAD:DP\t20,0:20\t1,19:20\n'
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'a.vcf';p.write_text(text)
            with self.assertRaises(ValueError):m.parse_pileup(p,targets,'T','N')
            with self.assertRaises(ValueError):m.parse_pileup(p,targets,'WRONG','N')

    def test_only_biallelic_snps_and_unique_positions(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'targets.vcf';p.write_text('chr1\t1\t.\tA\tT,G\n')
            with self.assertRaises(ValueError):m.read_targets(p)
            p.write_text('chr1\t1\t.\tA\tT\nchr1\t1\t.\tA\tG\n')
            with self.assertRaises(ValueError):m.read_targets(p)

    def test_command_is_restricted_to_target_regions_without_call(self):
        a=argparse.Namespace(container=Path('/synthetic/bcf.sif'),reference=Path('/synthetic/ref.fa'),tumor_bam=Path('/synthetic/T.bam'),normal_bam=Path('/synthetic/N.bam'))
        cmd=m.command(a,Path('/synthetic/100.bed'),Path('/synthetic/out.vcf'))
        self.assertIn('mpileup',cmd);self.assertNotIn('call',cmd)
        self.assertEqual(cmd[cmd.index('--regions-file')+1],'/synthetic/100.bed')
        self.assertEqual(cmd[cmd.index('--ff')+1],'3844')


if __name__=='__main__':unittest.main()
