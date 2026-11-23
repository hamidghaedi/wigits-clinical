import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow'))
from prepare_context_seeds import exact_indel_context, prepare


class ContextSeedTests(unittest.TestCase):
    def test_context_must_encode_exact_allele_not_repeat_substring(self):
        self.assertTrue(exact_indel_context('chr1:13:A:AT', '10-3-CC-CA-TGG-4M1I2M'))
        self.assertFalse(exact_indel_context('chr1:13:A:AT', '10-3-CC-CA-TGG-7M'))
        self.assertFalse(exact_indel_context('chr1:14:A:AT', '10-3-CC-CA-TGG-4M1I2M'))
        self.assertFalse(exact_indel_context('chr1:13:A:AG', '10-3-CC-CA-TGG-4M1I2M'))
        self.assertTrue(exact_indel_context('chr1:13:ATT:A', '10-3-CC-CA-GG-4M2D2M'))

    def test_haplotype_only_observation_is_never_promoted(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            header = '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tT\n'
            source = p / 'source.vcf'
            source.write_text(header + 'chr1\t13\t.\tA\tAT\t20\tPASS\tTIER=LOW_CONFIDENCE\tGT\t0/1\n')
            obs = p / 'obs.tsv'
            fields = ['requested', 'observed', 'context', 'method', 'read_name', 'native_events']
            with obs.open('w') as f:
                w = csv.DictWriter(f, fieldnames=fields, delimiter='\t'); w.writeheader()
                for method, name in [('ALT_HAPLOTYPE_40', 'read1'), ('ALT_HAPLOTYPE_20', 'read2')]:
                    w.writerow(dict(requested='chr1:13:A:AT', observed='chr1:13:A:AT',
                                    context='10-3-CC-CA-TGG-7M', method=method,
                                    read_name=name, native_events=0))
            self.assertEqual(prepare(source, obs, source, p / 'out.vcf'), 0)

    def test_native_primary_resolves_tie_only_when_observed(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);source=p/'source.vcf';obs=p/'obs.tsv';primary=p/'primary.tsv'
            source.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tT\n'
                              'chr1\t13\t.\tA\tAT\t20\tPASS\tTIER=LOW_CONFIDENCE\tGT\t0/1\n')
            obs.write_text('requested\tobserved\tcontext\tmethod\tread_name\tnative_events\n'
                           'chr1:13:A:AT\tchr1:13:A:AT\t10-3-CC-CA-TGG-4M1I2M\tCIGAR\tr1\t1\n'
                           'chr1:13:A:AT\tchr1:13:A:AT\t10-3-CC-CA-TGA-4M1I2M\tCIGAR\tr2\t1\n')
            self.assertEqual(prepare(source,obs,source,p/'tied.vcf'),0)
            primary.write_text('requested\tcontext\trole\nchr1:13:A:AT\t10-3-CC-CA-TGG-4M1I2M\tPRIMARY\n')
            self.assertEqual(prepare(source,obs,source,p/'selected.vcf',primary),1)
            audit=json.loads((p/'selected.vcf.json').read_text())
            self.assertEqual(audit['chr1:13:A:AT']['selection'],'NATIVE_ALTCONTEXT')
            primary.write_text(primary.read_text().replace('TGG','TGC'))
            with self.assertRaisesRegex(ValueError,'lacks exact observed'):
                prepare(source,obs,source,p/'bad.vcf',primary)


if __name__ == '__main__':
    unittest.main()
