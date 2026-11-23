import importlib.util
from pathlib import Path
import tempfile
import gzip
import unittest

spec = importlib.util.spec_from_file_location('audit', Path(__file__).resolve().parents[1] / 'workflow/audit_native_somatic_calls.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class NativeAuditTests(unittest.TestCase):
    def test_missing_is_not_zero_and_native_af_is_not_ad_fraction(self):
        s = m.support({'AD': '90,5', 'DP': '100', 'AF': '.08'})
        self.assertEqual(s['native_af'], .08)
        self.assertEqual(s['alt_over_dp'], .05)
        self.assertIsNone(s['AFMQ'])
        self.assertIsNone(m.support({'AD': '1,2,3'})['alt'])

    def test_complements_and_tcw_are_symmetric(self):
        self.assertEqual(m.substitution('C', 'T', 'TCA'), ('C>T', True))
        self.assertEqual(m.substitution('G', 'A', 'TGA'), ('C>T', True))
        self.assertEqual(m.substitution('G', 'T', 'TGA'), ('C>A', False))
        self.assertEqual(m.substitution('C', 'T', 'ACA'), ('C>T', False))
        self.assertEqual(m.substitution('CT', 'C', 'TCA'), (None, None))

    def test_missing_metrics_do_not_trigger_flags(self):
        self.assertEqual(m.diagnostics(m.support({}), m.support({})), [])

    def test_exact_pass_only_and_selected_filtered_retained(self):
        header = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tN\tT\n'
        fmt = 'AD:DP:AF:FSB:RSB:AFBQ:AFMQ:MUC'
        normal = '98,0:100:0:0.5,.:0.5,.:.:.:.'
        tumor = '80,8:100:.1:0.5,0:0.5,1:30:60:7'
        rows = [f'chr1\t{pos}\t.\tC\tT\t100\t{filt}\tTNC=TCA;TIER=PANEL\t{fmt}\t{normal}\t{tumor}\n'
                for pos, filt in [(1, 'PASS'), (2, '.'), (3, 'minTumorQual')]]
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'a.vcf.gz'
            with gzip.open(p, 'wt') as f: f.write(header + ''.join(rows))
            result = m.scan_vcf(p, 'T', 'N', {'chr1:3:C:T'})
        self.assertEqual(result['counts']['all_records'], 3)
        self.assertEqual(result['counts']['pass_records'], 1)
        self.assertEqual(result['counts']['tcw_ct_or_cg'], 1)
        self.assertEqual(result['vaf_bins']['0.10-<0.20'], 1)
        self.assertEqual(result['diagnostic_counts']['fsb_alt_exclusive'], 1)
        self.assertEqual(result['selected_records']['chr1:3:C:T']['filter'], 'minTumorQual')


if __name__ == '__main__': unittest.main()
