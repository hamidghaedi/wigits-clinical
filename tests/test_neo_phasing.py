"""Source-backed phase translation, including unsafe independent LPS collisions."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow'))
import neo_phasing as phase

HEADER = '''##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=LPS,Number=.,Type=Integer,Description="Local Phase Set">
##INFO=<ID=LPS_RC,Number=.,Type=Integer,Description="Local Phase Set Read Count">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOUR
'''


def row(pos, lps='7', counts='4', gt='0/1', filt='PASS', chrom='chr1', extra=''):
    info = f'LPS={lps};LPS_RC={counts}' if lps else '.'
    return f'{chrom}\t{pos}\t.\tA\tT\t50\t{filt}\t{info}{extra}\tGT:AD\t0/0:20,0\t{gt}:10,9\n'


class PhasingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def export(self, somatic, germline='', source=None, anchors=None):
        paths = []
        for name, content in [('somatic', somatic), ('germline', germline), ('source', source if source is not None else somatic)]:
            p = self.root / (name + '.vcf')
            p.write_text(HEADER + content)
            paths.append(str(p))
        output, evidence = self.root/'context.vcf', self.root/'evidence.json'
        counts = phase.export_context(paths[0], paths[1], paths[2], 'TUMOUR',
                     anchors or {phase.key(f) for l in somatic.splitlines() for f in [l.split('\t')] if f[0]=='chr1' and abs(int(f[1])-100)<=100 and f[6]=='PASS'}, output, evidence)
        records = {f[1]: f for line in output.read_text().splitlines() if not line.startswith('#')
                   for f in [line.split('\t')]}
        return records, json.loads(evidence.read_text()), counts

    def test_supported_cis_anchors_and_measured_tags_are_retained(self):
        rows, evidence, counts = self.export(row(100) + row(103))
        self.assertEqual(set(rows), {'100','103'})
        self.assertEqual(rows['100'][7], 'LPS=7;LPS_RC=4')
        self.assertEqual(rows['103'][8:], ['GT:AD:HP','0/1:10,9:SAGE_chr1_7'])
        self.assertEqual(counts['anchors'], 2)
        self.assertEqual(evidence['records'][0]['phase']['read_counts'], [4])

    def test_germline_lps_number_alone_cannot_establish_phase(self):
        rows, evidence, _ = self.export(row(100), row(103))
        self.assertNotIn('HP', rows['103'][8].split(':'))
        self.assertEqual(evidence['records'][1]['status'], 'NO_SOMATIC_PHASE_RECORD')

    def test_exact_germline_allele_can_use_tumour_somatic_phase_source(self):
        rows, evidence, _ = self.export(row(100), row(103, lps='99'), row(100)+row(103, filt='LOW_VAF'))
        self.assertTrue(rows['103'][8].endswith(':HP'))
        self.assertEqual(rows['103'][7], 'LPS=99;LPS_RC=4')
        self.assertEqual(evidence['records'][1]['phase']['filter'], 'LOW_VAF')
        self.assertEqual(rows['103'][6], 'PASS')  # original germline FILTER, never relabelled

    def test_overlapping_sets_are_unresolved_even_on_unexported_record(self):
        rows, evidence, _ = self.export(row(100)+row(103), source=row(100)+row(103)+row(500,'7,8','2,3'))
        self.assertTrue(all('HP' not in f[8].split(':') for f in rows.values()))
        self.assertTrue(all(r['status']=='AMBIGUOUS_OR_ABSENT_LPS' for r in evidence['records']))

    def test_missing_or_zero_read_support_does_not_invent_hp(self):
        for support in ['.', '0', '3,4']:
            rows, _, _ = self.export(row(100, counts=support))
            self.assertNotIn('HP', rows['100'][8].split(':'))

    def test_no_phase_stays_missing_and_different_groups_stay_different(self):
        rows, _, _ = self.export(row(100)+row(103, lps='8')+row(106,lps=''))
        self.assertNotEqual(rows['100'][-1], rows['103'][-1])
        self.assertNotIn('HP', rows['106'][8].split(':'))

    def test_native_sage_missing_gt_uses_measured_support_and_purple_genotype(self):
        rows, evidence, _ = self.export(row(100),source=row(100,gt='./.'))
        self.assertEqual(rows['100'][8:],['GT:AD:HP','0/1:10,9:SAGE_chr1_7'])
        self.assertEqual(evidence['records'][0]['phase']['gt'],'./.')
        self.assertEqual(evidence['records'][0]['phase']['ad'],[10,9])

    def test_no_alt_reads_or_conflicting_source_gt_cannot_link(self):
        for source in [row(100,gt='./.').replace('10,9','19,0'),row(100,gt='0/0')]:
            rows, _, _ = self.export(row(100),source=source)
            self.assertNotIn('HP', rows['100'][8].split(':'))

    def test_chromosomes_are_separate_namespaces_and_proximity_is_bounded(self):
        rows, _, _ = self.export(row(100)+row(103,chrom='chr2')+row(300))
        self.assertEqual(set(rows), {'100'})

    def test_nonselected_somatic_is_not_misrepresented_as_germline(self):
        rows, evidence, _ = self.export(row(100)+row(103),anchors={('chr1','100','A','T')})
        self.assertEqual(set(rows),{'100'})
        self.assertEqual(evidence['omitted_context'][0]['status'],'NONSELECTED_SOMATIC_CONTEXT_NOT_EXPORTED')

    def test_duplicate_phase_source_is_unresolved(self):
        rows, evidence, _ = self.export(row(100), source=row(100)+row(100,'8'))
        self.assertNotIn('HP',rows['100'][8].split(':'))
        self.assertEqual(evidence['records'][0]['status'],'DUPLICATE_PHASE_SOURCE')

    def test_homozygote_genotype_and_filtered_calls_are_not_rewritten(self):
        rows, _, _ = self.export(row(100)+row(106,filt='LOW_VAF'),row(103,gt='1/1'))
        self.assertEqual(set(rows),{'100','103'})
        self.assertEqual(rows['103'][8:],['GT:AD','1/1:10,9'])

    def test_missing_tumour_sample_does_not_use_normal_genotype(self):
        (self.root/'somatic.vcf').write_text(HEADER+row(100))
        germ = HEADER.replace('\tNORMAL\tTUMOUR','\tNORMAL') + row(103).rsplit('\t',1)[0]+'\n'
        (self.root/'germline.vcf').write_text(germ)
        out, ledger=self.root/'out.vcf', self.root/'ledger.json'
        phase.export_context(str(self.root/'somatic.vcf'), str(self.root/'germline.vcf'),
            str(self.root/'somatic.vcf'),'TUMOUR',{('chr1','100','A','T')},out,ledger)
        line=[l for l in out.read_text().splitlines() if l.startswith('chr1\t103')][0]
        self.assertEqual(line.split('\t')[-1],'.:.')
        self.assertEqual(json.loads(ledger.read_text())['missing_tumour_genotypes'],1)

    def test_germline_tumour_column_is_resolved_by_name_in_reversed_sample_order(self):
        somatic=self.root/'somatic.vcf';somatic.write_text(HEADER+row(100))
        fields=row(103,gt='1/1').rstrip('\n').split('\t')
        germline=self.root/'germline.vcf'
        germline.write_text(HEADER.replace('\tNORMAL\tTUMOUR','\tTUMOUR\tNORMAL')+
            '\t'.join(fields[:9]+[fields[10],fields[9]])+'\n')
        output=self.root/'context.vcf'
        phase.export_context(str(somatic),str(germline),str(somatic),'TUMOUR',
            {('chr1','100','A','T')},output,self.root/'evidence.json')
        line=next(l for l in output.read_text().splitlines() if l.startswith('chr1\t103'))
        self.assertEqual(line.split('\t')[-1],'1/1:10,9')

    def test_existing_hp_is_preserved_as_evidence_not_trusted_across_callers(self):
        s = row(100).replace('GT:AD','GT:AD:HP').replace('0/0:20,0','0/0:20,0:.').replace('0/1:10,9','0/1:10,9:old')
        rows, evidence, _ = self.export(s)
        self.assertIn('SW_ORIGINAL_HP', rows['100'][8])
        self.assertEqual(evidence['records'][0]['original_hp'],'old')

    def test_usage_records_actual_converter_rows_not_available_hp_counts(self):
        self.export(row(100))
        table=self.root/'prediction/MHC_Class_I/subject.proximal_variants.tsv'
        table.parent.mkdir(parents=True)
        table.write_text('main_somatic_variant\ttype\n1.GENE.TR.MT.1A/T\tsomatic\n')
        usage=phase.summarise_usage(self.root/'prediction',self.root/'evidence.json',self.root/'usage.json')
        self.assertEqual(usage['converter_tables'][0]['rows'],1)
        self.assertEqual(usage['converter_tables'][0]['types'],{'somatic':1})
        self.assertEqual(usage['source_phase_counts']['SUPPORTED_SOMATIC_LPS'],1)


if __name__ == '__main__':
    unittest.main()
