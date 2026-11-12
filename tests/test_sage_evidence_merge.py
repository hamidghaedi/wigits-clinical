from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflow'))
from merge_sage_evidence import merge, consumer_view
from deepsomatic_adapter import read_rows, info


HEADER='''##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele counts">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tN\tT
'''
NATIVE=HEADER.replace('#CHROM', '''##sageVersion=5.0.2
##INFO=<ID=RC_INFO,Number=1,Type=String,Description="Native read context">
##FORMAT=<ID=RC_CNT,Number=6,Type=Integer,Description="Native fragment support">
#CHROM''')


class EvidenceMergeTests(unittest.TestCase):
    def setup_files(self,root):
        source=root/'ds.vcf';native=root/'native.vcf'
        source.write_text(HEADER+
            'chr17\t10\t.\tA\tC\t30\tPASS\t.\tGT:AD:DP\t./.:20,1:21\t0/1:25,15:40\n'+
            'chr17\t20\t.\tA\tG\t5\tNoCall\t.\tGT:AD:DP\t./.:20,0:20\t0/1:30,2:32\n')
        native.write_text(NATIVE+
            'chr17\t10\t.\tA\tC\t9\tminTumorQual\tRC_INFO=actual_read\tGT:AD:DP:RC_CNT\t./.:17,1:19:1,0,0,0,17,19\t./.:20,7:29:4,2,1,0,20,29\n')
        return source,native

    def test_native_filter_never_drops_ds_allele_and_counts_stay_separate(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);source,native=self.setup_files(p);out=p/'merged.vcf'
            result=merge(source,native,out)
            rows=list(read_rows(out));self.assertEqual(len(rows),2)
            self.assertEqual(rows[0][1][5:7],['30','PASS'])
            self.assertEqual(info(rows[0][1][7])['SAGE_FILTER'],'minTumorQual')
            self.assertEqual(info(rows[1][1][7])['DS_EVIDENCE_STATUS'],'MISSING_CONTEXT')
            fields=dict(zip(rows[0][1][8].split(':'),rows[0][1][10].split(':')))
            self.assertEqual(fields['AD'],'25,15');self.assertEqual(fields['SAGE_AD'],'20,7')
            self.assertEqual(fields['DP'],'40');self.assertEqual(fields['SAGE_DP'],'29')
            self.assertTrue(result['preservation']['alleles_filters_qual_and_genotypes_preserved'])
            self.assertEqual(consumer_view(out,p/'neo.vcf',neo=True),1)
            _,neo=next(read_rows(p/'neo.vcf'))
            fields=dict(zip(neo[8].split(':'),neo[10].split(':')))
            self.assertEqual(fields['DP'],'29');self.assertEqual(fields['DS_DP'],'40')

    def test_ambiguous_context_preserved_and_consumer_not_sent_guessed_context(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);source,native=self.setup_files(p)
            with native.open('a') as f:f.write(native.read_text().splitlines()[-1]+'\n')
            out=p/'merged.vcf';merge(source,native,out)
            rows=list(read_rows(out));self.assertEqual(len(rows),2)
            self.assertEqual(info(rows[0][1][7])['DS_EVIDENCE_STATUS'],'AMBIGUOUS_CONTEXT')
            with self.assertRaises(ValueError):consumer_view(out,p/'unsafe.vcf')

    def test_malformed_context_vector_is_not_replaced_with_ad(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);source,native=self.setup_files(p)
            native.write_text(native.read_text().replace('4,2,1,0,20,29','20,7'))
            with self.assertRaises(ValueError):merge(source,native,p/'merged.vcf')

    def test_multiple_contexts_use_native_primary_only_with_explicit_record(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);source,native=self.setup_files(p)
            text=native.read_text().replace('actual_read','1-1-L-A-R-3M')
            second=text.splitlines()[-1].replace('1-1-L-A-R-3M','1-1-L-C-R-3M')
            native.write_text(text+second+'\n')
            candidates=p/'candidates.tsv'
            candidates.write_text('Chromosome\tPosition\tRef\tAlt\tCoreBases\tCoreCigar\n'
                                  'chr17\t10\tA\tC\tC\t3M\nchr17\t10\tA\tC\tA\t3M\n')
            out=p/'merged.vcf';merge(source,native,out,candidates)
            _,row=next(read_rows(out));attrs=info(row[7])
            self.assertEqual(attrs['SAGE_CONTEXT_COUNT'],'2')
            self.assertEqual(attrs['RC_INFO'],'1-1-L-C-R-3M')
            self.assertEqual(attrs['DS_EVIDENCE_STATUS'],'MEASURED')

    def test_supplement_cannot_overwrite_discovery_or_fake_indel_context(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);source,native=self.setup_files(p)
            supplement=p/'append.vcf';supplement.write_text(native.read_text())
            with self.assertRaisesRegex(ValueError,'must not replace'):
                merge(source,native,p/'merged.vcf',supplement=supplement)
            supplement.write_text(native.read_text().replace('chr17\t10', 'chr17\t13')
                .replace('\tA\tC\t', '\tA\tAT\t')
                .replace('RC_INFO=actual_read', 'RC_INFO=10-3-CC-CA-TGG-7M;DS_SEED_METHOD=OBSERVED_BAM_READ')
                .replace('#CHROM', '##INFO=<ID=DS_SEED_METHOD,Number=1,Type=String,Description="Observed source">\n#CHROM'))
            with self.assertRaisesRegex(ValueError,'exact allele'):
                merge(source,native,p/'merged.vcf',supplement=supplement)

    def test_exact_cigar_supplement_preserves_caller_and_measured_counts(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);source,native=self.setup_files(p)
            with source.open('a') as f:
                f.write('chr17\t13\t.\tA\tAT\t11\tPASS\t.\tGT:AD:DP\t./.:22,0:22\t0/1:20,12:32\n')
            supplement=p/'append.vcf'
            supplement.write_text(native.read_text().replace('chr17\t10', 'chr17\t13')
                .replace('\tA\tC\t', '\tA\tAT\t')
                .replace('RC_INFO=actual_read', 'RC_INFO=10-3-CC-CA-TGG-4M1I2M;DS_SEED_METHOD=OBSERVED_BAM_READ')
                .replace('#CHROM', '##INFO=<ID=DS_SEED_METHOD,Number=1,Type=String,Description="Observed source">\n#CHROM'))
            result=merge(source,native,p/'merged.vcf',supplement=supplement)
            self.assertEqual(result['counts']['PASS_MEASURED'],2)
            self.assertTrue(result['preservation']['alleles_filters_qual_and_genotypes_preserved'])

    def test_consumer_uses_native_repeat_values_and_native_absence(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);source,native=self.setup_files(p)
            definition='##INFO=<ID=REP_C,Number=1,Type=Integer,Description="Repeat count">\n'
            source.write_text(source.read_text().replace('#CHROM',definition+'#CHROM')
                              .replace('\tPASS\t.\t','\tPASS\tREP_C=19\t'))
            native.write_text(native.read_text().replace('#CHROM',definition+'#CHROM')
                              .replace('RC_INFO=actual_read','RC_INFO=actual_read;REP_C=9'))
            merge(source,native,p/'merged.vcf');consumer_view(p/'merged.vcf',p/'consumer.vcf')
            _,row=next(read_rows(p/'consumer.vcf'));attrs=info(row[7])
            self.assertEqual(attrs['REP_C'],'9');self.assertEqual(attrs['DS_ORIGINAL_INFO_REP_C'],'19')
            consumer_view(p/'consumer.vcf',p/'twice.vcf')
            self.assertEqual(info(next(read_rows(p/'twice.vcf'))[1][7]),attrs)
            native.write_text(native.read_text().replace(';REP_C=9',''))
            merge(source,native,p/'absent.vcf');consumer_view(p/'absent.vcf',p/'absent.consumer.vcf')
            attrs=info(next(read_rows(p/'absent.consumer.vcf'))[1][7])
            self.assertNotIn('REP_C',attrs);self.assertEqual(attrs['DS_ORIGINAL_INFO_REP_C'],'19')
