import importlib.util
from pathlib import Path
import tempfile
import unittest

spec=importlib.util.spec_from_file_location('adapter',Path(__file__).resolve().parents[1]/'workflow/deepsomatic_adapter.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def row(alt='G',sample='0/1:20:14,6:0.3:18:17,1:0.0556',filt='PASS'):
    return ['chr1','100','.','A',alt,'25',filt,'.','GT:DP:AD:VAF:NDP:NAD:NAF',sample]


class AdapterTests(unittest.TestCase):
    def test_real_normal_support_and_original_tumour(self):
        original=row();out,state=m.map_record(original)
        self.assertEqual(out[:8],original[:8]);self.assertEqual(state,'MEASURED')
        self.assertEqual(out[-1],original[-1]+':0.3')
        normal=dict(zip(out[8].split(':'),out[9].split(':')))
        self.assertEqual(normal['AD'],'17,1');self.assertEqual(normal['DP'],'18')
        self.assertEqual(normal['GT'],'./.');self.assertAlmostEqual(float(normal['AF']),1/18)

    def test_missing_is_not_zero(self):
        out,state=m.map_record(row(sample='0/1:20:14,6:0.3:.:.:.'))
        self.assertEqual(state,'MISSING_NORMAL_COUNTS')
        self.assertEqual(out[9].split(':')[2],'.');self.assertEqual(out[9].split(':')[-1],'.')

    def test_measured_zero_depth_is_not_normal_genotype(self):
        out,state=m.map_record(row(sample='0/1:20:14,6:0.3:0:0,0:0'))
        self.assertEqual(state,'MEASURED');self.assertTrue(out[9].startswith('./.:0:0,0:'))
        self.assertEqual(out[9].split(':')[-1],'.')

    def test_multiallelic_counts_retained_until_native_projection(self):
        out,state=m.map_record(row(alt='G,T',sample='1/2:20:10,6,4:0.3,0.2:18:15,1,2:0.0556,0.1111'))
        self.assertEqual(out[9].split(':')[2],'15,1,2')
        self.assertEqual(out[-1].split(':')[2],'10,6,4')
        with self.assertRaises(ValueError):m.map_record(row(alt='G,T'))

    def test_filters_never_promoted(self):
        for filt in ['NoCall','GERMLINE','RefCall','.','PASS']:
            out,_=m.map_record(row(filt=filt));self.assertEqual(out[6],filt)

    def test_file_preservation_and_naf_declaration(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);source=p/'raw.vcf';out=p/'paired.vcf'
            text='##fileformat=VCFv4.2\n##FORMAT=<ID=NAF,Number=R,Type=Float,Description="normal">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tT\n'+'\t'.join(row())+'\n'
            source.write_text(text);result=m.build(source,out,'T','N')
            self.assertEqual(source.read_text(),text);self.assertEqual(result['records'],1)
            self.assertFalse(result['downstream_ready'])
            self.assertIn('ID=NAF,Number=.',out.read_text());self.assertIn('\tN\tT\n',out.read_text())
            with self.assertRaises(FileExistsError):m.build(source,out,'T','N')

    def test_explicit_somatic_view_does_not_promote_or_drop_nocalls(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);source=p/'raw.vcf';out=p/'paired.vcf'
            text='##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tT\n'
            text+=''.join('\t'.join(row(filt=f))+'\n' for f in ['RefCall','GERMLINE','NoCall','PASS'])
            source.write_text(text);result=m.build(source,out,'T','N',True)
            self.assertEqual(result['source_records'],4);self.assertEqual(result['records'],2)
            self.assertEqual(result['filters'],{'NoCall':1,'PASS':1});self.assertEqual(source.read_text(),text)

    def test_context_audit_keeps_ds_only_allele(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);a=p/'a.vcf';s=p/'s.vcf';o=p/'audit.tsv'
            head='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tN\tT\n'
            a.write_text(head+'chr1\t100\t.\tA\tG\t50\tPASS\tTIER=PANEL\tGT:AD:DP\t./.:18,0:18\t0/1:10,10:20\nchr1\t200\t.\tC\tT\t50\tPASS\tTIER=HIGH_CONFIDENCE\tGT:AD:DP\t./.:20,0:20\t0/1:5,15:20\n')
            s.write_text(head+'chr1\t100\t.\tA\tG\t50\tminTumorQual\tTIER=PANEL;RC_INFO=realcontext\tGT:AD:DP:RC_CNT\t0/0:18,0:18:0,0,0,0,18,18\t0/1:11,9:20:9,0,0,0,11,20\n')
            result=m.audit(a,s,o,'T','N')
            self.assertEqual(result['counts']['total'],2)
            self.assertEqual(result['counts']['read_context_available'],1)
            self.assertFalse(result['downstream_ready']);self.assertIn('minTumorQual',o.read_text())

    def test_multiple_sage_contexts_are_recorded_not_chosen(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);a=p/'a.vcf';s=p/'s.vcf';o=p/'audit.tsv'
            head='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tN\tT\n'
            data='chr1\t100\t.\tA\tG\t50\tPASS\tTIER=PANEL;RC_INFO=context\tGT:AD:DP:RC_CNT\t0/0:18,0:18:0,0,0,0,18,18\t0/1:11,9:20:9,0,0,0,11,20\n'
            a.write_text(head+data);s.write_text(head+data+data.replace('context','othercontext'))
            result=m.audit(a,s,o,'T','N')
            self.assertEqual(result['counts']['multiple_native_records'],1)
            self.assertEqual(result['counts']['read_context_available'],0)
            self.assertIn('othercontext',o.read_text())

    def test_annotation_preservation_compares_values_not_qual_formatting(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);a=p/'a.vcf';b=p/'b.vcf'
            head='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tT\n'
            data='chr1\t100\t.\tA\tG\t0.3\tPASS\t.\tGT:AD:DP\t0/1:10,10:20\n'
            a.write_text(head+data);b.write_text(head+data.replace('\t0.3\t','\t0.30\t'))
            self.assertEqual(m.verify_annotations(a,b)['records'],1)
            b.write_text(head+data.replace('10,10','11,9'))
            with self.assertRaises(ValueError):m.verify_annotations(a,b)

    def test_java_serialization_does_not_replace_source_payload(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);a=p/'source.vcf';b=p/'native.vcf';out=p/'final.vcf'
            head='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tT\n'
            line='chr1\t100\t.\tA\tG\t0.3\tNoCall\t.\tGT:AD:DP:NAD\t./.:10,10:20:.,.\n'
            a.write_text(head+line)
            b.write_text(head+line.replace('\t0.3\t','\t0.30\t').replace('\t.\tGT:','\tTIER=PANEL\tGT:').replace(':.,.',':.'))
            m.preserve_source_payload(a,b,out)
            self.assertEqual(m.verify_annotations(a,out)['records'],1)
            self.assertIn('TIER=PANEL',out.read_text());self.assertIn(':.,.',out.read_text())


if __name__=='__main__':unittest.main()
