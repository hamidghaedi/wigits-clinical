import gzip
import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('ds_test',ROOT/'workflow/test_deepsomatic.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


class DeepSomaticSeamTests(unittest.TestCase):
    def config(self):
        return dict(model_type='FFPE_WGS',image='/test/google-1.10.0.sif',reference='/test/native.fa',
                    tumor_bam='/test/t.bam',normal_bam='/test/n.bam',tumor_id='T',normal_id='N',
                    regions=['chr17:1-10000','chr19:1-10000'])

    def test_official_model_reference_sample_bindings_and_bounded_workers(self):
        cmd=module.command(self.config(),Path('/test/out'),16,'0')
        for arg in ['run_deepsomatic','--model_type=FFPE_WGS','--ref=/test/native.fa',
                    '--reads_tumor=/test/t.bam','--reads_normal=/test/n.bam',
                    '--sample_name_tumor=T','--sample_name_normal=N','--num_shards=16',
                    '--regions=chr17:1-10000 chr19:1-10000',
                    '--postprocess_variants_extra_args=cpus=16,num_partitions=48',
                    'CUDA_VISIBLE_DEVICES=0']:
            self.assertIn(arg,cmd)
        self.assertNotIn('pbrun',cmd)

    def test_missing_scope_or_gpu_is_not_silently_changed(self):
        for config,cuda in [(self.config(),''),({**self.config(),'regions':[]},'0')]:
            with self.assertRaises(ValueError):module.command(config,Path('/test/out'),16,cuda)

    def test_empty_callset_and_filters_are_retained_as_observations(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'calls.vcf.gz'
            header='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tT\n'
            with gzip.open(p,'wt') as h:h.write(header)
            self.assertEqual(module.summarize_vcf(p)['records'],0)
            with gzip.open(p,'wt') as h:
                h.write(header+'chr17\t10\t.\tA\tT\t10\tGERMLINE\t.\tGT:AD:DP:VAF\t0/0:9,1:10:0.1\n')
            result=module.summarize_vcf(p)
            self.assertEqual(result['filters'],{'GERMLINE':1})
            self.assertEqual(result['samples'],['T'])
