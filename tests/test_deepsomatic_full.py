import importlib.util
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

spec=importlib.util.spec_from_file_location('full',Path(__file__).resolve().parents[1]/'workflow/deepsomatic_full.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class FullCallerTests(unittest.TestCase):
    def test_full_scope_and_native_defaults(self):
        c=dict(scope='whole_genome',model_type='FFPE_WGS',image='google.sif',reference='r.fa',
            tumor_bam='t.bam',normal_bam='n.bam',tumor_id='T',normal_id='N')
        cmd=m.dry_command(c,Path('/out'))
        self.assertFalse(any('--regions' in x for x in cmd))
        self.assertIn('--num_shards=64',cmd)
        self.assertIn('--call_variants_extra_args=batch_size=256',cmd)
        self.assertIn('--postprocess_variants_extra_args=cpus=64,num_partitions=192',cmd)
        with self.assertRaises(ValueError):m.dry_command({**c,'regions':['chr17']},Path('/out'))

    def test_native_command_split(self):
        text='\n'.join('***** Running the command:*****\n'+cmd+'\n' for cmd in
            ('time seq 0 63 | parallel make_examples_somatic --task {}',
             'time call_variants --examples a@64.gz','time postprocess_variants --cpus 64'))
        commands=m.parse_commands(text)
        self.assertEqual(list(commands),list(m.STAGES))
        with self.assertRaises(ValueError):m.parse_commands(text.replace('--task','--regions'))
        with self.assertRaises(ValueError):m.parse_commands('')

    def test_gpu_only_for_inference(self):
        with patch.dict(os.environ,{'CUDA_VISIBLE_DEVICES':'2'}):
            self.assertIn('CUDA_VISIBLE_DEVICES=2',m.container({'image':'g.sif'},True))
            self.assertNotIn('--nv',m.container({'image':'g.sif'}))

    def test_resource_growth(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'r.fa.fai').write_text('chr1\t3200000000\t0\t60\t61\n')
            (p/'t').touch();(p/'n').touch()
            c=dict(reference=str(p/'r.fa'),tumor_bam=str(p/'t'),normal_bam=str(p/'n'))
            r=m.resources(c)
            self.assertEqual(r['examples']['cpus'],64)
            self.assertGreaterEqual(r['examples']['hours'],48)
            self.assertGreaterEqual(r['inference']['hours'],60)
