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

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('intake',ROOT/'workflow/bam_intake.py')
intake=importlib.util.module_from_spec(spec);spec.loader.exec_module(intake)


class HeaderTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.ref=self.root/'reference.fa';self.ref.write_text('>chr1\nAAAA\n')
        self.dictionary=self.root/'reference.dict'
        self.dictionary.write_text('@SQ\tSN:chr1\tLN:4\tM5:098890dde069e9abad63f19a0d9e1f32\n')
        self.h=self.dictionary.read_text()+'@RG\tID:rg\tSM:synthetic\tPL:ILLUMINA\n'
    def tearDown(self):self.tmp.cleanup()
    def test_matching_md5(self):
        self.assertIn('M5',intake.check_header(self.h,self.dictionary,'synthetic',self.ref))
    def test_reference_mismatch(self):
        for old,new in [('LN:4','LN:5'),('098890dde069e9abad63f19a0d9e1f32','wrong')]:
            with self.subTest(field=old),self.assertRaises(ValueError):
                intake.check_header(self.h.replace(old,new),self.dictionary,'synthetic',self.ref)
    def test_missing_mixed_sample_and_platform(self):
        for h in [self.h.split('@RG')[0],self.h.replace('SM:synthetic','SM:other'),
                  self.h+'@RG\tID:other\tSM:other\tPL:ILLUMINA\n',self.h.replace('PL:ILLUMINA','PL:ULTIMA')]:
            with self.subTest(header=h),self.assertRaises(ValueError):
                intake.check_header(h,self.dictionary,'synthetic',self.ref)
    def test_reference_origin_required_without_md5(self):
        h=self.h.replace('\tM5:098890dde069e9abad63f19a0d9e1f32','')
        with self.assertRaises(ValueError):intake.check_header(h,self.dictionary,'synthetic',self.ref)
        self.assertIn('source PG',intake.check_header(h+'@PG\tID:align\tCL:align --ref /released/reference.fa\n',self.dictionary,'synthetic',self.ref))
    def test_raw_bam_launch_rejected(self):
        sheet=self.root/'sheet.csv';sheet.write_text('sample_id,filetype,sequence_type,filepath\nsynthetic,bam,dna,/raw.bam\n')
        with self.assertRaises(ValueError):intake.verify_launch(sheet,None)
    def test_failed_intake_cannot_publish(self):
        sheet=self.root/'sheet.csv';sheet.write_text('sample_id,filetype,sequence_type,filepath\nsynthetic,bam,dna,/raw.bam\n')
        (self.root/'synthetic').mkdir();(self.root/'synthetic/receipt.json').write_text(json.dumps(dict(status='FAILED')))
        with self.assertRaises(ValueError):intake.publish(sheet,self.root,self.root/'ready.csv')
        self.assertFalse((self.root/'ready.csv').exists())
    def test_failed_intake_prevents_analytical_launcher(self):
        import shutil
        here=self.root/'workflow';(here/'conf').mkdir(parents=True)
        (here/'conf/frontenac_paths.env').write_text(f'RUN_ROOT={self.root}/runs\nSLURM_ACCOUNT_CPU=synthetic\n')
        (here/'conf/bam_intake.env').write_text('INTAKE_REFERENCE_FASTA=ref\nINTAKE_REFERENCE_DICT=dict\nINTAKE_REFERENCE_SHA256=bad\nINTAKE_GATK_COMMAND=gatk\n')
        (here/'bam_intake.py').write_text('raise SystemExit(1)\n')
        marker=self.root/'unexpected_analysis'
        (here/'build_input_resources.py').write_text(f'open({str(marker)!r},"w").write("bad")\n')
        (here/'run_oncoanalyser_frontenac.sh').write_text(f'#!/bin/bash\ntouch {marker}\n')
        sheet=self.root/'source.csv';sheet.write_text('sample_id,filetype\nsynthetic,bam\n')
        result=subprocess.run(['bash',str(ROOT/'workflow/run_bam_oncoanalyser_frontenac.sbatch'),'test',str(sheet)],
                              env={**os.environ,'WIGITS_WORKFLOW_DIR':str(here)},capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assertFalse(marker.exists())


@unittest.skipUnless(os.environ.get('INTAKE_TEST_GATK_COMMAND'),'Slurm container integration required')
class RealToolTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.ref=self.root/'reference.fa';self.ref.write_text('>chr1\n'+'A'*2000+'\n')
        self.digest=hashlib.sha256(self.ref.read_bytes()).hexdigest()
        self.dictionary=self.root/'reference.dict'
        self.dictionary.write_text('@SQ\tSN:chr1\tLN:2000\tM5:'+hashlib.md5(b'A'*2000).hexdigest()+'\n')
    def tearDown(self):self.tmp.cleanup()
    def make(self,name='source',mc=None,hard=False,bad_rg=False):
        rows=['@HD\tVN:1.6\tSO:unsorted',self.dictionary.read_text().strip(),
              '@RG\tID:rg\tSM:synthetic\tPL:ILLUMINA']
        def rec(q,flag,pos,cig,mate,tag=None):
            chrom='chr1' if pos else '*'
            f=[q,str(flag),chrom,str(pos),'60' if pos else '0',cig,
               '=' if pos else '*',str(mate),'0','A'*20,'I'*20,
               'RG:Z:'+('other' if bad_rg else 'rg'),'OQ:Z:'+'G'*20,'BC:Z:ACGT']
            if tag is not None:f.append('MC:Z:'+tag)
            rows.append('\t'.join(f))
        rec('pair',99,101,'20M',151,mc);rec('pair',147,151,'10M10S',101,mc)
        rec('pair',99+2048,301,'10H10M' if hard else '10S10M',151,mc)
        if hard:
            # Hard-clipped record carries only the unclipped sequence, valid SAM.
            rows[-1]=rows[-1].replace('A'*20,'A'*10).replace('I'*20,'I'*10).replace('G'*20,'G'*10)
        rec('pair',99+256,501,'20M',151,mc)
        rec('duplicate',99+1024,601,'20M',651,mc);rec('duplicate',147+1024,651,'20M',601,mc)
        rec('qcfail',99+512,701,'20M',751,mc);rec('qcfail',147+512,751,'20M',701,mc)
        rec('discordant',65,801,'20M',851,mc);rec('discordant',129,851,'20M',801,mc)
        rec('unmapped',77,0,'*',0);rec('unmapped',141,0,'*',0)
        rows=rows[:3]+list(reversed(rows[3:]))
        sam=self.root/f'{name}.sam';sam.write_text('\n'.join(rows)+'\n')
        bam=self.root/f'{name}.bam'
        subprocess.run(['samtools','view','--no-PG','-b','-o',str(bam),str(sam)],check=True)
        return bam
    def args(self,bam):
        return argparse.Namespace(bam=bam,sample='synthetic',outdir=self.root/'out',threads=8,
            reference=self.ref,dictionary=self.dictionary,reference_sha256=self.digest,
            gatk_command=os.environ['INTAKE_TEST_GATK_COMMAND'])
    def prepare_checked(self,bam):
        try:intake.prepare(self.args(bam))
        except Exception:
            log=self.root/'out/commands.log'
            if log.exists():print(log.read_text()[-3500:],flush=True)
            raise
    def test_missing_or_incorrect_mc_fixed_and_all_records_preserved(self):
        bam=self.make(mc='1M')
        before=hashlib.sha256(bam.read_bytes()).hexdigest()
        self.prepare_checked(bam)
        data=json.loads((self.root/'out/receipt.json').read_text())
        self.assertEqual(data['status'],'PASS')
        self.assertEqual(data['before']['audit']['records'],12)
        self.assertEqual(data['after']['audit']['records'],12)
        self.assertEqual(data['after']['audit']['missing_mc'],0)
        self.assertEqual(hashlib.sha256(bam.read_bytes()).hexdigest(),before)
        out=data['output']['path']
        lines=subprocess.check_output(['samtools','view',out],text=True).splitlines()
        left=[x for x in lines if x.startswith('pair\t99\t')][0]
        self.assertIn('MC:Z:10M10S',left)
        self.assertTrue(Path(out+'.bai').exists())
        # Passed intake is consumed by the only supported analytical launch path.
        sheet=self.root/'raw.csv';sheet.write_text('sample_id,sample_type,sequence_type,filetype,filepath\nsynthetic,tumor,dna,bam,'+str(bam)+'\n')
        receipts=self.root/'receipts';(receipts/'synthetic').mkdir(parents=True)
        (receipts/'synthetic/receipt.json').write_text(json.dumps(data))
        intake.publish(sheet,receipts,self.root/'ready.csv')
        intake.verify_launch(self.root/'ready.csv',receipts,self.ref,self.digest)
        with self.assertRaises(ValueError):intake.verify_launch(sheet,receipts,self.ref,self.digest)
        with self.assertRaises(ValueError):intake.verify_launch(self.root/'ready.csv',receipts,self.ref,'different-reference')
        with Path(out).open('ab') as h:h.write(b'changed')
        with self.assertRaises(ValueError):intake.verify_launch(self.root/'ready.csv',receipts,self.ref,self.digest)
    def test_missing_mc_and_broken_source_index(self):
        bam=self.make();Path(str(bam)+'.bai').write_bytes(b'bad index')
        self.prepare_checked(bam)
        self.assertEqual(json.loads((self.root/'out/receipt.json').read_text())['status'],'PASS')
    def test_valid_prepared_bam_retains_all_record_content(self):
        bam=self.make();self.prepare_checked(bam)
        first=json.loads((self.root/'out/receipt.json').read_text())
        good=Path(first['output']['path']);good_digest=hashlib.sha256(good.read_bytes()).hexdigest()
        args=self.args(good);args.outdir=self.root/'second'
        intake.prepare(args)
        self.assertEqual(hashlib.sha256(good.read_bytes()).hexdigest(),good_digest)
        second=json.loads((args.outdir/'receipt.json').read_text())
        command=['samtools','checksum','-F','0','-b','0xfff','-c','-P','-C','-M','-t','*','-T']
        old=subprocess.check_output(command+[str(good)],text=True)
        new=subprocess.check_output(command+[second['output']['path']],text=True)
        self.assertEqual(sorted(x for x in old.splitlines() if not x.startswith('#')),
                         sorted(x for x in new.splitlines() if not x.startswith('#')))
    def test_malformed_mc_rejected(self):
        bam=self.make(mc='badCIGAR')
        with self.assertRaises(ValueError):intake.prepare(self.args(bam))
        self.assertFalse((self.root/'out/synthetic.oncoanalyser.bam').exists())
    def test_hard_clipped_supplementary_rejected(self):
        bam=self.make(hard=True)
        with self.assertRaises(ValueError):intake.prepare(self.args(bam))
        self.assertFalse((self.root/'out/synthetic.oncoanalyser.bam').exists())
    def test_missing_record_read_group_rejected(self):
        bam=self.make(bad_rg=True)
        with self.assertRaises(ValueError):intake.prepare(self.args(bam))
    def test_truncated_bam_rejected(self):
        bam=self.make();bam.write_bytes(bam.read_bytes()[:-40])
        with self.assertRaises(subprocess.CalledProcessError):intake.prepare(self.args(bam))
    def test_reference_dictionary_really_matches_fasta(self):
        bam=self.make()
        self.dictionary.write_text(self.dictionary.read_text().replace('LN:2000','LN:2001'))
        with self.assertRaises(ValueError):intake.prepare(self.args(bam))
    def test_conservation_rejects_changed_quality_or_lost_record(self):
        for mutation in ['quality','original_quality','sequence','record_loss','alignment']:
            with self.subTest(mutation=mutation):
                bam=self.make(name=mutation)
                args=self.args(bam);args.outdir=self.root/mutation
                original_scan=intake.scan
                def changed_scan(path,sample,prefix,*rest,**kwargs):
                    if prefix.name=='prepared':
                        text=subprocess.check_output(['samtools','view','-h',str(path)],text=True)
                        lines=text.splitlines();idx=next(i for i,l in enumerate(lines) if not l.startswith('@'))
                        f=lines[idx].split('\t')
                        if mutation=='record_loss':lines.pop(idx)
                        else:
                            if mutation=='quality':f[10]='H'+f[10][1:]
                            elif mutation=='original_quality':f=[x.replace('OQ:Z:G','OQ:Z:H',1) for x in f]
                            elif mutation=='sequence':f[9]='C'+f[9][1:]
                            elif mutation=='alignment':f[4]='20'
                            lines[idx]='\t'.join(f)
                        temp=path.with_suffix('.changed.bam')
                        subprocess.run(['samtools','view','--no-PG','-b','-o',str(temp),'-'],input='\n'.join(lines)+'\n',text=True,check=True)
                        os.replace(temp,path)
                    return original_scan(path,sample,prefix,*rest,**kwargs)
                with patch.object(intake,'scan',changed_scan),self.assertRaises(ValueError):intake.prepare(args)
                self.assertEqual(json.loads((args.outdir/'receipt.json').read_text())['status'],'FAILED')


if __name__=='__main__':unittest.main()
