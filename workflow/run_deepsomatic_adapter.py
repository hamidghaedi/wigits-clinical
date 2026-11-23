#!/usr/bin/env python3
"""Run measured count mapping and native reference annotation under Slurm.

No SAGE-compatible caller claim, RNA/NEO disabling or downstream submission.
PAVE/PURPLE consumption awaits the native read-context increment. Uses existing
SAGE 5.0.2 and bcftools images only; no installations or reference downloads.
"""
import argparse
from datetime import datetime,timezone
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from deepsomatic_adapter import verify_annotations,preserve_source_payload


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--source',type=Path,required=True);a=p.parse_args()
    c=json.loads(a.config.read_text());out=a.out.resolve();src=a.source.resolve()
    out.mkdir(parents=True,exist_ok=False)
    record=dict(version='0.1.3',job=os.environ.get('SLURM_JOB_ID'),started_utc=datetime.now(timezone.utc).isoformat(),
                config=c,steps=[],downstream_ready=False)
    def save(): (out/'execution.json').write_text(json.dumps(record,indent=2)+'\n')
    def run(cmd,name):
        start=time.monotonic();item=dict(name=name,command=list(map(str,cmd)));record['steps'].append(item);save()
        with (out/(name+'.log')).open('w') as log:
            r=subprocess.run(item['command'],stdout=log,stderr=subprocess.STDOUT)
        item.update(exit_code=r.returncode,seconds=time.monotonic()-start);save();r.check_returncode()
    def image(key):return ['apptainer','exec','--cleanenv','--bind','/global',c[key]]
    adapter=src/'deepsomatic_adapter.py'
    try:
        run([sys.executable,adapter,'counts','--vcf',c['vcf'],'--out',out/'paired.vcf',
             '--tumour',c['tumor_id'],'--normal',c['normal_id'],'--somatic-candidates'],'counts')
        run(image('bcftools_image')+['bcftools','norm','-m','-any','--old-rec-tag','DS_ORIGINAL_ALLELE',
                                  '-Ov','-o',out/'normalized.vcf',out/'paired.vcf'],'normalize')
        jar=subprocess.check_output(image('sage_image')+['bash','-c','d=$(dirname "$(readlink -f "$(command -v sage)")"); printf "%s/sage.jar" "$d"'],text=True).strip()
        java=image('sage_image')+['java','-Xmx6g','--class-path',jar]
        run(java+[src/'deepsomatic_native/ReferenceRepeatAnnotator.java','--self-test'],'native_repeat_selftest')
        res=Path(c['resources'])
        run(java+[src/'deepsomatic_native/NativeSageRegions.java',out/'normalized.vcf',out/'regions.vcf',
                  '-tumor',c['tumor_id'],'-tumor_bam',c['tumor_bam'],
                  '-reference',c['normal_id'],'-reference_bam',c['normal_bam'],
                  '-output_vcf',out/'unused_caller_output.vcf.gz','-ref_genome',c['reference'],'-ref_genome_version','38',
                  '-hotspots',res/'dna/variants/KnownHotspots.somatic.38.vcf.gz',
                  '-high_confidence_bed',c['high_confidence_bed'],
                  '-driver_gene_panel',res/'common/DriverGenePanel.38.tsv',
                  '-ensembl_data_dir',res/'common/ensembl_data'],'native_regions')
        run(java+[src/'deepsomatic_native/ReferenceRepeatAnnotator.java',out/'regions.vcf',out/'native_annotations.vcf',c['reference']],'native_reference_context')
        preserve_source_payload(out/'normalized.vcf',out/'native_annotations.vcf',out/'annotated.vcf')
        record['preservation']=verify_annotations(out/'normalized.vcf',out/'annotated.vcf');save()
        run(image('bcftools_image')+['bcftools','view','-Oz','-o',out/'annotated.vcf.gz',out/'annotated.vcf'],'compress')
        run(image('bcftools_image')+['bcftools','index','-t',out/'annotated.vcf.gz'],'index')
        run([sys.executable,adapter,'audit','--vcf',out/'annotated.vcf.gz','--sage',c['sage_vcf'],
             '--out',out/'native_context_audit.tsv','--tumour',c['tumor_id'],'--normal',c['normal_id']],'native_context_audit')
        record['exit_code']=0
    except Exception as e:
        record.update(exit_code=1,error=str(e));raise
    finally:
        record['finished_utc']=datetime.now(timezone.utc).isoformat();save()


if __name__=='__main__':main()
