#!/usr/bin/env python3
"""Bounded native-consumer tests; the complete adapter VCF is never filtered.

The probe explicitly tests measured PASS alleles and records unresolved alleles
separately. It cannot establish feature completeness when that list is nonempty.
PAVE and RNA append are real executions; NEO is its installed count reader,
not a full neoantigen prediction. Full PURPLE needs the full caller output.
"""
import argparse
from collections import Counter
from datetime import datetime,timezone
import json
import os
from pathlib import Path
import subprocess
from deepsomatic_adapter import read_rows,allele_key,info
from merge_sage_evidence import headers,merge,consumer_view


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['config','source','evidence','candidates','out']:
        p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--supplement',type=Path)
    a=p.parse_args();c=json.loads(a.config.read_text());out=a.out.resolve();out.mkdir(exist_ok=False)
    src=Path(__file__).resolve().parent
    rec=dict(started_utc=datetime.now(timezone.utc).isoformat(),job=os.environ['SLURM_JOB_ID'],steps=[],full_downstream_ready=False)
    def save(): (out/'execution.json').write_text(json.dumps(rec,indent=2)+'\n')
    def run(cmd,name):
        rec['steps'].append(dict(name=name,command=list(map(str,cmd))));save()
        with (out/(name+'.log')).open('w') as f:
            result=subprocess.run(list(map(str,cmd)),stdout=f,stderr=subprocess.STDOUT)
        rec['steps'][-1]['exit_code']=result.returncode;save();result.check_returncode()
    def image(key):return ['apptainer','exec','--cleanenv','--bind','/global',c[key]]
    def jar(key,tool):
        return subprocess.check_output(image(key)+['bash','-c',
            'd=$(dirname "$(readlink -f "$(command -v '+tool+')")"); printf "%s/'+tool+'.jar" "$d"'],text=True).strip()
    try:
        merged=out/'all_with_evidence.vcf';rec['merge']=merge(a.source,a.evidence,merged,a.candidates,a.supplement);save()
        pass_rows=[(s,r) for s,r in read_rows(merged) if r[6]=='PASS']
        unresolved=[allele_key(r) for s,r in pass_rows if info(r[7])['DS_EVIDENCE_STATUS']!='MEASURED']
        rec['unresolved_pass_alleles']=unresolved
        rows=[r for s,r in pass_rows if info(r[7])['DS_EVIDENCE_STATUS']=='MEASURED']
        rec['tested_pass_alleles']=len(rows);save()
        raw_probe=out/'measured_pass_source.vcf'
        with raw_probe.open('w') as f:
            f.write('\n'.join(headers(merged))+'\n')
            for row in rows:f.write('\t'.join(row)+'\n')
        probe=out/'measured_pass_probe.vcf'
        consumer_view(raw_probe,probe)
        sage_jar=jar('sage_image','sage')
        run(image('sage_image')+['java','-Xmx4g','--class-path',sage_jar,
            src/'deepsomatic_native/NativeContextCheck.java',probe,c['reference']],'native_context_roundtrip')
        run(image('bcftools_image')+['bcftools','view','-Oz','-o',out/'probe.vcf.gz',probe],'compress')
        run(image('bcftools_image')+['bcftools','index','-t',out/'probe.vcf.gz'],'index')
        run(image('sage_image')+['sage','-Xmx24g','com.hartwig.hmftools.sage.append.SageAppendApplication',
            '-input_vcf',out/'probe.vcf.gz','-max_read_depth','100000',
            '-reference',c['rna_id'],'-reference_bam',c['rna_bam'],
            '-ref_genome',c['reference'],'-ref_genome_version','38','-sequencing_type','ILLUMINA',
            '-write_frag_lengths','-skip_bqr','-skip_msi_jitter','-threads',os.environ['SLURM_CPUS_PER_TASK'],
            '-output_vcf',out/'rna.append.vcf.gz'],'native_rna_append')
        appended=list(read_rows(out/'rna.append.vcf.gz'))
        if [allele_key(r) for s,r in appended]!=[allele_key(r) for r in rows]:
            raise ValueError('Native RNA append changed tested allele list')
        rec['rna_alleles_retained']=len(appended);save()
        consumer_view(out/'rna.append.vcf.gz',out/'neo.context_depth.vcf',neo=True)
        run(image('neo_image')+['java','-Xmx4g','--class-path',jar('neo_image','neo'),
            src/'deepsomatic_native/NativeNeoCounts.java',out/'neo.context_depth.vcf',c['tumor_id'],c['rna_id'],out/'neo.counts.tsv'],
            'native_neo_reader')
        r=Path(c['resources'])
        run(image('pave_image')+['pave','-Xmx24g','-sample',c['tumor_id'],'-input_vcf',out/'probe.vcf.gz',
            '-ref_genome',c['reference'],'-ref_genome_version','38',
            '-pon_file',r/'dna/variants/hmf_wgs_sage_pon_1000.38.tsv.gz','-gnomad_freq_dir',r/'dna/variants/gnomad',
            '-clinvar_vcf',r/'dna/variants/clinvar.38.vcf.gz','-driver_gene_panel',r/'common/DriverGenePanel.38.tsv',
            '-mappability_bed',r/'dna/variants/mappability_150.38.bed.gz','-ensembl_data_dir',r/'common/ensembl_data',
            '-sequencing_type','ILLUMINA','-threads',os.environ['SLURM_CPUS_PER_TASK'],
            '-output_vcf',out/'pave.vcf.gz'],'native_pave')
        paved=list(read_rows(out/'pave.vcf.gz'))
        rec['pave_records']=len(paved)
        if Counter(allele_key(r) for s,r in paved)!=Counter(allele_key(r) for r in rows):
            raise ValueError('Native PAVE changed the tested allele multiset')
        rec['pave_alleles_retained']=True
        rec['exit_code']=0
    except Exception as e:
        rec.update(exit_code=1,error=str(e));raise
    finally:
        rec['finished_utc']=datetime.now(timezone.utc).isoformat();save()


if __name__=='__main__':main()
