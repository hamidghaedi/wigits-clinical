#!/usr/bin/env python3
"""Measure whether unmodified SAGE can supply context for pilot DS PASS alleles.

This is an isolated diagnostic, not a new selected callset or a forced genotype.
Use native SAGE 5.0.2 specific_var_file (CSV) with the original reference,
resources, reads and scientific defaults. It may fail to discover DS-only alleles.
Do not use this restricted pass to estimate sample-wide contamination or phasing.
"""
import argparse
import csv
from datetime import datetime,timezone
import json
import hashlib
import os
from pathlib import Path
import subprocess
from deepsomatic_adapter import read_rows,audit


def stage_redux_sidecars(config,out):
    """Native SAGE uses one jitter_bqr_dir for both samples, even with separate BAM dirs."""
    target=out/'redux_sidecars';target.mkdir()
    sources={}
    for role in ('normal','tumor'):
        sample=config[role+'_id'];parent=Path(config[role+'_bam']).parent
        for suffix in ('bqr.tsv','jitter_params.tsv','ms_table.tsv.gz'):
            source=parent/f'{sample}.redux.{suffix}'
            if not source.is_file():raise FileNotFoundError(source)
            (target/source.name).symlink_to(source.resolve())
            sources[source.name]=str(source.resolve())
    return target,sources


def select_probe(vcf):
    rows=[row for _,row in read_rows(vcf) if row[6]=='PASS']
    if any(',' in row[4] for row in rows):raise ValueError('Normalize alleles first')
    if not 0<len(rows)<=1000:raise ValueError('Bounded pilot requires 1–1000 PASS alleles')
    return rows


def evidence_command(native, jar):
    """Use actual SAGE evidence code with retention overlay; never synthesize context.

    Restricted-scope TINC is not a sample-wide measurement and must not rewrite
    this evidence export. The production/control caller command is unaffected.
    """
    cmd=list(native)
    i=cmd.index('sage')
    cmd[i:i+2]=['java','-Xmx24g','-Ddeepsomatic.evidence_only=true',
                 '-cp',str(jar),'com.hartwig.hmftools.sage.SageApplication']
    return [x for x in cmd if x not in ('-run_tinc','-write_fit_variants')]+['-log_candidates']


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',type=Path,required=True);p.add_argument('--adapted',type=Path,required=True)
    p.add_argument('--evidence-jar',type=Path,help='Isolated SAGE evidence-retention build, not the caller image')
    p.add_argument('--out',type=Path,required=True);a=p.parse_args();c=json.loads(a.config.read_text())
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False);rows=select_probe(a.adapted)
    sidecar_dir,sidecars=stage_redux_sidecars(c,out)
    with (out/'requested.csv').open('w') as f:
        w=csv.writer(f);w.writerow(['Chromosome','Position','Ref','Alt'])
        w.writerows([row[0],row[1],row[3],row[4]] for row in rows)
    r=Path(c['resources']);cpus=int(os.environ['SLURM_CPUS_PER_TASK'])
    cmd=['apptainer','exec','--cleanenv','--bind','/global',c['sage_image'],'sage','-Xmx24g',
         '-reference',c['normal_id'],'-reference_bam',c['normal_bam'],'-ref_sample_count','1',
         '-tumor',c['tumor_id'],'-tumor_bam',c['tumor_bam'],
         '-jitter_bqr_dir',str(sidecar_dir)+'/',
         '-ref_genome',c['reference'],'-ref_genome_version','38',
         '-hotspots',str(r/'dna/variants/KnownHotspots.somatic.38.vcf.gz'),
         '-driver_gene_panel',str(r/'common/DriverGenePanel.38.tsv'),
         '-high_confidence_bed',c['high_confidence_bed'],'-ensembl_data_dir',str(r/'common/ensembl_data'),
         '-sequencing_type','ILLUMINA','-include_mt','-run_tinc','-write_fit_variants',
         '-pon_file',str(r/'dna/variants/hmf_wgs_sage_pon_1000.38.tsv.gz'),
         '-gnomad_freq_dir',str(r/'dna/variants/gnomad'),'-threads',str(cpus),'-log_level','DEBUG',
         '-specific_var_file',str(out/'requested.csv'),'-output_vcf',str(out/'probe.sage.somatic.vcf.gz')]
    evidence_build=None
    if a.evidence_jar:
        jar=a.evidence_jar.resolve()
        evidence_build=dict(path=str(jar),sha256=hashlib.sha256(jar.read_bytes()).hexdigest())
        cmd=evidence_command(cmd,jar)
    record=dict(version='0.2.0',job=os.environ['SLURM_JOB_ID'],requested_pass_alleles=len(rows),command=cmd,redux_sidecars=sidecars,
                evidence_build=evidence_build,
                started_utc=datetime.now(timezone.utc).isoformat(),
                limits='Restricted native evidence test; no full-sample TINC/phasing conclusions, no substituted outputs')
    path=out/'execution.json';path.write_text(json.dumps(record,indent=2)+'\n')
    with (out/'sage.log').open('w') as log:
        run=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
    record.update(exit_code=run.returncode,finished_utc=datetime.now(timezone.utc).isoformat())
    path.write_text(json.dumps(record,indent=2)+'\n');run.check_returncode()
    audit(a.adapted,out/'probe.sage.somatic.vcf.gz',out/'probe_audit.tsv',c['tumor_id'],c['normal_id'])


if __name__=='__main__':main()
