#!/usr/bin/env python3
"""Indexed read-only support at explicit SNPs, not a caller or identity verdict.

The supplied SNPCheck VCF fixes exact loci/REF/ALT. bcftools mpileup queries
only those one-base BED regions on the supplied indexed BAMs and reference;
no bcftools call, full BAM scan, new variant selection or report modification.
Operational settings: MQ>=50, BQ>=30, BAQ disabled, depth cap 10000 per BAM,
exclude UNMAP/SECONDARY/QCFAIL/DUP/SUPPLEMENTARY; retain bcftools' default
mate-overlap handling. These counts are independent QC observations, not
identical to AMBER or SAGE fragment definitions and not unique molecules.
Coarse AF bands (<=.1 homoREF, >=.9 homoALT, .2-.8 het) at DP>=10 are only
exploratory comparisons; they are not calibrated genotype/identity calls.
"""
import argparse
import collections
import csv
import datetime
import gzip
import hashlib
import json
from pathlib import Path
import subprocess

VERSION = '0.1.0'


def read_targets(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    targets = []
    seen = set()
    with opener(path, 'rt') as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            row = line.rstrip().split('\t')
            target = dict(chrom=row[0], pos=int(row[1]), ref=row[3], alt=row[4])
            if any(target[k] not in 'ACGT' or len(target[k]) != 1 for k in ('ref','alt')):
                raise ValueError('Only explicit biallelic SNP targets are supported')
            key = (target['chrom'], target['pos'])
            if key in seen:
                raise ValueError('Duplicate SNP position')
            seen.add(key)
            targets.append(target)
    if not targets:
        raise ValueError('No SNP targets')
    return targets


def counts(sample, ref, alt, alleles):
    if sample.get('DP', '.') == '.' or sample.get('AD', '.') == '.':
        return dict(ref=None, alt=None, dp=None, allele_fraction=None, other=None)
    dp = int(sample['DP'])
    ad = [int(n) for n in sample['AD'].split(',')]
    if len(ad) != len(alleles):
        raise ValueError('AD cardinality differs from pileup alleles')
    ref_count = ad[0]
    alt_count = ad[alleles.index(alt)] if alt in alleles else 0
    return dict(ref=ref_count, alt=alt_count, dp=dp,
                allele_fraction=alt_count/(ref_count+alt_count) if ref_count+alt_count else None,
                other=sum(ad)-ref_count-alt_count)


def parse_pileup(path, targets, tumor_id, normal_id):
    target_by_position = {(t['chrom'],t['pos']):t for t in targets}
    observations = {}
    with Path(path).open() as stream:
        for line in stream:
            if line.startswith('##'):
                continue
            row = line.rstrip().split('\t')
            if line.startswith('#CHROM'):
                columns = row
                if set(columns[9:]) != {tumor_id, normal_id}:
                    raise ValueError('Pileup sample IDs differ from supplied pair')
                continue
            key = (row[0],int(row[1]))
            if key not in target_by_position:
                raise ValueError('Pileup returned an unrequested locus')
            target = target_by_position[key]
            if len(row[3]) != 1:
                continue  # Ignore an additional indel record at the same locus.
            if row[3] != target['ref']:
                raise ValueError('Pileup reference differs from target REF')
            if key in observations:
                raise ValueError('Multiple SNP pileup records at one target')
            alleles = [row[3]] + row[4].split(',')
            fmt = row[8].split(':')
            samples = {s:dict(zip(fmt,row[columns.index(s)].split(':'))) for s in (tumor_id,normal_id)}
            observations[key] = {role:counts(samples[s],target['ref'],target['alt'],alleles)
                                 for role,s in [('tumor',tumor_id),('normal',normal_id)]}
    results = []
    for target in targets:
        pair = observations.get((target['chrom'],target['pos']))
        missing = dict(ref=None,alt=None,dp=None,allele_fraction=None,other=None)
        results.append(dict(target, status='MEASURED' if pair else 'NO_PILEUP_RECORD',
                            tumor=pair['tumor'] if pair else missing.copy(),
                            normal=pair['normal'] if pair else missing.copy()))
    return results


def band(support):
    af = support['allele_fraction']
    if support['dp'] is None or support['dp'] < 10 or af is None:
        return 'UNINFORMATIVE'
    if af <= .1:
        return 'REF_HOMOZYGOUS_LIKE'
    if af >= .9:
        return 'ALT_HOMOZYGOUS_LIKE'
    if .2 <= af <= .8:
        return 'HETEROZYGOUS_LIKE'
    return 'UNCERTAIN'


def summarize(rows):
    pairs = collections.Counter((band(r['tumor']),band(r['normal'])) for r in rows)
    informative = [(a,b,n) for (a,b),n in pairs.items()
                   if a not in ('UNINFORMATIVE','UNCERTAIN') and b not in ('UNINFORMATIVE','UNCERTAIN')]
    opposite = {('REF_HOMOZYGOUS_LIKE','ALT_HOMOZYGOUS_LIKE'),('ALT_HOMOZYGOUS_LIKE','REF_HOMOZYGOUS_LIKE')}
    return dict(targets=len(rows), measured=sum(r['status']=='MEASURED' for r in rows),
                informative=sum(n for _,_,n in informative),
                coarse_concordant=sum(n for a,b,n in informative if a==b),
                opposing_homozygous_like=sum(n for a,b,n in informative if (a,b) in opposite),
                categories=[dict(tumor=a,normal=b,count=n) for (a,b),n in sorted(pairs.items())],
                disposition='EXPLORATORY_ONLY_NO_CALIBRATED_IDENTITY_VERDICT')


def command(args, bed, output):
    return ['apptainer','exec','--cleanenv','--bind','/global',str(args.container),
            'bcftools','mpileup','--fasta-ref',str(args.reference),'--regions-file',str(bed),
            '--annotate','FORMAT/AD,FORMAT/DP','--min-MQ','50','--min-BQ','30',
            '--no-BAQ','--max-depth','10000','--ff','3844',
            '--output-type','v','--output',str(output),str(args.tumor_bam),str(args.normal_bam)]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for field in ('targets','reference','tumor-bam','normal-bam','container','out'):
        p.add_argument('--'+field, type=Path, required=True)
    p.add_argument('--tumor-id',required=True)
    p.add_argument('--normal-id',required=True)
    args=p.parse_args()
    targets=read_targets(args.targets)
    args.out.mkdir(parents=True,exist_ok=True)
    bed=args.out/'targets.bed'
    bed.write_text(''.join(f"{t['chrom']}\t{t['pos']-1}\t{t['pos']}\n" for t in targets))
    output=args.out/'support.mpileup.vcf'
    cmd=command(args,bed,output)
    version=subprocess.check_output(['apptainer','exec','--cleanenv',str(args.container),'bcftools','--version'],text=True)
    with (args.out/'mpileup.stderr').open('w') as err:
        subprocess.run(cmd,check=True,stderr=err)
    rows=parse_pileup(output,targets,args.tumor_id,args.normal_id)
    result=dict(version=VERSION,created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                command=cmd,bcftools_version=version,target_source=str(args.targets),
                target_sha256=hashlib.sha256(args.targets.read_bytes()).hexdigest(),
                inputs={k:str(v) for k,v in vars(args).items()},summary=summarize(rows),sites=rows)
    (args.out/'support.json').write_text(json.dumps(result,indent=2)+'\n')
    with (args.out/'support.tsv').open('w',newline='') as out:
        fields=['chrom','pos','ref','alt','status']+[r+'_'+k for r in ('tumor','normal') for k in ('ref','alt','dp','allele_fraction','other')]
        w=csv.DictWriter(out,fieldnames=fields,delimiter='\t');w.writeheader()
        for r in rows:
            item={k:r[k] for k in fields if k in r}
            for role in ('tumor','normal'):
                item.update({role+'_'+k:v for k,v in r[role].items()})
            w.writerow(item)
    print(json.dumps(result['summary']))


if __name__=='__main__':main()
