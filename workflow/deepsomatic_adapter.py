#!/usr/bin/env python3
"""Experimental paired-count adapter; no selection or scientific filter changes.

Component 0.1.3. Map real DeepSomatic FORMAT/NAD,NDP to normal AD,DP;
normal GT stays missing and AF is explicitly derived from NAD/NDP. Tumour
genotype fields, QUAL, FILTER and INFO remain unchanged. Add AF from tumour
VAF, keeping VAF. With explicit --somatic-candidates exclude only RefCall and
GERMLINE from the derivative, retaining the original file and counting exclusions.
Every other FILTER class remains unchanged. Native bcftools projects multiallelic records; no allele
realignment is applied. NAF has an inconsistent Number=R declaration in Google
1.10 despite ALT-only values in observed records: retain its raw vector with
Number=. and an explicit original-allele description. It is not used for counts.

The native SAGE evidence audit never supplies fabricated RC_INFO/RC_CNT or LPS.
It records exact allele matches and missing context separately. This intermediate
VCF is not released as a feature-complete downstream sage_dir.
"""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re

VERSION='0.1.3'


def reader(path):
    return gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)


def info(text):
    return dict(x.split('=',1) if '=' in x else (x,True) for x in text.split(';') if x and x!='.')


def allele_key(row):
    return tuple(row[i] for i in (0,1,3,4))


def fractions(counts, depth):
    if depth in ('','.'): return '.'
    depth=int(depth)
    if depth<0:raise ValueError('Negative depth')
    if depth==0:return '.'
    return ','.join('.' if x=='.' else format(int(x)/depth,'.8g') for x in counts[1:])


def map_record(row):
    fields=row[8].split(':'); values=row[9].split(':')
    if len(fields)!=len(values):raise ValueError('Malformed FORMAT payload')
    original=dict(zip(fields,values));nalleles=len(row[4].split(','))+1
    for tag in ['AD','NAD']:
        value=original.get(tag,'.')
        if value!='.' and len(value.split(','))!=nalleles:
            raise ValueError(f'{tag} does not match allele count at {row[0]}:{row[1]}')
        if value!='.' and any(x!='.' and int(x)<0 for x in value.split(',')):
            raise ValueError('Negative allele count')
    for tag in ['DP','NDP']:
        if original.get(tag,'.')!='.' and int(original[tag])<0:raise ValueError('Negative depth')
    for tag in ['GT','AD','DP']:
        if tag not in original:raise ValueError('Missing tumour field '+tag)
    if 'AF' in original:raise ValueError('Source already has AF; specify mapping rather than overwrite')
    output_fields=fields+['AF']
    tumour_af=original.get('VAF','.')
    if tumour_af!='.' and len(tumour_af.split(','))!=nalleles-1:
        raise ValueError('VAF does not match ALT count')
    normal={k:'.' for k in output_fields}
    normal.update(GT='./.',AD=original.get('NAD','.'),DP=original.get('NDP','.'))
    normal['AF']=fractions(normal['AD'].split(','),normal['DP']) if normal['AD']!='.' else '.'
    output=row[:8]+[':'.join(output_fields),':'.join(normal[k] for k in output_fields),row[9]+':'+tumour_af]
    state='MEASURED' if normal['AD']!='.' and normal['DP']!='.' else 'MISSING_NORMAL_COUNTS'
    return output,state


def build(source,out,tumour,normal,candidates_only=False):
    if tumour==normal:raise ValueError('Tumour and normal IDs must differ')
    counts=Counter();filters=Counter();source_filters=Counter();header_seen=False
    with reader(source) as inp,open(out,'x') as dest:
        for line in inp:
            if line.startswith('##FORMAT=<ID=NAF,'):
                dest.write('##FORMAT=<ID=NAF,Number=.,Type=Float,Description="Original DeepSomatic NAF vector in original allele order before normalization; preserved verbatim; normal AF derives from NAD/NDP">\n')
            elif line.startswith('#CHROM'):
                columns=line.rstrip('\n').split('\t')
                if columns[9:]!=[tumour]:raise ValueError('Expected one named DeepSomatic tumour column')
                dest.write('##FORMAT=<ID=AF,Number=A,Type=Float,Description="Tumour: original VAF; normal: measured NAD alternate count divided by NDP, not a normal genotype call">\n')
                dest.write(f'##deepsomaticCountAdapter={VERSION}\n')
                dest.write('##deepsomaticCountAdapterStatus=experimental_not_downstream_ready\n')
                dest.write('\t'.join(columns[:9]+[normal,tumour])+'\n');header_seen=True
            elif line.startswith('#'):dest.write(line)
            else:
                if not header_seen:raise ValueError('Missing VCF header')
                row=line.rstrip('\n').split('\t')
                if len(row)!=10:raise ValueError('Expected one tumour sample')
                source_filters[row[6]]+=1
                if candidates_only and set(row[6].split(';')) & {'RefCall','GERMLINE'}:
                    continue
                mapped,state=map_record(row);dest.write('\t'.join(mapped)+'\n')
                counts[state]+=1;filters[row[6]]+=1
                if ',' in row[4]:counts['multiallelic_records']+=1
    summary=dict(version=VERSION,records=sum(filters.values()),normal_counts=dict(counts),filters=dict(filters),
                 source_records=sum(source_filters.values()),source_filters=dict(source_filters),
                 candidate_scope='exclude_RefCall_and_GERMLINE_only' if candidates_only else 'all_records',
                 source=str(Path(source).resolve()),output=str(Path(out).resolve()),
                 source_sha256=file_digest(source),
                 native_read_context='unresolved; no synthetic fields emitted',downstream_ready=False)
    Path(str(out)+'.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary


def file_digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        while chunk:=f.read(1024*1024):h.update(chunk)
    return h.hexdigest()


def read_rows(path):
    samples=[]
    with reader(path) as f:
        for line in f:
            if line.startswith('#CHROM'):samples=line.rstrip().split('\t')[9:]
            elif not line.startswith('#'):yield samples,line.rstrip('\n').split('\t')


def verify_annotations(before,after):
    """Annotations must not change allele/QUAL/FILTER/genotypes or drop records."""
    from itertools import zip_longest
    total=0
    for left,right in zip_longest(read_rows(before),read_rows(after)):
        if left is None or right is None:raise ValueError('Annotation changed record count')
        ls,lr=left;rs,rr=right
        same_qual=lr[5]==rr[5] or (lr[5]!='.' and rr[5]!='.' and float(lr[5])==float(rr[5]))
        if ls!=rs or lr[:5]!=rr[:5] or lr[6]!=rr[6] or not same_qual:
            raise ValueError('Annotation changed allele, QUAL or FILTER')
        for i in range(len(ls)):
            a=dict(zip(lr[8].split(':'),lr[9+i].split(':')))
            b=dict(zip(rr[8].split(':'),rr[9+i].split(':')))
            for k,v in a.items():
                if b.get(k,'.')!=v:raise ValueError('Annotation changed genotype field '+k)
        total+=1
    return dict(records=total,alleles_filters_qual_and_genotypes_preserved=True)


def preserve_source_payload(source,annotations,out):
    """Take only INFO annotations; retain exact source QUAL/FILTER/FORMAT values.

    HTSJDK may reformat QUAL or collapse missing vectors. These are not new
    measurements; use the normalized source payload verbatim instead.
    """
    from itertools import zip_longest
    with open(out,'x') as f:
        with reader(annotations) as annotated:
            for line in annotated:
                if not line.startswith('#'):break
                f.write(line)
        for left,right in zip_longest(read_rows(source),read_rows(annotations)):
            if left is None or right is None:raise ValueError('Annotation changed record count')
            ls,lr=left;rs,rr=right
            if ls!=rs or allele_key(lr)!=allele_key(rr):raise ValueError('Annotation changed allele/order')
            result=lr[:];result[7]=rr[7]
            f.write('\t'.join(result)+'\n')


def audit(adapted,sage,out,tumour,normal):
    targets=list(read_rows(adapted));keys={allele_key(row) for _,row in targets}
    if any(',' in key[-1] for key in keys):raise ValueError('Audit normalized biallelic records')
    matches={}
    for samples,row in read_rows(sage):
        key=allele_key(row)
        if key in keys:matches.setdefault(key,[]).append((samples,row))
    rows=[];counts=Counter()
    for samples,row in targets:
        key=allele_key(row)
        ai=info(row[7]);ds=dict(zip(row[8].split(':'),row[9+samples.index(tumour)].split(':')))
        rec=dict(chrom=key[0],pos=key[1],ref=key[2],alt=key[3],ds_filter=row[6],
                 exact_sage_match=key in matches,ds_tumour_AD=ds.get('AD','.'),
                 ds_tier=ai.get('TIER','.'),tier_agreement='NA',sage_filter='.',
                 rc_info=False,tumour_rc_cnt=False,normal_rc_cnt=False,lps=False,
                 sage_tumour_AD='.',sage_normal_AD='.',read_context_available=False,
                 native_record_count=len(matches.get(key,[])),native_records_json='[]')
        native_records=[]
        for ss,sr in matches.get(key,[]):
            si=info(sr[7]);formats=sr[8].split(':')
            st=dict(zip(formats,sr[9+ss.index(tumour)].split(':')))
            sn=dict(zip(formats,sr[9+ss.index(normal)].split(':')))
            native_records.append(dict(filter=sr[6],tier=si.get('TIER'),rc_info=si.get('RC_INFO'),
                                       tumour_rc_cnt=st.get('RC_CNT'),normal_rc_cnt=sn.get('RC_CNT'),
                                       lps=si.get('LPS'),tumour_AD=st.get('AD'),normal_AD=sn.get('AD')))
        rec['native_records_json']=json.dumps(native_records,separators=(',',':'))
        if len(native_records)==1:
            nr=native_records[0]
            rec.update(sage_filter=nr['filter'],rc_info=bool(nr['rc_info']),
                       tumour_rc_cnt=nr['tumour_rc_cnt'] not in (None,'.'),normal_rc_cnt=nr['normal_rc_cnt'] not in (None,'.'),
                       lps=bool(nr['lps']),sage_tumour_AD=nr['tumour_AD'],sage_normal_AD=nr['normal_AD'],
                       tier_agreement=str(nr['tier']==ai.get('TIER')))
            rec['read_context_available']=rec['rc_info'] and rec['tumour_rc_cnt'] and rec['normal_rc_cnt']
        if len(native_records)>1:
            counts['multiple_native_records']+=1
            if row[6]=='PASS':counts['ds_pass_multiple_native_records']+=1
        counts['total']+=1
        for flag in ['exact_sage_match','read_context_available','lps']:
            counts[flag]+=int(rec[flag])
            if row[6]=='PASS':counts['ds_pass_'+flag]+=int(rec[flag])
        if row[6]=='PASS':counts['ds_pass_total']+=1
        rows.append(rec)
    if rows:
        with open(out,'x') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    summary=dict(counts=counts,downstream_ready=False,
                 remaining=['Authentic context for DeepSomatic-only alleles','Phasing completeness and MNV handling',
                            'Actual SAGE_APPEND/NEO consumer tests','TINC diagnostic/recovery semantics'],
                 note='Exact SAGE context availability is not caller agreement or approval to copy classifications. No DS allele was removed.')
    Path(str(out)+'.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['counts','audit'])
    p.add_argument('--vcf',required=True);p.add_argument('--out',required=True)
    p.add_argument('--tumour',required=True);p.add_argument('--normal',required=True);p.add_argument('--sage')
    p.add_argument('--somatic-candidates',action='store_true',help='Explicit derivative excluding RefCall/GERMLINE; original VCF retained, other FILTERs unchanged')
    a=p.parse_args()
    if a.mode=='counts':result=build(a.vcf,a.out,a.tumour,a.normal,a.somatic_candidates)
    else:
        if not a.sage:p.error('audit needs --sage')
        result=audit(a.vcf,a.sage,a.out,a.tumour,a.normal)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
