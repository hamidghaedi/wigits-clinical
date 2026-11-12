#!/usr/bin/env python3
"""Join measured native evidence onto DeepSomatic alleles, never intersect callers.

DeepSomatic QUAL/FILTER/INFO and all existing FORMAT values remain authoritative.
Native collisions are stored as SAGE_* rather than overwriting caller values.
Native SAGE genotype fields that do not collide (RC_CNT etc.) retain their names
for native readers. The actual sageVersion header describes the context producer,
not the variant caller. Missing/ambiguous contexts remain explicit, not filled.
"""
import argparse
import csv
from collections import Counter, defaultdict
import json
from pathlib import Path
import re
from deepsomatic_adapter import reader, read_rows, allele_key, info, verify_annotations

VERSION='0.1.2'
DEFINITION=re.compile(r'##(INFO|FORMAT)=<ID=([^,]+),')
CONTEXT_TAGS=('TIER','TNC','REP_C','REP_S','MH')


def headers(path):
    out=[]
    with reader(path) as f:
        for line in f:
            if not line.startswith('#'):break
            out.append(line.rstrip('\n'))
    return out


def native_primary_contexts(path):
    """Native RefContextCache emits AltContext.candidate before secondCandidate.

    AltContext ranks full matches, core matches, then minimum event count via
    ReadContextCandidate.compareTo. Preserve this existing native choice and
    retain the entire multi-context evidence VCF; never add overlapping counts.
    """
    primary={}
    with reader(path) as f:
        for row in csv.DictReader(f,delimiter='\t'):
            key=tuple(row[k] for k in ('Chromosome','Position','Ref','Alt'))
            primary.setdefault(key,(row['CoreBases'],row['CoreCigar']))
    return primary


def merge(source, evidence, out, candidates=None, supplement=None):
    original=headers(source); native=headers(evidence)
    if supplement:
        primary_keys={allele_key(r) for _,r in read_rows(evidence)}
        extra_keys={allele_key(r) for _,r in read_rows(supplement)}
        if primary_keys & extra_keys:
            raise ValueError('Supplement must not replace native-discovery evidence')
        definitions={DEFINITION.match(x).groups():x for x in native if DEFINITION.match(x)}
        for line in headers(supplement):
            match=DEFINITION.match(line)
            if not match:continue
            key=match.groups()
            if key in definitions and definitions[key]!=line:
                raise ValueError('Conflicting evidence definition: '+str(key))
            if key not in definitions:
                native.append(line);definitions[key]=line
        if [x for x in headers(supplement) if x.startswith('##sageVersion=')]!=[x for x in native if x.startswith('##sageVersion=')]:
            raise ValueError('Supplement must use the same native evidence version')
    known={m.groups() for line in original if (m:=DEFINITION.match(line))}
    mappings={}
    added=[]
    for line in native:
        m=DEFINITION.match(line)
        if not m:continue
        kind,tag=m.groups()
        dest='SAGE_'+tag if (kind,tag) in known else tag
        if (kind,dest) in known and dest!=tag:
            raise ValueError('Input already contains native evidence '+dest)
        mappings[kind,tag]=dest
        added.append(line.replace('ID='+tag+',','ID='+dest+',',1))
    versions=[x for x in native if x.startswith('##sageVersion=')]
    if len(versions)!=1:raise ValueError('Actual native SAGE version header required')
    matches=defaultdict(list)
    for samples,row in read_rows(evidence):matches[allele_key(row)].append((samples,row))
    if supplement:
        for samples,row in read_rows(supplement):
            if info(row[7]).get('DS_SEED_METHOD')!='OBSERVED_BAM_READ':
                raise ValueError('Supplement requires explicit observed-context provenance')
            from prepare_context_seeds import exact_indel_context
            if not exact_indel_context(':'.join(allele_key(row)),str(info(row[7]).get('RC_INFO',''))):
                raise ValueError('Supplement context does not encode the exact allele')
            matches[allele_key(row)].append((samples,row))
    primary=native_primary_contexts(candidates) if candidates else {}
    counts=Counter()
    with open(out,'x') as f:
        f.write('\n'.join(x for x in original if not x.startswith('#CHROM'))+'\n')
        f.write('\n'.join(added+versions)+'\n')
        f.write(f'##deepsomaticEvidenceAdapter={VERSION}\n')
        f.write('##deepsomaticEvidenceScope=restricted_native_evidence_not_sample_wide_TINC\n')
        f.write('##INFO=<ID=DS_EVIDENCE_STATUS,Number=1,Type=String,Description="Measured native context availability; not caller confidence">\n')
        f.write('##INFO=<ID=SAGE_FILTER,Number=.,Type=String,Description="Native evidence producer filters; do not replace DeepSomatic FILTER">\n')
        f.write('##INFO=<ID=SAGE_QUAL,Number=1,Type=Float,Description="Native evidence quality; not DeepSomatic QUAL">\n')
        f.write('##INFO=<ID=SAGE_INFO_FIELDS,Number=.,Type=String,Description="INFO tags actually present in the native evidence record; preserves native optional-field absence">\n')
        f.write('##INFO=<ID=SAGE_CONTEXT_COUNT,Number=1,Type=Integer,Description="Number of native contexts; additional contexts retained in evidence VCF; primary is native AltContext candidate, counts are not summed">\n')
        f.write(next(x for x in original if x.startswith('#CHROM'))+'\n')
        for samples,row in read_rows(source):
            rs=matches.get(allele_key(row),[])
            context_count=len(rs)
            if len(rs)>1 and allele_key(row) in primary:
                core,cigar=primary[allele_key(row)]
                preferred=[]
                for ss,sr in rs:
                    parts=str(info(sr[7]).get('RC_INFO','')).split('-',5)
                    if len(parts)==6 and (parts[3],parts[5])==(core,cigar):preferred.append((ss,sr))
                if len(preferred)==1:
                    rs=preferred;counts['NATIVE_PRIMARY_MULTIPLE_CONTEXTS']+=1
            status='MEASURED' if len(rs)==1 else ('AMBIGUOUS_CONTEXT' if rs else 'MISSING_CONTEXT')
            values=info(row[7]); values['DS_EVIDENCE_STATUS']=status
            values['SAGE_CONTEXT_COUNT']=context_count
            result=row[:]
            if status=='MEASURED':
                ns,nr=rs[0]
                if set(ns)!=set(samples):raise ValueError('Native evidence sample identities differ')
                ni=info(nr[7])
                values['SAGE_INFO_FIELDS']=','.join(sorted(ni))
                nf=nr[8].split(':');sf=row[8].split(':')
                if not ni.get('RC_INFO'):raise ValueError('Matched native record lacks read context')
                for tag,value in ni.items():
                    values[mappings['INFO',tag]]=value
                values['SAGE_FILTER']=nr[6].replace(';',',')
                values['SAGE_QUAL']=nr[5]
                result[8]=':'.join(sf+[mappings['FORMAT',tag] for tag in nf])
                for i,sample in enumerate(samples):
                    nv=nr[9+ns.index(sample)].split(':')
                    if len(nv)!=len(nf):raise ValueError('Malformed native FORMAT')
                    fields=dict(zip(nf,nv))
                    rc=fields.get('RC_CNT','').split(',')
                    if len(rc)!=6 or any(int(x)<0 for x in rc):
                        raise ValueError('Native read-support vector must contain six measured counts')
                    result[9+i]=row[9+i]+':'+':'.join(nv)
            result[7]=';'.join(k if v is True else k+'='+str(v) for k,v in values.items())
            f.write('\t'.join(result)+'\n')
            counts[status]+=1
            if row[6]=='PASS':counts['PASS_'+status]+=1
    preservation=verify_annotations(source,out)
    result=dict(version=VERSION,counts=dict(counts),preservation=preservation,
                full_downstream_ready=False,
                note='Caller payload retained; native support depth is SAGE_DP. Native NEO must use that measured depth with RC_CNT, not mixed caller DP.')
    Path(str(out)+'.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def consumer_view(source,out,neo=False):
    """PASS-only consumer test, separate from the complete retained variant file.

    NEO computes support from RC_CNT and denominator from DP. Its isolated
    input therefore uses the corresponding SAGE_DP; caller DP remains DS_DP.
    No caller counts change in the master or PAVE/PURPLE inputs.
    """
    h=headers(source)
    rows=[(s,r) for s,r in read_rows(source) if r[6]=='PASS']
    for samples,row in rows:
        if info(row[7]).get('DS_EVIDENCE_STATUS')!='MEASURED':
            raise ValueError('PASS allele without unique native evidence: '+str(allele_key(row)))
    with open(out,'x') as f:
        for line in h:
            if line.startswith('#CHROM') and not any(x.startswith('##INFO=<ID=DS_CONTEXT_PROJECTED,') for x in h):
                for definition in h:
                    m=DEFINITION.match(definition)
                    if m and m.group(1)=='INFO' and m.group(2) in CONTEXT_TAGS:
                        tag=m.group(2)
                        f.write(definition.replace('ID='+tag+',','ID=DS_ORIGINAL_INFO_'+tag+',',1)+'\n')
                f.write('##INFO=<ID=DS_CONTEXT_PROJECTED,Number=0,Type=Flag,Description="Native context annotations used by this consumer; original annotations retained in DS_ORIGINAL_INFO tags and unchanged master">\n')
            if line.startswith('#CHROM') and neo:
                f.write('##FORMAT=<ID=DS_DP,Number=1,Type=Integer,Description="Original DeepSomatic depth before NEO-only context-depth projection">\n')
            f.write(line+'\n')
        for samples,row in rows:
            result=row[:]
            attrs=info(row[7])
            if 'DS_CONTEXT_PROJECTED' not in attrs:
                if 'SAGE_INFO_FIELDS' not in attrs:
                    raise ValueError('Native INFO presence metadata required; regenerate with adapter 0.1.2')
                present=set(str(attrs['SAGE_INFO_FIELDS']).split(','))
                for tag in CONTEXT_TAGS:
                    native_tag='SAGE_'+tag if 'SAGE_'+tag in attrs else tag
                    native_value=attrs.get(native_tag)
                    if tag in attrs:attrs['DS_ORIGINAL_INFO_'+tag]=attrs[tag]
                    if tag in present:
                        if native_value is None:raise ValueError('Native annotation payload missing: '+tag)
                        attrs[tag]=native_value
                    else:
                        attrs.pop(tag,None)
                attrs['DS_CONTEXT_PROJECTED']=True
                result[7]=';'.join(k if v is True else k+'='+str(v) for k,v in attrs.items())
            if neo:
                tags=row[8].split(':')
                if 'DS_DP' in tags:raise ValueError('Already projected')
                result[8]+=':DS_DP'
                for i in range(len(samples)):
                    values=row[9+i].split(':'); data=dict(zip(tags,values))
                    # RNA was measured by native append and has its own valid DP;
                    # the inherited SAGE_DP tag may be missing for this new sample.
                    measured=data.get('SAGE_DP','.')
                    depth=data['DP'] if measured in ('.','') else measured
                    values[tags.index('DP')]=depth
                    result[9+i]=':'.join(values+[data['DP']])
            f.write('\t'.join(result)+'\n')
    return len(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True);p.add_argument('--evidence',type=Path)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--candidates',type=Path,help='Native ordered CandidateWriter output, required to resolve multiple contexts')
    p.add_argument('--supplement',type=Path,help='Native append measurements from observed exact-CIGAR contexts missing from discovery')
    p.add_argument('--consumer-view',action='store_true');p.add_argument('--neo',action='store_true')
    a=p.parse_args()
    if a.consumer_view:print(consumer_view(a.source,a.out,a.neo))
    elif a.evidence:print(json.dumps(merge(a.source,a.evidence,a.out,a.candidates,a.supplement),indent=2))
    else:p.error('--evidence or --consumer-view required')


if __name__=='__main__':main()
