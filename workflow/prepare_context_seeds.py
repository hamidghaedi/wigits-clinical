#!/usr/bin/env python3
"""Select observed native-builder contexts for discovery-missing alleles.

Only exact CIGAR events qualify. A substring match in a repeat can occur in a
reference read and is diagnostic only, even if the native context builder accepts
it. Count distinct original read names, select maximum support then minimum
native event count; an unresolved tie stays unresolved. Never construct bases,
alter qualities, or replace caller genotypes. All contexts remain in the audit.
"""
import argparse
from collections import defaultdict
import csv
import json
import re
from pathlib import Path
from deepsomatic_adapter import read_rows,allele_key,info
from merge_sage_evidence import headers


def exact_indel_context(key, context):
    """Require the serialized observed CIGAR to encode this exact allele."""
    chrom, pos, ref, alt = key.split(':')
    parts = context.split('-', 5)
    if len(parts) != 6:
        return False
    start, _, left, core, right, cigar = parts
    bases = left + core + right
    rp, ri = int(start), 0
    for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
        length = int(length)
        if rp - 1 == int(pos) and ri > 0:
            if op == 'I' and ref == alt[:1] and len(ref) == 1:
                if bases[ri - 1:ri + length] == alt:
                    return True
            if op == 'D' and len(ref) == length + 1 and alt == ref[:1]:
                if bases[ri - 1:ri] == alt:
                    return True
        if op in 'MIS=X':
            ri += length
        if op in 'MDN=X':
            rp += length
    return False


def prepare(source, observations, header_source, out, primary=None):
    targets={':'.join(allele_key(r)):r for s,r in read_rows(source) if r[6]=='PASS'}
    candidates=defaultdict(lambda:defaultdict(lambda:dict(reads=set(),events=[],methods=set())))
    with observations.open() as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['requested']!=r['observed'] or not r['context'][0].isdigit():continue
            if r['method']!='CIGAR':continue
            if not exact_indel_context(r['requested'], r['context']):continue
            item=candidates[r['requested']][r['context']]
            item['reads'].add(r['read_name']);item['events'].append(int(r['native_events']));item['methods'].add(r['method'])
    preferred={}
    if primary:
        with primary.open() as f:
            for r in csv.DictReader(f,delimiter='\t'):
                if r['role']=='PRIMARY':preferred[r['requested']]=r['context']
    audit={};selected=[]
    for key,contexts in candidates.items():
        choices=[dict(context=k,reads=len(v['reads']),native_events=min(v['events']),methods=sorted(v['methods'])) for k,v in contexts.items()]
        choices.sort(key=lambda x:(-x['reads'],x['native_events'],x['context']))
        best=choices[0];tie=sum((v['reads'],v['native_events'])==(best['reads'],best['native_events']) for v in choices)>1
        native=[v for v in choices if v['context']==preferred.get(key)]
        if key in preferred and len(native)!=1:
            raise ValueError('Native primary lacks exact observed-CIGAR evidence: '+key)
        if native:best=native[0];tie=False
        audit[key]=dict(contexts=choices,status='TIED_CONTEXTS' if tie else 'SELECTED_OBSERVED_CONTEXT',
                       selection='NATIVE_ALTCONTEXT' if native else 'OBSERVED_READ_COUNT_THEN_EVENTS')
        if not tie:selected.append((targets[key],best))
    with out.open('x') as f:
        f.write('\n'.join(x for x in headers(header_source) if not x.startswith('#CHROM'))+'\n')
        f.write('##INFO=<ID=DS_SEED_READS,Number=1,Type=Integer,Description="Distinct observed read names supporting the exported native-builder context">\n')
        f.write('##INFO=<ID=DS_SEED_METHOD,Number=1,Type=String,Description="Context from observed BAM read; not a native SAGE discovery or classification">\n')
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for row,best in selected:
            attrs=dict(TIER=info(row[7])['TIER'],RC_INFO=best['context'],RC_NM=best['native_events'],
                       DS_SEED_READS=best['reads'],DS_SEED_METHOD='OBSERVED_BAM_READ')
            f.write('\t'.join(row[:5]+['.','.',';'.join(k+'='+str(v) for k,v in attrs.items())])+'\n')
    Path(str(out)+'.json').write_text(json.dumps(audit,indent=2)+'\n')
    return len(selected)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['source','observations','header-source','out']:p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--primary',type=Path,help='Native AltContext primary selection from the same inspected reads')
    a=p.parse_args();print(prepare(a.source,a.observations,a.header_source,a.out,a.primary))
