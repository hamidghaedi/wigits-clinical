#!/usr/bin/env python3
"""Compare measured WgsMetrics and depth histograms without changing either.

No averaging/merging, acceptance tolerance, scientific threshold or relabelling.
All numeric values are compared as decimals; differences remain explicit.
"""
import argparse
import csv
from decimal import Decimal,InvalidOperation
import json
from pathlib import Path


def sections(path):
    metric=None;hist={};lines=Path(path).read_text().splitlines()
    for i,line in enumerate(lines):
        if line.startswith('GENOME_TERRITORY\t'):
            metric=dict(zip(line.split('\t'),lines[i+1].split('\t')))
        if line.startswith('coverage\t'):
            for row in lines[i+1:]:
                if not row or row.startswith('#'):break
                k,v=row.split('\t')[:2];hist[k]=v
    if metric is None or not hist:raise ValueError('WgsMetrics row/depth histogram missing: '+str(path))
    return metric,hist


def equal(a,b):
    try:return Decimal(a)==Decimal(b)
    except (InvalidOperation,TypeError,ValueError):return a==b


def compare(baseline,candidate):
    old,oh=sections(baseline);new,nh=sections(candidate)
    differences=[dict(field=k,picard=old.get(k),parabricks=new.get(k))
                 for k in sorted(set(old)|set(new)) if not equal(old.get(k),new.get(k))]
    bins=[dict(depth=k,picard=oh.get(k),parabricks=nh.get(k))
          for k in sorted(set(oh)|set(nh),key=int) if not equal(oh.get(k),nh.get(k))]
    return dict(picard=str(baseline),parabricks=str(candidate),
                metrics_equal=not differences,histogram_equal=not bins,
                metric_differences=differences,histogram_differences=bins,
                adoption='not decided; timing, memory and differences require review')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--baseline',required=True);p.add_argument('--candidate',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();Path(a.out).write_text(json.dumps(compare(a.baseline,a.candidate),indent=2)+'\n')

if __name__=='__main__':main()
