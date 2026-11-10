#!/usr/bin/env python3
"""Bounded Google DeepSomatic experiment on native-reference tumour/normal BAMs.

Not a production caller switch. Inputs and output root are explicit; original
BAMs and native runs are read-only. Retain the raw VCF, including caller FILTERs.
No invented normal genotypes, SAGE read context, phasing or side-output files.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import json
import os
from pathlib import Path
import subprocess

VERSION = '0.1.0'


def command(config, out, cpus, cuda):
    if config['model_type'] not in ('WGS','FFPE_WGS'):
        raise ValueError('This paired-WGS experiment requires explicit WGS or FFPE_WGS')
    if not config['regions']:
        raise ValueError('This first test requires an explicit bounded region scope')
    if not cuda:
        raise ValueError('Submit the experimental job with one allocated GPU')
    return ['apptainer','exec','--cleanenv','--nv','--bind','/global',
        '--env',f'CUDA_VISIBLE_DEVICES={cuda}',
        str(config['image']),'run_deepsomatic',
        f"--model_type={config['model_type']}",f"--ref={config['reference']}",
        f"--reads_tumor={config['tumor_bam']}",f"--reads_normal={config['normal_bam']}",
        f"--sample_name_tumor={config['tumor_id']}",f"--sample_name_normal={config['normal_id']}",
        f'--output_vcf={out}/deepsomatic.vcf.gz',f'--num_shards={cpus}',
        f'--intermediate_results_dir={out}/intermediate',f'--logging_dir={out}/logs',
        f"--regions={' '.join(config['regions'])}",
        f'--postprocess_variants_extra_args=cpus={cpus},num_partitions={3*cpus}']


def summarize_vcf(path):
    counts, samples, formats, infos, versions = Counter(), [], set(), set(), []
    with gzip.open(path,'rt') as handle:
        for line in handle:
            if line.startswith('##') and ('version' in line.lower() or 'source=' in line):
                versions.append(line.rstrip())
            if line.startswith('#CHROM'):
                samples=line.rstrip().split('\t')[9:]
            elif not line.startswith('#'):
                fields=line.rstrip().split('\t')
                counts[fields[6]]+=1
                formats.update(fields[8].split(':'))
                infos.update(v.split('=',1)[0] for v in fields[7].split(';') if v != '.')
    return dict(records=sum(counts.values()),filters=dict(counts),samples=samples,
                format_fields=sorted(formats),info_fields=sorted(infos),version_headers=versions)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    config=json.loads(args.config.read_text())
    cpus=int(os.environ['SLURM_CPUS_PER_TASK'])
    cmd=command(config,args.out.resolve(),cpus,os.environ.get('CUDA_VISIBLE_DEVICES',''))
    args.out.mkdir(parents=True,exist_ok=False)
    record=dict(component_version=VERSION,caller='Google DeepSomatic',
        expected_caller_version='1.10.0',config=config,command=cmd,
        started_utc=datetime.now(timezone.utc).isoformat(),slurm_job=os.environ['SLURM_JOB_ID'])
    (args.out/'execution.json').write_text(json.dumps(record,indent=2)+'\n')
    with (args.out/'runner.log').open('w') as log:
        code=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT).returncode
    record.update(exit_code=code,finished_utc=datetime.now(timezone.utc).isoformat())
    if code==0:
        record['vcf_summary']=summarize_vcf(args.out/'deepsomatic.vcf.gz')
    (args.out/'execution.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(dict(exit_code=code,vcf_summary=record.get('vcf_summary')),indent=2))
    raise SystemExit(code)


if __name__=='__main__':
    main()
