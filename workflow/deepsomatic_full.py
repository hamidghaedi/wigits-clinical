#!/usr/bin/env python3
"""Whole-genome Google DeepSomatic using native generated commands in three jobs.

Google run_deepsomatic --dry_run supplies all model arguments. CPU preparation,
GPU inference and CPU postprocessing share its intermediates; no region filter,
caller changes, BAM rewriting or fabricated normal evidence is introduced.
"""
import argparse
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import re
import subprocess

VERSION = '0.2.1'
STAGES = ('examples', 'inference', 'postprocess')


def resources(config):
    bases = sum(int(line.split('\t')[1]) for line in Path(str(config['reference'])+'.fai').read_text().splitlines())
    bam_bytes = sum(Path(config[k]).stat().st_size for k in ('tumor_bam', 'normal_bam'))
    scale = bases / 20_000_000
    # Same-case 20Mb pilot timings. CPU scaling is an estimate, with 2x
    # uncertainty margin; this is not a whole-genome runtime prediction.
    def hours(seconds, speedup, floor):
        return max(floor, math.ceil(seconds * scale * 2 / speedup / 3600 / 12)*12)
    memory = max(256, math.ceil((bam_bytes / 2**30) / 400)*128)
    return dict(reference_bases=bases, bam_bytes=bam_bytes,
        basis='Same-case 20Mb pilot; genome-length proxy; 2x margin; CPU 16-to-64 scaling estimated at 4x; not runtime-validated',
        examples=dict(cpus=64, mem_gb=memory, hours=hours(2293.329,4,24)),
        inference=dict(cpus=16, mem_gb=128, hours=hours(733.125,1,24)),
        postprocess=dict(cpus=64, mem_gb=memory, hours=hours(112.872,4,12)))


def container(config, gpu=False):
    cmd=['apptainer','exec','--cleanenv','--bind','/global']
    if gpu:
        cuda=os.environ.get('CUDA_VISIBLE_DEVICES')
        if not cuda:
            raise ValueError('Inference requires the Slurm GPU mask')
        cmd += ['--nv','--env',f'CUDA_VISIBLE_DEVICES={cuda}']
    return cmd + [config['image']]


def dry_command(config, out):
    if config.get('regions') or config.get('scope') != 'whole_genome':
        raise ValueError('Whole-genome config must explicitly select whole_genome without regions')
    if config['model_type'] not in ('FFPE_WGS','WGS'):
        raise ValueError('Paired WGS model required')
    return container(config)+['run_deepsomatic','--dry_run=true',
        f"--model_type={config['model_type']}",f"--ref={config['reference']}",
        f"--reads_tumor={config['tumor_bam']}",f"--reads_normal={config['normal_bam']}",
        f"--sample_name_tumor={config['tumor_id']}",f"--sample_name_normal={config['normal_id']}",
        f'--output_vcf={out}/deepsomatic.vcf.gz','--num_shards=64',
        f'--intermediate_results_dir={out}/intermediate',
        '--call_variants_extra_args=batch_size=256',
        '--postprocess_variants_extra_args=cpus=64,num_partitions=192']


def parse_commands(text):
    commands=re.findall(r'\*\*\*\*\* Running the command:\*\*\*\*\*\s*\n([^\n]+)',text)
    if len(commands)!=3:
        raise ValueError(f'Expected three native commands, received {len(commands)}')
    for cmd,tool in zip(commands,('make_examples_somatic','call_variants','postprocess_variants')):
        if tool not in cmd or re.search(r'--regions(?:=|\s)',cmd):
            raise ValueError('Unexpected native command or restricted region scope')
    return dict(zip(STAGES,commands))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--stage',choices=STAGES,required=True)
    a=p.parse_args(); config=json.loads(a.config.read_text()); out=a.out.resolve()
    expected=16 if a.stage=='inference' else 64
    if int(os.environ['SLURM_CPUS_PER_TASK'])!=expected:
        raise ValueError(f'{a.stage} requires {expected} CPUs')
    if a.stage=='examples':
        out.mkdir(parents=True,exist_ok=False)
        with (out/'native_dry_run.log').open('w') as log:
            subprocess.run(dry_command(config,out),stdout=log,stderr=subprocess.STDOUT,check=True)
        commands=parse_commands((out/'native_dry_run.log').read_text())
        (out/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
    commands=json.loads((out/'commands.json').read_text())
    cmd=container(config,a.stage=='inference')+['bash','-o','pipefail','-c',commands[a.stage]]
    record=dict(version=VERSION,stage=a.stage,job=os.environ['SLURM_JOB_ID'],
        started_utc=datetime.now(timezone.utc).isoformat(),command=cmd,scope='whole_genome')
    record_path=out/f'{a.stage}.json'
    record_path.write_text(json.dumps(record,indent=2)+'\n')
    with (out/f'{a.stage}.log').open('w') as log:
        code=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT).returncode
    record.update(exit_code=code,finished_utc=datetime.now(timezone.utc).isoformat())
    record_path.write_text(json.dumps(record,indent=2)+'\n')
    raise SystemExit(code)


if __name__=='__main__':
    main()
