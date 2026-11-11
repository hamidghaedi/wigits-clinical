#!/usr/bin/env python3
"""Replay native whole-genome example shards on two nodes, without editing inputs.

Component 0.1.0. Read the completed single-node native command; preserve every
scientific argument. Only output locations, task distribution and optional native
runtime TSVs differ. One Slurm allocation with two tasks on distinct nodes starts
both groups together. No GPU, inference or postprocessing job is submitted here.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import shlex
import socket
import subprocess
import time

VERSION = '0.1.0'


def now():
    return datetime.now(timezone.utc).isoformat()


def shard_ids(rank, groups, shards):
    if groups < 1 or shards < groups or not 0 <= rank < groups:
        raise ValueError('Invalid shard distribution')
    return list(range(rank, shards, groups))


def native_command(command, baseline, destination, shards=64):
    tokens = shlex.split(command)
    prefix = ['time', 'seq', '0', str(shards-1), '|', 'parallel', '-q', '--halt', '2', '--line-buffer']
    if tokens[:len(prefix)] != prefix:
        raise ValueError('Expected the pinned native GNU parallel command')
    argv = tokens[len(prefix):]
    if argv[0] != '/opt/deepvariant/bin/make_examples_somatic':
        raise ValueError('Expected native make_examples_somatic')
    if argv[-2:] != ['--task', '{}'] or any(x.startswith('--regions') for x in argv):
        raise ValueError('Expected whole-genome scope with native task placeholder')
    output_index = argv.index('--examples') + 1
    expected = str(baseline / 'intermediate' / f'make_examples_somatic.tfrecord@{shards}.gz')
    if argv[output_index] != expected:
        raise ValueError('Unexpected native shard pattern')
    argv[output_index] = str(destination / 'intermediate' / f'make_examples_somatic.tfrecord@{shards}.gz')
    argv += ['--runtime_by_region', str(destination / 'runtime' / f'regions@{shards}.tsv')]
    return argv


def prepare(config_path, baseline, out):
    baseline = baseline.resolve(); out = out.resolve()
    if out == baseline or baseline in out.parents or out in baseline.parents:
        raise ValueError('Use a separate output tree')
    config = json.loads(config_path.read_text())
    completed = json.loads((baseline/'examples.json').read_text())
    if completed.get('exit_code') != 0 or config.get('scope') != 'whole_genome':
        raise ValueError('Requires a completed whole-genome example baseline')
    source = json.loads((baseline/'commands.json').read_text())['examples']
    argv = native_command(source, baseline, out)
    for flag, key in [('--ref','reference'), ('--reads_tumor','tumor_bam'), ('--reads_normal','normal_bam')]:
        if argv[argv.index(flag)+1] != config[key]:
            raise ValueError('Baseline/config input mismatch: '+key)
    out.mkdir(parents=True, exist_ok=False)
    for name in ['intermediate','runtime','logs','workers']:
        (out/name).mkdir()
    elapsed = (datetime.fromisoformat(completed['finished_utc']) - datetime.fromisoformat(completed['started_utc'])).total_seconds()
    plan = dict(version=VERSION, config=config, baseline=str(baseline), out=str(out),
                native_argv=argv, groups=2, shards=64, cpus_per_node=32,
                baseline_seconds=elapsed, baseline_job=completed['job'],
                created_utc=now(), timing_caveats=['Different nodes/cache/load; runtime instrumentation added',
                                               'Same 64 shards tests node distribution, not more total CPUs'])
    (out/'plan.json').write_text(json.dumps(plan, indent=2)+'\n')
    return plan


def worker(out):
    plan = json.loads((out/'plan.json').read_text())
    rank = int(os.environ['SLURM_PROCID'])
    cpus = int(os.environ['SLURM_CPUS_PER_TASK'])
    if cpus != plan['cpus_per_node'] or int(os.environ['SLURM_NTASKS']) != plan['groups']:
        raise ValueError('Allocation does not match the experiment')
    ids = shard_ids(rank, plan['groups'], plan['shards'])
    rec = dict(rank=rank, host=socket.gethostname(), job=os.environ['SLURM_JOB_ID'],
               started_utc=now(), shard_ids=ids, cpus=cpus)
    path = out/'workers'/f'{rank}.json'
    path.write_text(json.dumps(rec, indent=2)+'\n')
    def run(task):
        argv = [str(task) if x == '{}' else x for x in plan['native_argv']]
        cmd = ['apptainer','exec','--cleanenv','--bind','/global',plan['config']['image']] + argv
        start = time.monotonic()
        item = dict(task=task, command=cmd, started_utc=now())
        with (out/'logs'/f'shard-{task:05d}.log').open('x') as log:
            item['exit_code'] = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT).returncode
        item.update(seconds=time.monotonic()-start, finished_utc=now())
        (out/'workers'/f'shard-{task:05d}.json').write_text(json.dumps(item, indent=2)+'\n')
        return item
    with ThreadPoolExecutor(max_workers=cpus) as pool:
        results = list(pool.map(run, ids))
    code = 0 if all(x['exit_code']==0 for x in results) else 1
    rec.update(exit_code=code, finished_utc=now(), shards=results)
    path.write_text(json.dumps(rec, indent=2)+'\n')
    return code


def digest(path):
    h = hashlib.sha256()
    opener = gzip.open if path.suffix == '.gz' else open
    size = 0
    with opener(path, 'rb') as handle:
        while chunk := handle.read(8*1024*1024):
            h.update(chunk); size += len(chunk)
    return dict(sha256=h.hexdigest(), uncompressed_bytes=size)


def compare(out):
    plan = json.loads((out/'plan.json').read_text()); baseline = Path(plan['baseline'])
    workers = [json.loads((out/'workers'/f'{i}.json').read_text()) for i in range(plan['groups'])]
    if any(x.get('exit_code') != 0 for x in workers):
        raise ValueError('Example groups not completed successfully')
    expected = []
    for i in range(plan['shards']):
        for stem in ['make_examples_somatic.tfrecord','make_examples_somatic_call_variant_outputs.tfrecord']:
            expected.append(f'{stem}-{i:05d}-of-{plan["shards"]:05d}.gz')
        expected.append(f'make_examples_somatic.tfrecord-{i:05d}-of-{plan["shards"]:05d}.gz.example_info.json')
    def compare_file(name):
        left = digest(baseline/'intermediate'/name)
        right = digest(out/'intermediate'/name)
        return dict(file=name, baseline=left, distributed=right, identical=left==right)
    with ThreadPoolExecutor(max_workers=int(os.environ.get('SLURM_CPUS_PER_TASK','4'))) as pool:
        files = list(pool.map(compare_file, expected))
    seconds = (max(datetime.fromisoformat(x['finished_utc']) for x in workers) -
               min(datetime.fromisoformat(x['started_utc']) for x in workers)).total_seconds()
    report = dict(version=VERSION, baseline_seconds=plan['baseline_seconds'], distributed_seconds=seconds,
                  speed_ratio=plan['baseline_seconds']/seconds, nodes=[x['host'] for x in workers],
                  identical_files=sum(x['identical'] for x in files), total_files=len(files), files=files,
                  caveats=plan['timing_caveats'],
                  interpretation='Decompressed byte equality is strong replay evidence; differences need semantic review, not automatic scientific disagreement.')
    (out/'comparison.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='files'}, indent=2))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=['prepare','worker','compare'])
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--config', type=Path); p.add_argument('--baseline', type=Path)
    a=p.parse_args()
    if a.mode=='prepare':
        if not a.config or not a.baseline: p.error('prepare needs --config and --baseline')
        prepare(a.config,a.baseline,a.out)
    elif a.mode=='worker': raise SystemExit(worker(a.out))
    else: compare(a.out)


if __name__=='__main__': main()
