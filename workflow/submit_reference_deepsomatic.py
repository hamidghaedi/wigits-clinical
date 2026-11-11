#!/usr/bin/env python3
"""Submit the established full caller after a reference control has completed.

The native control owns alignment/REDUX and all native downstream analyses.
This dispatch sizes three caller jobs from the actual published BAM sizes.
It does not claim that the experimental downstream adapter is complete.
Use a frozen copy beside deepsomatic_full.py and deepsomatic_full.sbatch.
"""
import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import subprocess

VERSION = '0.1.0'
BASELINE_BYTES = 433_692_554_748


def resources(bam_bytes):
    # Preserve the issued whole-genome limits; grow for larger inputs. The
    # example stage measured 19h48 for the baseline; 60h remains conservative.
    # Inference is still running, so 72h is a limit, not a measured duration.
    scale = max(1, bam_bytes / BASELINE_BYTES)
    memory = max(256, math.ceil((bam_bytes / 2**30) / 400) * 128)
    stages = {
        'examples': dict(cpus=64, mem_gb=memory, hours=math.ceil(60*scale/12)*12),
        'inference': dict(cpus=16, mem_gb=128, hours=math.ceil(72*scale/12)*12),
        'postprocess': dict(cpus=64, mem_gb=memory, hours=math.ceil(12*scale/12)*12),
    }
    if memory > 768 or any(v['hours'] > 336 for v in stages.values()):
        raise ValueError('Measured workload exceeds the established node/time envelope')
    return dict(version=VERSION, bam_bytes=bam_bytes, baseline_bam_bytes=BASELINE_BYTES,
                scale=scale, stages=stages,
                basis='Measured BAM bytes scaled against issued full-genome limits; '
                      'not a fitted runtime model; inference baseline not yet complete')


def commands(source, config_path, out, record):
    """Build independent CPU/GPU commands; caller wires afterok dependencies."""
    result = {}
    for stage, r in record['stages'].items():
        gpu = stage == 'inference'
        cmd = ['sbatch', '--parsable', '--kill-on-invalid-dep=yes',
               '--account=def-hpcg6049_' + ('gpu' if gpu else 'cpu'),
               f'--job-name=SEQC2-DS-{stage}', f"--cpus-per-task={r['cpus']}",
               f"--mem={r['mem_gb']}G", f"--time={r['hours']}:00:00",
               f'--output={out.parent}/{stage}.%j.log']
        if gpu:
            cmd += ['--gres=gpu:1']
        cmd += [str(source/'deepsomatic_full.sbatch'), str(source/'deepsomatic_full.py'),
                str(config_path), str(out), stage]
        result[stage] = cmd
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    root = a.out.resolve()
    root.mkdir(parents=True, exist_ok=True)
    submission = root/'submission.json'
    if submission.exists():
        raise FileExistsError('Submission exists; inspect recorded jobs instead of submitting twice')
    config = json.loads(a.config.read_text())
    if config['scope'] != 'whole_genome' or config['model_type'] != 'FFPE_WGS':
        raise ValueError('This comparison requires whole-genome FFPE_WGS')
    paths = [Path(config[key]) for key in ('tumor_bam', 'normal_bam')]
    for path in paths:
        if not path.is_file() or not Path(str(path)+'.bai').is_file():
            raise FileNotFoundError(f'Published REDUX BAM/index missing: {path}')
    r = resources(sum(path.stat().st_size for path in paths))
    (root/'resources.json').write_text(json.dumps(r, indent=2)+'\n')
    planned = commands(Path(__file__).resolve().parent, a.config.resolve(), root/'results', r)
    record = dict(version=VERSION, started_utc=datetime.now(timezone.utc).isoformat(), jobs={}, commands={})
    submission.write_text(json.dumps(record, indent=2)+'\n')
    previous = None
    for stage, cmd in planned.items():
        if previous:
            cmd.insert(1, '--dependency=afterok:'+previous)
        record['commands'][stage] = cmd
        submission.write_text(json.dumps(record, indent=2)+'\n')
        previous = subprocess.check_output(cmd, text=True).strip().split(';')[0]
        record['jobs'][stage] = previous
        submission.write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps(record['jobs']))


if __name__ == '__main__':
    main()
