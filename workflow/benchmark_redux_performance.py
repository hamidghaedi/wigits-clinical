#!/usr/bin/env python3
"""Isolated performance experiment; never edits a native task or production defaults.

Replay an existing native REDUX task with64 CPUs, identical heap/inputs/scientific
options, varying only native partitioning/concatenation controls. Execute under
Slurm using benchmark_redux_performance.sbatch. All outputs stay in a new directory.
This measures speed; it does not automatically approve a variant for production.
"""
import argparse
import json
import os
from pathlib import Path
import re
import resource
import subprocess
import time

MODES = {'baseline': '', 'partition8': '-partition_ratio 8',
         'partition8_parallel': '-partition_ratio 8 -parallel_concat',
         'baseline_repeat': ''}  # bound cache/order effects with a second baseline


def experiment_script(original, mode):
    text, count = re.subn(r'-threads\s+\d+', '-threads 64', original)
    if count != 1:
        raise ValueError('Expected one native REDUX thread argument')
    text, count = re.subn(r'-Xmx\d+', '-Xmx128849018880', text)
    if count != 1 or '-partition_ratio' in text or '-parallel_concat' in text:
        raise ValueError('Expected an unmodified native task with one heap argument')
    return text.replace('-threads 64', '-threads 64 ' + MODES[mode])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native-work', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--image', type=Path, required=True)
    args = parser.parse_args()
    if int(os.environ.get('SLURM_CPUS_PER_TASK', '0')) != 64:
        parser.error('Run this experiment under a64-CPU Slurm allocation')
    args.out.mkdir(parents=True, exist_ok=False)
    source = (args.native_work/'.command.sh').read_text()
    results = dict(native_work=str(args.native_work.resolve()), image=str(args.image.resolve()),
                   threads=64, heap_bytes=128849018880, trials=[])
    for mode in MODES:
        work = args.out/mode
        work.mkdir()
        for path in args.native_work.iterdir():
            if path.is_symlink() and not path.name.startswith('.'):
                (work/path.name).symlink_to(path.resolve(), target_is_directory=path.is_dir())
        (work/'command.sh').write_text(experiment_script(source, mode))
        command = ['apptainer','exec','--cleanenv','--no-home','--bind','/global',
                   str(args.image),'bash','command.sh']
        started = time.monotonic()
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        with (work/'command.log').open('w') as log:
            code = subprocess.run(command, cwd=work, stdout=log, stderr=subprocess.STDOUT).returncode
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        elapsed = time.monotonic()-started
        cpu = after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime
        results['trials'].append(dict(mode=mode, flags=MODES[mode], exit_code=code,
            wall_seconds=elapsed, child_cpu_seconds=cpu, average_busy_cpus=cpu/elapsed,
            cpu_measurement='local child accounting; final Slurm accounting is authoritative'))
        (args.out/'measurements.json').write_text(json.dumps(results,indent=2)+'\n')
        print(json.dumps(results['trials'][-1]), flush=True)
        if code:
            raise SystemExit(code)


if __name__ == '__main__':
    main()
