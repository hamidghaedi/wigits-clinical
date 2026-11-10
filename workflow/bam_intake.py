#!/usr/bin/env python3
"""BAM preparation v0.3.0: YAML inputs, optional Picard repair, native execution.

No custom audits, checksums, reference/header validators, record probes or
receipt-based launch gates. Standard tool failures stop the job. Reference
and sample assignments are supplied by the operator.
"""
import argparse
import csv
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import shlex
import subprocess
import time

import yaml

VERSION = '0.3.0'
FIELDS = ['group_id', 'subject_id', 'sample_id', 'sample_type',
          'sequence_type', 'filetype', 'info', 'filepath']


def write_json(path, value):
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(value, indent=2) + '\n')
    os.replace(tmp, path)


def prepare(args):
    """Write a separate mate-prepared BAM and index; never modify the input."""
    args.outdir.mkdir(parents=True, exist_ok=False)
    work = args.outdir / 'sorting'
    work.mkdir()
    output = args.outdir / f'{args.sample}.oncoanalyser.bam'
    fixed = work / 'mate_fixed.bam'
    state = dict(version=VERSION, status='PREPARING', source=str(args.bam),
                 output=str(output), started_utc=datetime.now(timezone.utc).isoformat(), commands=[])
    write_json(args.outdir / 'preparation.json', state)
    with (args.outdir / 'commands.log').open('w') as log:
        name_sort = ['samtools', 'sort', '-N', '-@', str(max(1, args.threads-8)),
                     '-m', '1G', '-l', '0', '-T', str(work/'name'), '-O', 'BAM', str(args.bam)]
        fix = shlex.split(args.gatk_command) + [
            '--java-options', '-Xmx16g -XX:ParallelGCThreads=4', 'FixMateInformation',
            '--INPUT', '/dev/stdin', '--OUTPUT', str(fixed), '--ASSUME_SORTED', 'true',
            '--SORT_ORDER', 'unsorted', '--ADD_MATE_CIGAR', 'true',
            '--IGNORE_MISSING_MATES', 'false', '--COMPRESSION_LEVEL', '1', '--TMP_DIR', str(work)]
        coordinate_sort = ['samtools', 'sort', '-@', str(max(1, args.threads-1)),
                           '-m', '1G', '-T', str(work/'coordinate'), '-o', str(output), str(fixed)]
        index = ['samtools', 'index', '-@', str(max(1, args.threads-1)), str(output)]
        try:
            state['commands'] = [name_sort, fix, coordinate_sort, index]
            write_json(args.outdir/'preparation.json', state)
            print(f'{args.sample}: name sorting and Picard mate preparation', flush=True)
            first = subprocess.Popen(name_sort, stdout=subprocess.PIPE, stderr=log)
            try:
                # GATK stdout includes completion text; its BAM output must be a file.
                second = subprocess.Popen(fix, stdin=first.stdout, stdout=log, stderr=log)
                first.stdout.close()
                codes = [second.wait(), first.wait()]
                if any(codes):
                    raise RuntimeError(f'Name sort/Picard failed: {codes}; see commands.log')
            finally:
                if first.poll() is None:
                    first.terminate()
            print(f'{args.sample}: coordinate sorting and indexing', flush=True)
            subprocess.run(coordinate_sort, check=True, stderr=log)
            subprocess.run(index, check=True, stderr=log)
            state.update(status='PREPARED', completed_utc=datetime.now(timezone.utc).isoformat())
            write_json(args.outdir/'preparation.json', state)
            fixed.unlink()
            print(f'{args.sample}: PREPARED {output}', flush=True)
        except Exception as error:
            state.update(status='FAILED', error=str(error))
            write_json(args.outdir/'preparation.json', state)
            raise


def batch(args):
    """Prepare YAML rows marked prepare_mates, then write the native CSV."""
    config = yaml.safe_load(args.config.read_text())
    args.outdir.mkdir(parents=True, exist_ok=True)
    state = dict(version=VERSION, config=str(args.config), status='PREPARING', jobs={})
    rows = []
    for sample in config['samples']:
        row = {key:sample.get(key, '') for key in FIELDS}
        if sample.get('prepare_mates', False):
            name = sample['sample_id']
            destination = args.outdir / name
            size = Path(sample['filepath']).stat().st_size / (1024**3)
            hours = math.ceil(max(24, size/4*1.5)/4)*4
            command = ['sbatch', '--parsable', '--account='+args.account,
                       '--cpus-per-task=64', '--mem=192G', f'--time={hours}:00:00',
                       '--output='+str(args.outdir/f'{name}.%j.log'),
                       str(Path(__file__).with_name('bam_intake_worker.sbatch')),
                       str(Path(__file__).resolve()), sample['filepath'], name,
                       str(destination), args.gatk_command]
            job = subprocess.check_output(command, text=True).strip().split(';')[0]
            state['jobs'][name] = dict(job_id=job, command=command)
            row['filepath'] = str(destination/f'{name}.oncoanalyser.bam')
            write_json(args.outdir/'dispatch.json', state)
        rows.append(row)
    while state['jobs']:
        statuses = {name:subprocess.check_output(
            ['sacct', '-n', '-X', '-j', job['job_id'], '--format=State,ExitCode', '-P'],
            text=True).strip() for name,job in state['jobs'].items()}
        if any(any(failure in value for failure in
                   ['FAILED', 'CANCELLED', 'TIMEOUT', 'OUT_OF_MEMORY', 'NODE_FAIL', 'PREEMPTED'])
               for value in statuses.values()):
            raise RuntimeError('BAM preparation job failed: '+str(statuses))
        if all(value.startswith('COMPLETED|0:0') for value in statuses.values()):
            break
        time.sleep(30)
    with args.prepared_samplesheet.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    state.update(status='PREPARED', samplesheet=str(args.prepared_samplesheet))
    write_json(args.outdir/'dispatch.json', state)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest='mode', required=True)
    worker = modes.add_parser('prepare')
    worker.add_argument('--bam', type=Path, required=True)
    worker.add_argument('--sample', required=True)
    worker.add_argument('--outdir', type=Path, required=True)
    worker.add_argument('--threads', type=int, required=True)
    worker.add_argument('--gatk-command', required=True)
    batch_parser = modes.add_parser('batch')
    batch_parser.add_argument('--config', type=Path, required=True)
    batch_parser.add_argument('--outdir', type=Path, required=True)
    batch_parser.add_argument('--prepared-samplesheet', type=Path, required=True)
    batch_parser.add_argument('--gatk-command', required=True)
    batch_parser.add_argument('--account', required=True)
    args = parser.parse_args()
    prepare(args) if args.mode == 'prepare' else batch(args)


if __name__ == '__main__':
    main()
