#!/usr/bin/env python3
"""Recover all primary read pairs, without region/quality/duplicate filtering.

OQ is preferred where available; otherwise use BAM QUAL. Output is BAM-derived,
not an assertion of original sequencer FASTQ identity. Independently count
primary source records with flagstat, validate every FASTQ pair, and require
complete primary-record accounting. Singletons/other records are preserved but
prevent paired-only workflow release. All bulk work belongs in a Slurm job.
"""
import argparse
from datetime import datetime, timezone
import gzip
import itertools
import json
from pathlib import Path
import subprocess


def reads(path):
    with gzip.open(path, 'rt') as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            seq, plus, quality = [handle.readline().rstrip('\r\n') for _ in range(3)]
            if not header.startswith('@') or not plus.startswith('+') or len(seq) != len(quality):
                raise ValueError(f'Malformed FASTQ in {path}')
            yield header.rstrip().split()[0]


def validate_fastqs(outputs):
    pairs = 0
    for a, b in itertools.zip_longest(reads(outputs['R1']), reads(outputs['R2'])):
        if a is None or b is None or a != b:
            raise ValueError('FASTQ pair counts or names differ')
        pairs += 1
    return dict(read_pairs=pairs,
                singleton_reads=sum(1 for _ in reads(outputs['singletons'])),
                other_reads=sum(1 for _ in reads(outputs['other'])))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bam', type=Path, required=True)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--threads', type=int, required=True)
    args = parser.parse_args()
    if args.threads < 4:
        parser.error('At least four CPUs required')
    if not args.bam.is_file():
        parser.error('Source BAM missing')
    args.outdir.mkdir(parents=True, exist_ok=False)
    outputs = {key: args.outdir / f'{args.sample}.full_{suffix}.fastq.gz'
               for key, suffix in [('R1', 'R1'), ('R2', 'R2'),
                                   ('singletons', 'singletons'), ('other', 'other')]}
    commands = []
    with (args.outdir / 'commands.log').open('w') as log:
        def run(cmd, **kwargs):
            commands.append(cmd)
            print(json.dumps(cmd), file=log, flush=True)
            return subprocess.run(cmd, stderr=log, check=True, **kwargs)

        version = run(['samtools', '--version'], stdout=subprocess.PIPE).stdout.decode().splitlines()[0]
        run(['samtools', 'quickcheck', '-v', str(args.bam)])
        print(f'{args.sample}: counting all primary source records', flush=True)
        flagstat = json.loads(run(['samtools', 'flagstat', '-@', str(args.threads - 1),
                                  '-O', 'json', str(args.bam)], stdout=subprocess.PIPE).stdout)
        (args.outdir / 'source_flagstat.json').write_text(json.dumps(flagstat, indent=2) + '\n')
        primary_records = sum(flagstat[k]['primary'] for k in ['QC-passed reads', 'QC-failed reads'])
        workers = (args.threads - 2) // 2
        collate = ['samtools', 'collate', '-u', '-O', '-@', str(workers), '-T',
                   str(args.outdir / 'collate'), str(args.bam)]
        fastq = ['samtools', 'fastq', '-@', str(workers), '-n', '-O', '-F', '0x900',
                 '-1', str(outputs['R1']), '-2', str(outputs['R2']),
                 '-s', str(outputs['singletons']), '-0', str(outputs['other']), '-']
        print(f'{args.sample}: recovering every primary read pair', flush=True)
        for command in [collate, fastq]:
            commands.append(command)
            print(json.dumps(command), file=log, flush=True)
        with subprocess.Popen(collate, stdout=subprocess.PIPE, stderr=log) as process:
            try:
                subprocess.run(fastq, stdin=process.stdout, stderr=log, check=True)
            finally:
                process.stdout.close()
            if process.wait() != 0:
                raise RuntimeError('samtools collate failed')
        print(f'{args.sample}: validating all pairs and source accounting', flush=True)
        counts = validate_fastqs(outputs)
        if primary_records != 2 * counts['read_pairs'] + counts['singleton_reads'] + counts['other_reads']:
            raise ValueError('FASTQ outputs do not account for all primary source records')
        summary = dict(sample_id=args.sample, source_bam=str(args.bam.resolve()),
                       source_bam_bytes=args.bam.stat().st_size, regions='ALL',
                       downsampling=False, excluded_flags='0x900',
                       quality_source='OQ when present, otherwise BAM QUAL',
                       samtools_version=version, primary_records=primary_records,
                       **counts, outputs={k: str(p.resolve()) for k, p in outputs.items()},
                       output_bytes={k: p.stat().st_size for k, p in outputs.items()},
                       commands=commands, completed_at=datetime.now(timezone.utc).isoformat())
        (args.outdir / 'recovery.json').write_text(json.dumps(summary, indent=2) + '\n')
        if counts['singleton_reads'] or counts['other_reads'] or not counts['read_pairs']:
            raise ValueError('Paired-only release requires nonzero pairs and no omitted primary reads')
        (args.outdir / 'completion.json').write_text(json.dumps(summary, indent=2) + '\n')
        print(json.dumps(counts), flush=True)


if __name__ == '__main__':
    main()
