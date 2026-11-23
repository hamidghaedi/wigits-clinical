#!/usr/bin/env python3
"""Export all chr17/chr19 primary reads and their mates; never downsample.

Uses indexed samtools view --fetch-pairs, name collation, then samtools fastq.
Duplicate/QC-failed/low-MAPQ reads are retained; secondary/supplementary records
are excluded to emit each sequencing read once. OQ is preferred where present,
otherwise BAM QUAL is used. These are BAM-derived, not original raw FASTQs.
Singletons and reads without a mate designation are retained in separate files.
"""
import argparse
import gzip
import itertools
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def reads(path):
    with gzip.open(path, 'rt') as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            sequence, plus, quality = [handle.readline().rstrip('\n') for _ in range(3)]
            if not header.startswith('@') or not plus.startswith('+') or len(sequence) != len(quality):
                raise ValueError(f'Malformed FASTQ record in {path}')
            yield header.rstrip('\n').split()[0]


def validate_fastqs(r1, r2, singles, other):
    pairs = 0
    for left, right in itertools.zip_longest(reads(r1), reads(r2)):
        if left is None or right is None or left != right:
            raise ValueError('R1/R2 read names or counts differ')
        pairs += 1
    return {'read_pairs': pairs,
            'singleton_reads': sum(1 for _ in reads(singles)),
            'other_reads': sum(1 for _ in reads(other))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bam', required=True, type=Path)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--outdir', required=True, type=Path)
    parser.add_argument('--threads', type=int, default=8)
    args = parser.parse_args()
    if args.threads < 4:
        parser.error('Use at least four allocated CPUs')
    if not args.bam.is_file() or not Path(str(args.bam) + '.bai').is_file():
        parser.error('The source BAM and .bam.bai index must exist')
    # A fresh directory prevents accidental replacement of a previous test set.
    args.outdir.mkdir(parents=True, exist_ok=False)
    prefix = args.outdir / (args.sample + '.chr17_chr19')
    selected = args.outdir / 'selected.bam'
    outputs = {k: Path(str(prefix) + suffix) for k, suffix in {
        'R1': '_R1.fastq.gz', 'R2': '_R2.fastq.gz',
        'singletons': '_singletons.fastq.gz', 'other': '_other.fastq.gz'}.items()}
    commands = []
    with (args.outdir / 'commands.log').open('w') as log:
        def run(command, **kwargs):
            commands.append(command)
            print(json.dumps(command), file=log, flush=True)
            return subprocess.run(command, check=True, stderr=log, **kwargs)

        version = run(['samtools', '--version'], capture_output=False,
                      stdout=subprocess.PIPE).stdout.decode().splitlines()[0]
        header = run(['samtools', 'view', '-H', str(args.bam)],
                     stdout=subprocess.PIPE).stdout.decode()
        refs = {field[3:] for line in header.splitlines() if line.startswith('@SQ\t')
                for field in line.split('\t') if field.startswith('SN:')}
        if not {'chr17', 'chr19'} <= refs:
            raise ValueError('Source BAM must contain chr17 and chr19')
        print(f'{args.sample}: extracting complete chromosome scope and mates', flush=True)
        run(['samtools', 'view', '-b', '-M', '--fetch-pairs', '-F', '0x900',
             '-@', str(args.threads - 1), '-o', str(selected), str(args.bam), 'chr17', 'chr19'])
        print(f'{args.sample}: converting to paired FASTQs', flush=True)
        workers = max(1, (args.threads - 2) // 2)
        collate = ['samtools', 'collate', '-u', '-O', '-@', str(workers),
                   '-T', str(args.outdir / 'collate'), str(selected)]
        fastq = ['samtools', 'fastq', '-@', str(workers), '-n', '-O',
                 '-1', str(outputs['R1']), '-2', str(outputs['R2']),
                 '-s', str(outputs['singletons']), '-0', str(outputs['other']), '-']
        commands.extend([collate, fastq])
        print(json.dumps(collate), file=log, flush=True)
        print(json.dumps(fastq), file=log, flush=True)
        with subprocess.Popen(collate, stdout=subprocess.PIPE, stderr=log) as process:
            try:
                subprocess.run(fastq, stdin=process.stdout, stderr=log, check=True)
            finally:
                process.stdout.close()
            if process.wait() != 0:
                raise RuntimeError('samtools collate failed')
        print(f'{args.sample}: checking paired read names and counts', flush=True)
        counts = validate_fastqs(outputs['R1'], outputs['R2'],
                                outputs['singletons'], outputs['other'])
        selected_records = int(run(['samtools', 'view', '-c', '-@',
                                    str(args.threads - 1), str(selected)],
                                   stdout=subprocess.PIPE).stdout)
        if selected_records != 2 * counts['read_pairs'] + counts['singleton_reads'] + counts['other_reads']:
            raise ValueError('FASTQs do not account for every selected primary BAM record')
        if counts['read_pairs'] == 0:
            raise ValueError('No paired reads were extracted')
        summary = {'sample_id': args.sample, 'source_bam': str(args.bam.resolve()),
                   'source_bam_bytes': args.bam.stat().st_size,
                   'regions': ['chr17', 'chr19'], 'include_mates_outside_regions': True,
                   'downsampling': False, 'excluded_flags': '0x900',
                   'quality_source': 'OQ when present, otherwise BAM QUAL',
                   'samtools_version': version, 'selected_primary_records': selected_records,
                   **counts, 'outputs': {k: str(v.resolve()) for k, v in outputs.items()},
                   'output_bytes': {k: v.stat().st_size for k, v in outputs.items()},
                   'completed_at': datetime.now(timezone.utc).isoformat(), 'commands': commands}
        (args.outdir / 'completion.json').write_text(json.dumps(summary, indent=2) + '\n')
        selected.unlink()  # Generated intermediate only; source BAM remains unchanged.
        print(json.dumps(counts), flush=True)


if __name__ == '__main__':
    main()
