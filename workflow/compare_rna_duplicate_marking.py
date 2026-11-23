#!/usr/bin/env python3
"""Small real-Picard comparison; run under Slurm, never on real case reads.

The baseline and RNA setting use the same installed GATK container. Compare
record identities, ordinary duplicate flags and duplication metrics. Optical
counts are expected to differ and are deliberately not claimed as measured in
the RNA setting. No libraries or software are installed.
"""
import argparse
import json
from pathlib import Path
import subprocess
import time


def synthetic_sam():
    lines = ['@HD\tVN:1.6\tSO:coordinate', '@SQ\tSN:chr1\tLN:10000',
             '@RG\tID:rg1\tSM:SYNTHETIC\tLB:lib1\tPL:ILLUMINA']
    records = []
    for group, count in [(0, 12), (1, 1), (2, 2)]:
        start = 100 + group * 1000
        for n in range(count):
            name = f'INST:1:FLOW:1:1101:{1000 + n}:{1000 + n}'
            # Keep group names distinct while retaining real coordinate fields.
            name = name.replace('FLOW:', f'FLOW{group}:')
            for flag, pos, mate, length in [(99, start, start+200, 250),
                                          (147, start+200, start, -250)]:
                records.append((pos, f'{name}\t{flag}\tchr1\t{pos}\t60\t50M\t=\t{mate}\t{length}\t'
                                      + 'A'*50 + '\t' + 'I'*50 + '\tRG:Z:rg1\tMC:Z:50M'))
    return '\n'.join(lines + [line for _, line in sorted(records)]) + '\n'


def sam_records(path):
    return sorted(tuple(line.split('\t')[:11]) for line in path.read_text().splitlines()
                  if not line.startswith('@'))


def metrics(path):
    lines = path.read_text().splitlines()
    for i, line in enumerate(lines):
        if line.startswith('LIBRARY\t'):
            return dict(zip(line.split('\t'), lines[i+1].split('\t')))
    raise ValueError('Picard duplication metrics missing')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--container', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    (args.out/'input.sam').write_text(synthetic_sam())
    evidence = {}
    for mode, extra in [('baseline', []), ('rna', ['--READ_NAME_REGEX', 'null'])]:
        command = ['apptainer', 'exec', '--cleanenv', '-B', str(args.out.resolve()),
                   str(args.container.resolve()), 'gatk', '--java-options', '-Xmx2g',
                   'MarkDuplicates', '--INPUT', str(args.out/'input.sam'),
                   '--OUTPUT', str(args.out/f'{mode}.bam'), '--METRICS_FILE',
                   str(args.out/f'{mode}.metrics'), '--TMP_DIR', str(args.out), *extra]
        started = time.monotonic()
        with (args.out/f'{mode}.log').open('w') as log:
            subprocess.run(command, check=True, stdout=log, stderr=log)
        elapsed = time.monotonic()-started
        convert = ['apptainer', 'exec', '--cleanenv', '-B', str(args.out.resolve()),
                   str(args.container.resolve()), 'gatk', '--java-options', '-Xmx2g',
                   'SamFormatConverter', '--INPUT', str(args.out/f'{mode}.bam'),
                   '--OUTPUT', str(args.out/f'{mode}.sam')]
        with (args.out/f'{mode}.convert.log').open('w') as log:
            subprocess.run(convert, check=True, stdout=log, stderr=log)
        evidence[mode] = dict(command=command, seconds=elapsed,
                              metrics=metrics(args.out/f'{mode}.metrics'))
    assert sam_records(args.out/'baseline.sam') == sam_records(args.out/'rna.sam'), 'Read or flag mismatch'
    baseline, rna = evidence['baseline']['metrics'], evidence['rna']['metrics']
    for field in ['READ_PAIRS_EXAMINED', 'READ_PAIR_DUPLICATES', 'PERCENT_DUPLICATION']:
        assert baseline[field] == rna[field], field
    assert int(baseline['READ_PAIR_OPTICAL_DUPLICATES']) > 0
    assert int(rna['READ_PAIR_OPTICAL_DUPLICATES']) == 0
    evidence.update(record_comparison='identical first 11 SAM fields for all 30 records',
                    duplicate_flags_preserved=True, optical_count_assessed_in_rna=False)
    (args.out/'result.json').write_text(json.dumps(evidence, indent=2)+'\n')
    print('PASS: same reads and ordinary duplicate flags; optical counts intentionally omitted')


if __name__ == '__main__':
    main()
