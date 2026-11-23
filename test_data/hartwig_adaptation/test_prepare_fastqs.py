#!/usr/bin/env python3
"""Small offline regression for chromosome selection and FASTQ pair recovery."""
import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class ChromosomeFastqTest(unittest.TestCase):
    def test_complete_pairs_quality_and_unpaired_reads(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sam, bam = root / 'input.sam', root / 'input.bam'
            records = [
                ('cross', 65, 'chr17', 100, 60, 'chr1', 300, 'OQ:Z:JJJJJJJJJJ'),
                ('cross', 129, 'chr1', 300, 60, 'chr17', 100, ''),
                ('unmapped', 73, 'chr19', 100, 60, '*', 0, ''),
                ('unmapped', 133, '*', 0, 0, 'chr19', 100, ''),
                ('duplicate', 1601, 'chr19', 200, 0, '=', 250, ''),
                ('duplicate', 1665, 'chr19', 250, 0, '=', 200, ''),
                ('cross', 321, 'chr19', 400, 60, 'chr1', 300, ''),
                ('outside', 65, 'chr1', 500, 60, '=', 550, ''),
                ('outside', 129, 'chr1', 550, 60, '=', 500, ''),
                ('unrelated_unmapped', 77, '*', 0, 0, '*', 0, ''),
                ('unrelated_unmapped', 141, '*', 0, 0, '*', 0, ''),
                ('orphan', 65, 'chr17', 500, 60, 'chr1', 900, ''),
                ('other', 0, 'chr19', 600, 60, '*', 0, ''),
            ]
            with sam.open('w') as handle:
                handle.write('@HD\tVN:1.6\tSO:unsorted\n')
                for chrom in ['chr17', 'chr19', 'chr1']:
                    handle.write(f'@SQ\tSN:{chrom}\tLN:1000\n')
                for name, flag, chrom, pos, mapq, mate, matepos, tag in records:
                    cigar = '*' if flag & 4 else '10M'
                    line = f'{name}\t{flag}\t{chrom}\t{pos}\t{mapq}\t{cigar}\t{mate}\t{matepos}\t0\tACGTACGTAA\tIIIIIIIIII'
                    handle.write(line + ('\t' + tag if tag else '') + '\n')
            subprocess.run(['samtools', 'sort', '-o', str(bam), str(sam)], check=True)
            subprocess.run(['samtools', 'index', str(bam)], check=True)
            out = root / 'outputs'
            subprocess.run([sys.executable, str(Path(__file__).with_name('prepare_fastqs.py')),
                            '--bam', str(bam), '--sample', 'synthetic',
                            '--outdir', str(out), '--threads', '4'], check=True)
            result = json.loads((out / 'completion.json').read_text())
            self.assertEqual(result['selected_primary_records'], 8)
            self.assertEqual(result['read_pairs'], 3)
            self.assertEqual(result['singleton_reads'], 1)
            self.assertEqual(result['other_reads'], 1)
            with gzip.open(result['outputs']['R1'], 'rt') as handle:
                lines = handle.read().splitlines()
            reads = {lines[i][1:]: lines[i + 3] for i in range(0, len(lines), 4)}
            self.assertEqual(set(reads), {'cross', 'unmapped', 'duplicate'})
            self.assertEqual(reads['cross'], 'JJJJJJJJJJ')
            self.assertFalse((out / 'selected.bam').exists())
            rerun = subprocess.run([sys.executable, str(Path(__file__).with_name('prepare_fastqs.py')),
                                    '--bam', str(bam), '--sample', 'synthetic',
                                    '--outdir', str(out), '--threads', '4'],
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.assertNotEqual(rerun.returncode, 0)


if __name__ == '__main__':
    unittest.main()
