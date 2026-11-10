import gzip
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'workflow' / f'{name}.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ResourceTests(unittest.TestCase):
    def test_workload_increases_initial_resources(self):
        module = load('build_input_resources')
        small = [dict(role=r, read_pairs=50_000_000) for r in ['tumor', 'normal', 'rna']]
        large = [dict(role=r, read_pairs=3_000_000_000) for r in ['tumor', 'normal', 'rna']]
        low, high = module.plan(small), module.plan(large)
        for name in ['BWAMEM2_ALIGN', 'REDUX', 'ESVEE', 'SAGE_SOMATIC', 'ISOFOX']:
            self.assertEqual(high[name]['cpus'], 64)
            self.assertGreater(high[name]['memory_gb'], low[name]['memory_gb'])
        self.assertGreater(high['ESVEE']['time_hours'], low['ESVEE']['time_hours'])
        self.assertEqual(high['GATK4_MARKDUPLICATES']['cpus'], 8)

    def test_oversized_input_is_not_silently_clamped(self):
        with self.assertRaises(ValueError):
            load('build_input_resources').allocation(1000, 48, 2, 128, 6)

    def test_unpaired_or_subset_input_rejected(self):
        module = load('build_input_resources')
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'completion.json'
            path.write_text(json.dumps(dict(regions=['chr17', 'chr19'], downsampling=False, excluded_flags='0x900')))
            with self.assertRaises(ValueError):
                module.load_library('tumor', path)


@unittest.skipUnless(shutil.which('samtools'), 'samtools required for synthetic BAM integration')
class RecoveryTests(unittest.TestCase):
    def test_preserves_all_primary_pairs_and_original_quality(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            sam = tmp / 'input.sam'
            rows = ['@HD\tVN:1.6\tSO:unsorted', '@SQ\tSN:chr1\tLN:1000']
            for name, flags in [('mapped', [99, 147]), ('duplicate', [1123, 1171]), ('unmapped', [77, 141])]:
                for flag in flags:
                    unmapped = bool(flag & 4)
                    fields = [name, str(flag), '*' if unmapped else 'chr1', '0' if unmapped else '1',
                              '0' if unmapped else '60', '*' if unmapped else '10M',
                              '*' if unmapped else '=', '0' if unmapped else '1', '0',
                              'ACGTACGTAC', '!!!!!!!!!!', 'OQ:Z:IIIIIIIIII']
                    rows.append('\t'.join(fields))
            for flag in [355, 2147]:
                rows.append(f'mapped\t{flag}\tchr1\t1\t60\t10M\t=\t1\t0\tACGTACGTAC\t!!!!!!!!!!')
            sam.write_text('\n'.join(rows) + '\n')
            bam = tmp / 'input.bam'
            subprocess.run(['samtools', 'sort', '-o', str(bam), str(sam)], check=True)
            out = tmp / 'fastqs'
            subprocess.run(['python3', str(ROOT / 'workflow/prepare_full_fastqs.py'), '--bam', str(bam),
                            '--sample', 'synthetic', '--outdir', str(out), '--threads', '4'], check=True,
                           capture_output=True, text=True)
            data = json.loads((out / 'completion.json').read_text())
            self.assertEqual(data['read_pairs'], 3)
            self.assertEqual(data['primary_records'], 6)
            self.assertEqual(data['singleton_reads'], 0)
            with gzip.open(data['outputs']['R1'], 'rt') as handle:
                lines = handle.read().splitlines()
                self.assertEqual(lines[3::4], ['IIIIIIIIII'] * 3)


if __name__ == '__main__':
    unittest.main()
