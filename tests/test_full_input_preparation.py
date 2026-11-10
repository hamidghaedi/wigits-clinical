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
    def test_full_entry_uses_staged_lanes_and_native_hmf_reference(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            package = root/'package'
            (package/'workflow').mkdir(parents=True)
            # A sentinel proves this entry consumes the CSV, rather than
            # requiring BAM-recovery completion.json files for original FASTQs.
            (package/'samplesheet.csv').write_text('original lane inputs\n')
            (package/'workflow/build_input_resources.py').write_text(
                'import pathlib,sys\n'
                'p=pathlib.Path(sys.argv[sys.argv.index("--samplesheet")+1])\n'
                'assert p.read_text()=="original lane inputs\\n"\n'
                'for k in ["--out","--evidence"]:\n'
                ' pathlib.Path(sys.argv[sys.argv.index(k)+1]).write_text("sized")\n')
            marker = root/'native.json'
            (package/'workflow/run_oncoanalyser_frontenac.sh').write_text(
                '#!/usr/bin/env bash\n'
                'python3 - "$@" <<\'PY\'\n'
                'import json,os,pathlib,sys\n'
                f'pathlib.Path({str(marker)!r}).write_text(json.dumps(dict(args=sys.argv[1:],'
                'genome=os.environ["ONCO_GENOME"],source=os.environ["ONCO_SOURCE"])))\nPY\n')
            subprocess.run(['bash',str(ROOT/'workflow/run_full_native.sbatch'),str(package),'synthetic'],check=True)
            data = json.loads(marker.read_text())
            self.assertEqual(data['genome'],'GRCh38_hmf')
            self.assertEqual(data['source'],'nf-core/oncoanalyser')
            self.assertEqual(data['args'][:2],['synthetic',str(package/'samplesheet.csv')])
            self.assertEqual((root/'state/input_resources.config').read_text(),'sized')

    def test_existing_input_lanes_are_sized_without_inventing_read_counts(self):
        module = load('build_input_resources')
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            paths = []
            for i in range(4):
                p = tmp/f'lane{i}.fastq.gz'
                p.write_bytes(b'test')
                paths.append(str(p))
            sheet = tmp/'samplesheet.csv'
            sheet.write_text('sample_id,sample_type,sequence_type,filetype,filepath\n'
                + f'tumour,tumor,dna,fastq,{paths[0]};{paths[1]}\n'
                + f'tumour,tumor,dna,fastq,{paths[2]};{paths[3]}\n')
            result = module.load_samplesheet(sheet)
            self.assertEqual(len(result), 1)
            self.assertEqual(result[0]['input_bytes'], 16)
            self.assertIsNone(result[0]['read_pairs'])
            sheet.write_text(sheet.read_text() + f'tumour,tumor,dna,fastq,{paths[0]};{paths[1]}\n')
            with self.assertRaises(ValueError):
                module.load_samplesheet(sheet)

    def test_bam_byte_proxy_changes_allocation(self):
        module = load('build_input_resources')
        rows = [dict(role=r, read_pairs=None, input_bytes=200*1024**3) for r in ['tumor','normal']]
        result = module.plan(rows)
        self.assertEqual(result['REDUX']['cpus'], 64)
        self.assertEqual(result['ESVEE']['workload_units'], 40)

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
