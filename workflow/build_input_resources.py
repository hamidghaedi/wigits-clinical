#!/usr/bin/env python3
"""Generate initial site resources from measured existing inputs or recovery records.

Policy input-v2: units = ceil(primary reads / 250 million), or compressed
input GiB / 10 when a verified read count is unavailable. Heavy threaded
processes use 64 CPUs above 250 million input reads, otherwise 16.
Memory and walltime increase with stage-relevant units, with explicit
floors. This is conservative first-full-case extrapolation, not a fitted or
universally runtime-validated performance model. Native scientific code and
arguments remain unchanged. No retry-based sizing is substituted for this plan.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path


def allocation(units, floor_hours, hours_per_unit, base_gb, gb_per_unit, threaded=True):
    cpus = (64 if units >= 2 else 16) if threaded else 8
    hours = int(math.ceil(max(floor_hours, hours_per_unit * units) / 4) * 4)
    memory = int(math.ceil((base_gb + gb_per_unit * units + cpus) / 16) * 16)
    if hours > 336 or memory > 768:
        raise ValueError('Workload exceeds input-v2 eligible resource envelope; issue a reviewed plan')
    return dict(cpus=cpus, memory_gb=memory, time_hours=hours, workload_units=units)


def plan(libraries):
    def units(rows):
        return max(1, math.ceil(sum(
            2 * r['read_pairs'] / 250_000_000 if r.get('read_pairs') is not None
            else r['input_bytes'] / (10 * 1024**3) for r in rows)))
    dna = [r for r in libraries if r['role'] in ('tumor', 'normal')]
    tumor = [r for r in dna if r['role'] == 'tumor']
    normal = [r for r in dna if r['role'] == 'normal']
    rna = [r for r in libraries if r['role'] == 'rna']
    if len(tumor) != 1 or len(normal) != 1 or len(rna) > 1:
        raise ValueError('input-v2 expects one tumour/normal DNA pair and optional tumour RNA')
    d, t, n = units(dna), units(tumor), units(normal)
    r = units(rna)
    largest = max(units([x]) for x in libraries)
    specs = {
        'BWAMEM2_ALIGN': (max(t, n), 24, 2, 96, 2),
        'REDUX': (max(t, n), 16, 1.5, 64, 3),
        'SAGE_SOMATIC': (d, 24, 1.5, 96, 4),
        'SAGE_GERMLINE': (n, 16, 1.5, 96, 4),
        'ESVEE': (d, 48, 2, 128, 6),
        'AMBER': (d, 8, 0.5, 64, 2),
        'COBALT': (d, 8, 0.5, 64, 2),
        'VIRUSBREAKEND': (t, 24, 2, 96, 4),
        'BAMTOOLS': (largest, 8, 0.5, 32, 1.5),
        'CIDER': (largest, 12, 0.75, 32, 2),
        'LILAC': (d, 8, 0.5, 32, 2),
        'PURPLE': (d, 12, 0.5, 32, 2),
        'SAGE_APPEND.*': (max(t, r), 12, 0.75, 32, 2),
        'LINX_VISUALISER': (d, 8, 0.5, 32, 1),
    }
    if rna:
        specs.update({'STAR_ALIGN': (r, 12, 2, 64, 3),
                      'ISOFOX': (r, 12, 2, 64, 4),
                      'SAMTOOLS_SORT': (r, 8, 1, 32, 4)})
    out = {name: allocation(*args) for name, args in specs.items()}
    if rna:
        # Picard MarkDuplicates is not a 64-thread algorithm.
        out['GATK4_MARKDUPLICATES'] = allocation(r, 12, 1.5, 64, 4, threaded=False)
    return out


def load_library(role, path):
    data = json.loads(path.read_text())
    if data['regions'] != 'ALL' or data['downsampling'] or data['excluded_flags'] != '0x900':
        raise ValueError('Only complete, non-downsampled primary-read recovery is accepted')
    if data['singleton_reads'] or data['other_reads'] or data['read_pairs'] <= 0:
        raise ValueError('Incomplete paired input')
    if data['primary_records'] != 2 * data['read_pairs']:
        raise ValueError('Source record accounting mismatch')
    for mate in ['R1', 'R2']:
        if Path(data['outputs'][mate]).stat().st_size != data['output_bytes'][mate]:
            raise ValueError('FASTQ changed after completion')
    return dict(role=role, sample_id=data['sample_id'], read_pairs=data['read_pairs'],
                fastq_bytes=sum(data['output_bytes'][k] for k in ['R1', 'R2']),
                completion=str(path.resolve()),
                completion_sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def load_samplesheet(path):
    """Size existing BAM/FASTQ inputs; bytes are a proxy, never a read count.

    Ten GiB of compressed input per workload unit is a conservative alternative
    where no verified count is supplied. Existing names/paths remain unchanged;
    FASTQ lanes aggregate per sample, and duplicate source paths are refused.
    """
    groups, seen = {}, set()
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        if row['filetype'] in ('bai', 'crai'):
            continue
        if row['filetype'] not in ('bam', 'cram', 'fastq'):
            raise ValueError('This entry accepts alignment or FASTQ inputs only')
        role = 'rna' if row['sequence_type'] == 'rna' else row['sample_type']
        key = (role, row['sample_id'])
        data = groups.setdefault(key, dict(role=role, sample_id=row['sample_id'],
            read_pairs=None, input_bytes=0, files=[], workload_basis='compressed_input_bytes_proxy'))
        paths = row['filepath'].split(';')
        if row['filetype'] == 'fastq' and len(paths) != 2:
            raise ValueError('Every FASTQ row must specify both mates')
        for name in paths:
            p = Path(name)
            if str(p.resolve()) in seen:
                raise ValueError('Duplicate input path')
            seen.add(str(p.resolve()))
            if not p.is_file() or p.stat().st_size <= 0:
                raise ValueError(f'Input missing or empty: {p}')
            data['files'].append(dict(path=str(p.resolve()), bytes=p.stat().st_size))
            data['input_bytes'] += p.stat().st_size
    return list(groups.values())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    inputs = p.add_mutually_exclusive_group(required=True)
    inputs.add_argument('--library', nargs=2, action='append', metavar=('ROLE', 'COMPLETION_JSON'))
    inputs.add_argument('--samplesheet', type=Path)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--evidence', type=Path, required=True)
    a = p.parse_args()
    libraries = (load_samplesheet(a.samplesheet) if a.samplesheet else
                 [load_library(role, Path(path)) for role, path in a.library])
    resources = plan(libraries)
    lines = ['// input-v2: generated from measured input workload; see evidence JSON for basis.',
             'process {', '    resourceLimits = [cpus: 64, memory: 768.GB, time: 14.d]']
    for name, v in resources.items():
        lines.extend([f"    withName: '{name}' {{", f"        cpus = {v['cpus']}",
                      f"        memory = {v['memory_gb']}.GB", f"        time = {v['time_hours']}.h", '    }'])
    lines.append('}')
    a.out.write_text('\n'.join(lines) + '\n')
    a.evidence.write_text(json.dumps(dict(policy='input-v2', status='conservative extrapolation',
        created_at=datetime.now(timezone.utc).isoformat(), libraries=libraries, resources=resources), indent=2) + '\n')


if __name__ == '__main__':
    main()
