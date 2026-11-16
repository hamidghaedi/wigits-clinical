#!/usr/bin/env python3
"""Release reproducible bulk files after a completed native run.

Keep all published results (including final alignments), source manifests,
logs and commands. Save non-alignment results to the explicit release directory.
Only named restored inputs/archives and large intermediate BAM/CRAM/FASTQ files
inside this run are eligible. Nearline custody stays in the input plan/receipts.
No data hashing, archive recall or deletion outside the named run is performed.
"""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import stat
import subprocess

VERSION = '0.1.0'


def completed(job):
    reply = subprocess.run(['sacct', '-n', '-X', '-j', str(job), '--format=State,ExitCode', '-P'],
                           check=True, text=True, capture_output=True)
    rows = [line.split('|') for line in reply.stdout.splitlines() if line.strip()]
    if len(rows) != 1 or rows[0][:2] != ['COMPLETED', '0:0']:
        raise ValueError('Native head has not completed successfully: ' + str(job))


def identity(path):
    info = path.lstat()
    if not stat.S_ISREG(info.st_mode):
        raise ValueError('Expected a regular bulk file: ' + str(path))
    return dict(path=str(path), bytes=info.st_size, device=info.st_dev, inode=info.st_ino,
                mtime_ns=info.st_mtime_ns, allocated_bytes=info.st_blocks * 512, links=info.st_nlink)


def candidates(root, inputs):
    """Return explicit paths, with no symlink traversal into other runs."""
    found = {}
    for library in inputs['libraries']:
        receipt = json.loads(Path(library['receipt']).read_text())
        if receipt['status'] != 'complete':
            raise ValueError('Restoration is incomplete')
        paths = [(Path(library['archive']), None)]
        paths += [(Path(row['path']), row['bytes']) for row in receipt['files']]
        for path, expected in paths:
            if not path.exists():
                continue
            if not path.is_relative_to(root / 'inputs'):
                raise ValueError('Input outside the selected run: ' + str(path))
            row = identity(path)
            if expected is not None and row['bytes'] != expected:
                raise ValueError('Restored input size changed: ' + str(path))
            found[str(path)] = {**row, 'reason': 'restored_input_with_retained_nearline_archive'}
    for directory, dirs, files in os.walk(root / 'work', followlinks=False):
        dirs[:] = [name for name in dirs if not (Path(directory) / name).is_symlink()]
        for name in files:
            path = Path(directory) / name
            if path.is_symlink() or not name.endswith(('.bam', '.cram', '.ubam', '.fastq', '.fq', '.fastq.gz', '.fq.gz')):
                continue
            row = identity(path)
            if row['bytes'] >= 1024**3:
                found[str(path)] = {**row, 'reason': 'completed_run_bulk_intermediate'}
    return list(found.values())


def preserve_results(root, release, libraries):
    results = root / 'results'
    if not any(results.glob('*/orange/*.orange.pdf')) or not (results / 'multiqc_report.html').is_file():
        raise ValueError('Native report and MultiQC must be present')
    alignments = list(results.glob('*/alignments/*/*.bam')) + list(results.glob('*/alignments/*/*.cram'))
    if not alignments or any(p.is_symlink() for p in alignments):
        raise ValueError('Published final alignments must be independently retained regular files')
    if {row['sample_id'] for row in libraries} - {p.parent.name for p in alignments}:
        raise ValueError('A source library has no retained final alignment')
    release.mkdir(parents=True, exist_ok=True)
    shutil.copytree(results, release / 'results', dirs_exist_ok=True,
                    ignore=lambda directory, names: ['alignments'] if 'alignments' in names else [])
    return [identity(path) for path in alignments]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, required=True)
    parser.add_argument('--input-plan', type=Path, required=True)
    parser.add_argument('--head-job', required=True)
    parser.add_argument('--release-dir', type=Path, required=True)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--apply', action='store_true')
    args = parser.parse_args()
    root = args.run_dir.resolve(strict=True)
    completed(args.head_job)
    inputs = json.loads(args.input_plan.read_text())
    if Path(inputs['run_directory']).resolve() != root:
        raise ValueError('Input plan names another run')
    retained = preserve_results(root, args.release_dir, inputs['libraries'])
    if args.plan.exists():
        plan = json.loads(args.plan.read_text())
        if plan['run_directory'] != str(root) or str(plan['head_job']) != args.head_job:
            raise ValueError('Cleanup plan identifies another execution')
    else:
        plan = dict(version=VERSION, run_directory=str(root), head_job=args.head_job, created_at=datetime.now(timezone.utc).isoformat(),
                    input_plan=str(args.input_plan), retained_alignments=retained, files=candidates(root, inputs))
        args.plan.parent.mkdir(parents=True, exist_ok=True)
        args.plan.write_text(json.dumps(plan, indent=2) + '\n')
    print(json.dumps(dict(files=len(plan['files']), bytes=sum(r['bytes'] for r in plan['files']), apply=args.apply)), flush=True)
    if not args.apply:
        return
    journal = args.plan.with_suffix('.journal.jsonl')
    removed = []
    already_absent = []
    with journal.open('a') as handle:
        for row in plan['files']:
            path = Path(row['path'])
            if not path.is_relative_to(root / 'inputs') and not path.is_relative_to(root / 'work'):
                raise ValueError('Cleanup path outside the permitted run directories')
            if not path.exists():
                already_absent.append(str(path))
                continue
            observed = identity(path)
            if any(observed[k] != row[k] for k in ('bytes', 'device', 'inode', 'mtime_ns')):
                raise ValueError('Bulk file changed since planning: ' + str(path))
            path.unlink()
            removed.append(row)
            handle.write(json.dumps(dict(path=str(path), bytes=row['bytes'], removed_at=datetime.now(timezone.utc).isoformat())) + '\n')
            handle.flush()
    result = dict(completed_at=datetime.now(timezone.utc).isoformat(), plan=str(args.plan),
                  removed_files=len(removed), removed_logical_bytes=sum(r['bytes'] for r in removed),
                  already_absent=already_absent,
                  retained_alignments=retained, retained_results=str(args.release_dir / 'results'))
    args.plan.with_suffix('.complete.json').write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
