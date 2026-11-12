#!/usr/bin/env python3
"""Read-only full gzip/FASTQ/mate check for existing paired input files.

Run under Slurm. Pigz checks gzip CRC/truncation while Python checks every
four-line record and paired read name. This establishes file readability and
pair consistency, not agreement with an unavailable original archive checksum.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import tempfile


def read_name(header):
    name = header.split()[0]
    return name[:-2] if name.endswith((b'/1', b'/2')) else name


def check_pair(sample, paths):
    paths = [Path(p).resolve() for p in paths]
    before = [p.stat() for p in paths]
    processes = []
    with tempfile.TemporaryFile() as errors:
        try:
            for path in paths:
                processes.append(subprocess.Popen(
                    ['pigz', '-dc', '-p', '1', str(path)],
                    stdout=subprocess.PIPE, stderr=errors, bufsize=1024*1024))
            pairs = 0
            bases = [0, 0]
            while True:
                headers = [p.stdout.readline() for p in processes]
                if headers == [b'', b'']:
                    break
                if not all(headers):
                    raise ValueError(f'{sample}: unequal mate counts after {pairs} pairs')
                for i, process in enumerate(processes):
                    seq, plus, qual = [process.stdout.readline().rstrip(b'\r\n') for _ in range(3)]
                    if not headers[i].startswith(b'@') or not plus.startswith(b'+') or not seq or len(seq) != len(qual):
                        raise ValueError(f'{sample}: malformed mate {i+1} at record {pairs+1}')
                    bases[i] += len(seq)
                if read_name(headers[0]) != read_name(headers[1]):
                    raise ValueError(f'{sample}: different mate names at record {pairs+1}')
                pairs += 1
                if pairs % 10_000_000 == 0:
                    print(f'{sample}: {pairs:,} pairs checked', flush=True)
            codes = [p.wait() for p in processes]
            if any(codes):
                errors.seek(0)
                raise ValueError(f'{sample}: gzip integrity failed: {errors.read().decode(errors="replace")}')
            if not pairs:
                raise ValueError(f'{sample}: empty FASTQs')
            after = [p.stat() for p in paths]
            if any((a.st_ino, a.st_size, a.st_mtime_ns) != (b.st_ino, b.st_size, b.st_mtime_ns)
                   for a, b in zip(before, after)):
                raise ValueError(f'{sample}: input changed during scan')
            return dict(sample=sample, status='PASS', read_pairs=pairs, bases=bases,
                        files=[dict(path=str(p), bytes=s.st_size, mtime_ns=s.st_mtime_ns)
                               for p, s in zip(paths, after)])
        finally:
            for process in processes:
                if process.poll() is None:
                    process.terminate()
                if process.stdout:
                    process.stdout.close()
                process.wait()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samplesheet', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    with args.samplesheet.open() as handle:
        rows = list(csv.DictReader(handle))
    if not rows or any(r['filetype'] != 'fastq' or len(r['filepath'].split(';')) != 2 for r in rows):
        raise ValueError('Expected paired FASTQ rows')
    args.out.mkdir(parents=True, exist_ok=False)
    with ProcessPoolExecutor(max_workers=2) as pool:
        pending = [pool.submit(check_pair, r['sample_id'], r['filepath'].split(';')) for r in rows]
        results = []
        for row, future in zip(rows, pending):
            result = future.result()
            (args.out / (row['sample_id'] + '.json')).write_text(json.dumps(result, indent=2)+'\n')
            results.append(result)
    (args.out/'completion.json').write_text(json.dumps(dict(
        status='PASS', completed_utc=datetime.now(timezone.utc).isoformat(), libraries=results,
        checks=['full gzip decompression and CRC', 'FASTQ structure', 'all paired read names/counts'],
        archive_checksum_comparison='not performed; no bound original checksum supplied'), indent=2)+'\n')
    print('All input pairs PASS', flush=True)


if __name__ == '__main__':
    main()
