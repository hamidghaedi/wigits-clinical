#!/usr/bin/env python3
"""Transfer an exact released file list, one FASTQ through two hops at a time.

Only Globus metadata/API calls run here. Globus performs data transport. Resume
saved task IDs rather than submitting duplicates. Never delete payloads or
dispatch analysis. Both hops request checksum verification; destination sizes
are checked locally only after successful terminal transfer status.
"""
import argparse
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import subprocess
import time


def globus(*args):
    return json.loads(subprocess.check_output(['globus', *args, '--format', 'json'], text=True))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--state-dir', type=Path, required=True)
    a = p.parse_args()
    manifest = json.loads(a.manifest.read_text())
    a.state_dir.mkdir(parents=True, exist_ok=True)
    path = a.state_dir / 'state.json'
    with (a.state_dir/'controller.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        state = json.loads(path.read_text()) if path.exists() else {'files': {}}

        def save():
            state['updated_utc'] = datetime.now(timezone.utc).isoformat()
            tmp = path.with_suffix('.tmp')
            tmp.write_text(json.dumps(state, indent=2) + '\n')
            os.replace(tmp, path)

        try:
            state['status'] = 'ACTIVE'
            save()
            for i, row in enumerate(manifest['objects'], 1):
                file_state = state['files'].setdefault(row['source'], {})
                for hop, src_ep, dst_ep, src, dst in [
                    ('nas_to_fir', manifest['source_endpoint'], manifest['relay_endpoint'], row['source'], row['relay']),
                    ('fir_to_frontenac', manifest['relay_endpoint'], manifest['destination_endpoint'], row['relay'], row['destination']),
                ]:
                    if file_state.get(hop + '_complete'):
                        continue
                    key = hop + '_task'
                    if not file_state.get(key):
                        task = globus('transfer', src_ep + ':' + src, dst_ep + ':' + dst,
                                      '--sync-level', 'checksum', '--verify-checksum', '--preserve-mtime',
                                      '--label', f"{manifest['case']}-full-raw-{i:03d}-{hop}")
                        file_state[key] = task['task_id']
                        save()
                    task_id = file_state[key]
                    print(f'{i}/{len(manifest["objects"])} {Path(src).name} {hop} {task_id}', flush=True)
                    while True:
                        task = globus('task', 'show', task_id)
                        (a.state_dir/f'task_{task_id}.json').write_text(json.dumps(task, indent=2) + '\n')
                        if task['status'] == 'SUCCEEDED':
                            break
                        if task['status'] != 'ACTIVE':
                            raise RuntimeError(f'{task_id}: {task["status"]}')
                        time.sleep(30)
                    if task['source_endpoint_id'] != src_ep or task['destination_endpoint_id'] != dst_ep:
                        raise ValueError('Transfer endpoint mismatch')
                    if task['files_transferred'] + task.get('files_skipped', 0) != 1:
                        raise ValueError('Expected one transferred/checksum-matched file')
                    file_state[hop + '_complete'] = True
                    save()
                dest = Path(row['destination'])
                if not dest.is_file() or dest.is_symlink() or dest.stat().st_size != row['expected_bytes']:
                    raise ValueError(f'Destination identity/size mismatch: {dest}')
                file_state['status'] = 'VERIFIED_LOCAL'
                save()
            state['status'] = 'ALL_FASTQS_VERIFIED_LOCAL'
            save()
            print(state['status'], flush=True)
        except Exception as error:
            state['status'] = 'STOPPED_EXCEPTION'
            state['error'] = str(error)
            save()
            raise


if __name__ == '__main__':
    main()
