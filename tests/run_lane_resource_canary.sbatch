#!/usr/bin/env bash
#SBATCH --job-name=lane-resource-canary
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:20:00
# Metadata-only local tasks inside this head: no alignment, file scans or real data.
set -eo pipefail
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
module load nextflow/26.04.4
set -u
repo=$(readlink -f "${1:?repository required}")
out=$(readlink -m "${2:?new canary output directory required}")
mkdir -p "$out"
python3 - "$repo" "$out" <<'PY'
import csv, importlib.util, json, pathlib, subprocess, sys
repo, out = map(pathlib.Path, sys.argv[1:])
inputs = out / 'inputs'; inputs.mkdir(exist_ok=False)
rows = []
for sample, role, gib in [('N', 'normal', 5), ('T', 'tumor', 34)]:
    for lane in range(1, 9):
        paths = []
        for mate in (1, 2):
            name = f'{sample}{mate}.fastq.gz' if lane == 1 else f'{sample}_L{lane}_{mate}.fastq.gz'
            path = inputs / name
            with path.open('wb') as stream:
                stream.truncate(gib * 1024**3 // 2)  # Sparse: never read as FASTQ.
            paths.append(str(path))
        rows.append(dict(sample_id=sample, sample_type=role, sequence_type='dna',
                         filetype='fastq', filepath=';'.join(paths)))
with (out / 'samplesheet.csv').open('w') as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
subprocess.run([sys.executable, str(repo / 'workflow/build_input_resources.py'),
    '--samplesheet', str(out / 'samplesheet.csv'), '--out', str(out / 'resources.config'),
    '--evidence', str(out / 'resources.json')], check=True)
(out / 'canary.config').write_text('''
process.executor = 'local'
process.maxForks = 1
executor.cpus = 128
executor.memory = 1024.GB
''')
PY
cd "$out"
export NXF_OPTS='-Xms1g -Xmx4g'
export NXF_ANSI_LOG=false
nextflow run "$repo/tests/test_lane_resource_config.nf" \
    -c "$out/resources.config" -c "$out/canary.config" \
    --input_dir "$out/inputs" --outdir "$out/results" -work-dir "$out/work" \
    > "$out/nextflow_stdout.log" 2>&1
python3 - "$out" <<'PY'
import json, pathlib, sys
root = pathlib.Path(sys.argv[1]); plan = json.loads((root / 'resources.json').read_text())
expected = {}
for lane in plan['alignment_lanes']:
    expected[lane['sample_id'] + '.lane.tsv'] = lane
for sample, value in plan['redux_samples'].items():
    expected[sample + '.sample.tsv'] = value
observed = {}
for name, value in expected.items():
    sample, cpu, memory, hours = (root / 'results' / name).read_text().strip().split('\t')
    actual = [int(cpu), int(memory), int(hours)]
    wanted = [value[k] for k in ('cpus', 'memory_gb', 'time_hours')]
    assert actual == wanted, (name, actual, wanted)
    observed[name] = actual
(root / 'completion.json').write_text(json.dumps(dict(status='PASS', observed=observed,
    scope='real Nextflow directive resolution; metadata-only local tasks inside Slurm'), indent=2) + '\n')
for path in (root / 'inputs').iterdir():
    path.unlink()  # Only synthetic sparse test files created by this canary.
print(json.dumps(observed))
PY
