#!/usr/bin/env bash
#SBATCH --job-name=rna-resource-canary
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:20:00
set -eo pipefail
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
module load nextflow/26.04.4
set -u
repo=$(readlink -f "${1:?repository required}")
out=$(readlink -m "${2:?new output directory required}")
mkdir -p "$out"
python3 - "$repo" "$out" <<'PY'
import importlib.util,json,pathlib,sys
repo,out=map(pathlib.Path,sys.argv[1:]);inputs=out/'inputs';inputs.mkdir(exist_ok=False)
spec=importlib.util.spec_from_file_location('res',repo/'workflow/build_input_resources.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
expected={}
for sample,sizes in [('S',[6,100,40,50]),('L',[63,1633,167,299])]:
    for name,gib in zip(['R1','R2','unsorted.bam','sorted.bam','md.bam'],[sizes[0]/2,sizes[0]/2,*sizes[1:]]):
        with (inputs/(sample+'.'+name)).open('wb') as f:f.truncate(int(gib*1024**3))
    (inputs/(sample+'.bai')).write_bytes(b'x')
    for stage,gib in zip(m.RNA_RESOURCES,sizes):
        expected[sample+'.'+stage+'.tsv']=m.rna_allocation(stage,int(gib*1024**3))
(out/'expected.json').write_text(json.dumps(expected,indent=2)+'\n')
(out/'rna.config').write_text(m.render_rna_config())
(out/'canary.config').write_text("process.executor='local'\nprocess.maxForks=1\nexecutor.cpus=128\nexecutor.memory=1024.GB\n")
PY
cd "$out"
export NXF_OPTS='-Xms1g -Xmx4g' NXF_ANSI_LOG=false
nextflow run "$repo/tests/test_rna_resource_config.nf" -c "$out/rna.config" -c "$out/canary.config" \
    --inputs "$out/inputs" --outdir "$out/results" -work-dir "$out/work" > "$out/nextflow.log" 2>&1
python3 - "$out" <<'PY'
import json,pathlib,sys
out=pathlib.Path(sys.argv[1]);expected=json.loads((out/'expected.json').read_text());observed={}
for name,spec in expected.items():
    found=list(map(int,(out/'results'/name).read_text().split()))
    wanted=[spec[k] for k in ['cpus','memory_gb','time_hours']]
    assert found==wanted,(name,found,wanted)
    observed[name]=found
(out/'completion.json').write_text(json.dumps(dict(status='PASS',scope='metadata-only task directive resolution under Slurm',observed=observed),indent=2)+'\n')
for p in (out/'inputs').iterdir():p.unlink()  # Only this canary's sparse synthetic files.
print(json.dumps(observed))
PY
