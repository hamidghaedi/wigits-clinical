#!/usr/bin/env python3
"""Freeze a resume command from the run's original launch, preserving its cache.

No new samplesheet, input-resource policy, native source or alignment entry.
Native VirusBreakend can reuse only its completed Kraken read-name/report files.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import shutil

VERSION='0.4.0'


def resume_command(original, session, config, native_source=None):
    args=shlex.split(original)
    allowed = ['nf-core/oncoanalyser'] + ([str(native_source)] if native_source else [])
    if args[:2]!=['nextflow','run'] or len(args)<3 or args[2] not in allowed:
        raise ValueError('Expected the native pinned launcher')
    i=args.index('-resume')
    if i+1<len(args) and not args[i+1].startswith('-'):
        args[i+1]=session
    else:
        args.insert(i+1,session)
    args+=['-c',str(config)]
    return shlex.join(args)


def kraken_config(work, tumor_id, patched_jar=None):
    prefix=tumor_id+'.virusbreakend.vcf'
    destination=prefix+'.virusbreakend.working'
    commands=['module load apptainer/1.3.5','export TMPDIR="${SLURM_TMPDIR:-$PWD}"',
              'mkdir -p '+shlex.quote(destination)]
    for suffix in ('readnames.txt','kraken2.report.all.txt','viralreferencescontigs.txt'):
        f=work/(prefix+'.'+suffix)
        if not f.is_file():raise FileNotFoundError(f)
        commands.append('ln -sfn '+shlex.quote(str(f))+' '+shlex.quote(destination+'/'+f.name))
    jar=f"        ext.args = '--jar {patched_jar}'\n" if patched_jar else ''
    return ("// Reuse completed Kraken scan; regenerate incomplete downstream products.\nprocess {\n"
            "    withName: 'VIRUSBREAKEND' {\n"+jar+
            "        beforeScript = '''\n"+'\n'.join(commands)+"\n'''\n    }\n}\n")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--route',choices=('bam','native','rna','esvee'),required=True)
    p.add_argument('--session',required=True)
    p.add_argument('--native-source',type=Path,
                   help='Explicit existing native source path from the saved release; must match the original launch')
    p.add_argument('--patched-jar',type=Path)
    p.add_argument('--kraken-work',type=Path)
    p.add_argument('--tumor-id')
    p.add_argument('--rna-resources', action='store_true',
                   help='Use measured task-input RNA resources on an RNA recovery; preserve DNA allocations')
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    if a.rna_resources and a.route != 'rna':
        p.error('--rna-resources requires --route rna')
    here=Path(__file__).resolve().parent
    config=a.out/'recovery.config'
    if a.route=='rna':
        config.write_text((here/'conf/rna_duplicate_marking.config').read_text())
        if a.rna_resources:
            from build_input_resources import render_rna_config
            config.write_text(config.read_text() + render_rna_config())
    elif a.route in ('bam','esvee'):
        env_file=a.out/'esvee_jvm_recovery.env'
        shutil.copyfile(here/'conf/esvee_jvm_recovery.env',env_file)
        template='esvee_runtime_recovery.config' if a.route=='esvee' else 'bam_runtime_recovery.config'
        config.write_text((here/'conf'/template).read_text().replace('@ESVEE_ENV_FILE@',str(env_file.resolve())))
        if a.route=='esvee' and (a.kraken_work or a.tumor_id):
            if not all((a.kraken_work,a.tumor_id)):
                p.error('Kraken reuse requires both --kraken-work and --tumor-id')
            config.write_text(config.read_text()+kraken_config(a.kraken_work,a.tumor_id))
    else:
        if not all((a.patched_jar,a.kraken_work,a.tumor_id)):
            p.error('Native recovery requires patched jar, Kraken work directory and tumour ID')
        config.write_text(kraken_config(a.kraken_work,a.tumor_id,a.patched_jar))
    original=[l for l in (a.run/'launch.sh').read_text().splitlines() if l.startswith('nextflow run ')][-1]
    (a.out/'resume_command.sh').write_text('#!/usr/bin/env bash\nset -euo pipefail\nexec '+resume_command(original,a.session,config.resolve(),a.native_source)+'\n')
    (a.out/'run_directory.txt').write_text(str(a.run.resolve())+'\n')
    shutil.copyfile(here/'resume_completed_run.sbatch',a.out/'resume_completed_run.sbatch')
    evidence=dict(version=VERSION,route=a.route,session=a.session,rna_resource_policy='input-v5' if a.rna_resources else None,
        original_command=original,original_samplesheet_sha256=hashlib.sha256((a.run/'samplesheet.csv').read_bytes()).hexdigest(),
        files={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in a.out.iterdir() if f.is_file()})
    (a.out/'release.json').write_text(json.dumps(evidence,indent=2)+'\n')


if __name__=='__main__':main()
