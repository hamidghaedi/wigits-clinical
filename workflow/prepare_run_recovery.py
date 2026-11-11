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

VERSION='0.1.0'


def resume_command(original, session, config):
    args=shlex.split(original)
    if args[:3]!=['nextflow','run','nf-core/oncoanalyser']:
        raise ValueError('Expected the native pinned launcher')
    i=args.index('-resume')
    if i+1<len(args) and not args[i+1].startswith('-'):
        args[i+1]=session
    else:
        args.insert(i+1,session)
    args+=['-c',str(config)]
    return shlex.join(args)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--route',choices=('bam','native'),required=True)
    p.add_argument('--session',required=True)
    p.add_argument('--patched-jar',type=Path)
    p.add_argument('--kraken-work',type=Path)
    p.add_argument('--tumor-id')
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    here=Path(__file__).resolve().parent
    config=a.out/'recovery.config'
    if a.route=='bam':
        env_file=a.out/'esvee_jvm_recovery.env'
        shutil.copyfile(here/'conf/esvee_jvm_recovery.env',env_file)
        config.write_text((here/'conf/bam_runtime_recovery.config').read_text().replace('@ESVEE_ENV_FILE@',str(env_file.resolve())))
    else:
        if not all((a.patched_jar,a.kraken_work,a.tumor_id)):
            p.error('Native recovery requires patched jar, Kraken work directory and tumour ID')
        prefix=a.tumor_id+'.virusbreakend.vcf'
        destination=prefix+'.virusbreakend.working'
        commands=['module load apptainer/1.3.5','export TMPDIR="${SLURM_TMPDIR:-$PWD}"',
                  'mkdir -p '+shlex.quote(destination)]
        for suffix in ('readnames.txt','kraken2.report.all.txt','viralreferencescontigs.txt'):
            f=a.kraken_work/(prefix+'.'+suffix)
            if not f.is_file():raise FileNotFoundError(f)
            commands.append('ln -sfn '+shlex.quote(str(f))+' '+shlex.quote(destination+'/'+f.name))
        script='\n'.join(commands)
        config.write_text("// Recovery v0.1.0: preserve counts and reuse completed Kraken scan.\nprocess {\n"
            "    withName: 'VIRUSBREAKEND' {\n"
            f"        ext.args = '--jar {a.patched_jar}'\n"
            f"        beforeScript = '''\n{script}\n'''\n    }}\n}}\n")
    original=[l for l in (a.run/'launch.sh').read_text().splitlines() if l.startswith('nextflow run ')][-1]
    (a.out/'resume_command.sh').write_text('#!/usr/bin/env bash\nset -euo pipefail\nexec '+resume_command(original,a.session,config.resolve())+'\n')
    (a.out/'run_directory.txt').write_text(str(a.run.resolve())+'\n')
    shutil.copyfile(here/'resume_completed_run.sbatch',a.out/'resume_completed_run.sbatch')
    evidence=dict(version=VERSION,route=a.route,session=a.session,
        original_command=original,original_samplesheet_sha256=hashlib.sha256((a.run/'samplesheet.csv').read_bytes()).hexdigest(),
        files={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in a.out.iterdir() if f.is_file()})
    (a.out/'release.json').write_text(json.dumps(evidence,indent=2)+'\n')


if __name__=='__main__':main()
