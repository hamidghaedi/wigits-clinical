#!/usr/bin/env python3
"""Add an operational completed-output subscription to an isolated native source.

No native analysis module/configuration is changed. Cache one source checkout
per exact base commit and hook content. Original sources remain untouched.
"""
import argparse
import fcntl
import hashlib
import json
from pathlib import Path
import subprocess

HERE=Path(__file__).resolve().parent
ANCHOR='    NFCORE_ONCOANALYSER()\n'
HOOK='''
    // Completed BAM-stage output; no timer and no caller-completion dependency.
    if (System.getenv('WIGITS_LAB_QC_CONFIG')) {
        channel.topic('bamtools_metrics_dir').subscribe { record ->
            LabQcEvents.bamReady(record[1], log)
        }
    }
'''


def _prepare(source, destination_root):
    source=Path(source).resolve();destination_root=Path(destination_root)
    binding=source/'lab_qc_event_source.json'
    if binding.exists():return json.loads(binding.read_text())
    base=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()
    if subprocess.check_output(['git','-C',str(source),'status','--porcelain'],text=True).strip():
        raise ValueError('Use a clean released native source')
    helper=HERE/'overlays/lab-qc-events/lib/LabQcEvents.groovy'
    digest=hashlib.sha256((HOOK+helper.read_text()).encode()).hexdigest()
    target=destination_root/(base[:12]+'-lab-qc-'+digest[:12])
    if (target/'lab_qc_event_source.json').exists():return json.loads((target/'lab_qc_event_source.json').read_text())
    if target.exists():raise ValueError('Incomplete source preparation exists: '+str(target))
    destination_root.mkdir(parents=True,exist_ok=True)
    subprocess.run(['git','clone','--quiet','--shared','--no-checkout',str(source),str(target)],check=True)
    subprocess.run(['git','-C',str(target),'checkout','--quiet','--detach',base],check=True)
    main=target/'main.nf';text=main.read_text()
    if text.count(ANCHOR)!=1:raise ValueError('Native main workflow anchor is ambiguous')
    main.write_text(text.replace(ANCHOR,ANCHOR+HOOK))
    (target/'lib').mkdir(exist_ok=True)
    (target/'lib/LabQcEvents.groovy').write_text(helper.read_text())
    subprocess.run(['git','-C',str(target),'add','main.nf','lib/LabQcEvents.groovy'],check=True)
    subprocess.run(['git','-C',str(target),'diff','--cached','--check'],check=True)
    subprocess.run(['git','-C',str(target),'commit','--quiet','-m','Notify laboratory QC from completed BAMTOOLS output'],check=True)
    commit=subprocess.check_output(['git','-C',str(target),'rev-parse','HEAD'],text=True).strip()
    result=dict(source=str(target),commit=commit,base_source=str(source),base_commit=base,hook_sha256=digest)
    # Ignore the local binding receipt; it is runtime metadata, not source.
    with (target/'.git/info/exclude').open('a') as f:f.write('\n/lab_qc_event_source.json\n')
    (target/'lab_qc_event_source.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def prepare(source,destination_root):
    root=Path(destination_root);root.mkdir(parents=True,exist_ok=True)
    with (root/'.lab_qc_source.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        return _prepare(source,root)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',required=True);p.add_argument('--source-root',required=True)
    a=p.parse_args();print(json.dumps(prepare(a.source,a.source_root)))

if __name__=='__main__':main()
