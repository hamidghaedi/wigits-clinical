#!/usr/bin/env python3
"""Submit LAB_RUN_QC once per complete laboratory sequencing run, then stage it.

Called only by native completed-BAMTOOLS output events, never a timer. The
laboratory CSV defines complete DNA run membership; submitted cases never
shrink it. Consume actual completed task paths, not asynchronously published
copies. No SAGE/PURPLE/ORANGE or analytical-head completion dependency.
Read task-local versions and reference links, not final combined versions.
No BAM/FASTQ content reads, new fingerprints, hashes of raw data, or job changes.
Save submission intent before sbatch; ambiguous submissions are never repeated.
Existing externally managed QC is adopted explicitly and never resubmitted.
"""
import argparse
import csv
from datetime import datetime, timezone
from zoneinfo import ZoneInfo
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import tarfile

import yaml

VERSION = '2.2.1'
PHASES = {'Phase 1':'Phase_1','Phase_1':'Phase_1','Phase 1.5':'phase_1.5',
          'phase_1.5':'phase_1.5','Phase 2':'Phase_2','Phase_2':'Phase_2'}


def now():
    return datetime.now(timezone.utc).isoformat()


def save(path, data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(data,indent=2)+'\n');tmp.replace(path)


def read(path):
    return json.loads(Path(path).read_text())


def safe(value):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9._-]*',value):
        raise ValueError('Unsafe laboratory identifier: '+value)
    return value


def queues(config):
    tracker=read(config['native_tracker_config'])
    result={q['path']:q for q in tracker.get('queues',[])}
    for root in config.get('queue_roots',[]):
        for pointer in sorted(Path(root).glob('*/lab_qc_source.json')):
            q=read(pointer);result[q['path']]=q
    return list(result.values())


def candidates(config):
    result={}
    tracker=read(config['native_tracker_config'])
    for old in tracker.get('historical_runs',[]):
        package=Path(old['package']);plan=read(package/'input_plan.json');release=read(package/'release.json')
        result[plan['run_directory']]=dict(case=dict(package=str(package)),plan=plan,release=release,
                                          job=str(old['head_job']),laboratory=None)
    for q in queues(config):
        for case in read(q['path'])['cases']:
            state=Path(q['state_directory'])/(case['case_id']+'.json')
            if not state.exists() or not read(state).get('analysis_job'): continue
            package=Path(case['package'])
            if not (package/'input_plan.json').is_file(): continue
            plan=read(package/'input_plan.json');release=read(package/'release.json')
            if not release.get('native_source_version'): continue
            result[plan['run_directory']]=dict(case=case,plan=plan,release=release,
                               job=str(read(state)['analysis_job']),laboratory=q.get('laboratory'))
    return list(result.values())


def libraries(item):
    bound=item['plan']['libraries']
    if isinstance(bound,dict):
        return [dict(sample_id=x['sample_submitter_id'],run_id=Path(x['fastq_parent_dir']).parent.name)
                for x in bound.values() if x.get('analyte_type')=='DNA']
    return [dict(sample_id=x['sample_id'],run_id=x['run_id']) for x in bound if x['role'] in ('tumor','normal')]


def laboratory_rows(config, items):
    sources=set(config.get('laboratories',[]))
    for item in items:
        if item.get('laboratory'): sources.add(item['laboratory'])
        sources.update(p for p in item['release'].get('sources',{}) if p.endswith('.csv'))
    rows=[]
    for source in sorted(sources):
        with Path(source).open(newline='',encoding='utf-8-sig') as f:
            reader=csv.DictReader(f)
            if not {'sample_submitter_id','fastq_parent_dir','analyte_type','lanes','phase'} <= set(reader.fieldnames or []):continue
            for row in reader:
                if row['analyte_type']=='DNA':
                    rows.append(dict(row,source=source,run_id=Path(row['fastq_parent_dir']).parent.name))
    return rows


def task_data(run):
    latest={}
    for trace in sorted((Path(run)/'results/pipeline_info').glob('execution_trace_*.txt')):
        with trace.open() as f:
            for row in csv.DictReader(f,delimiter='\t'): latest[row['name']]=row
    return latest


def producing_task(tasks, name, sample):
    rows=[r for n,r in tasks.items() if n.split(' (')[0].split(':')[-1]==name
          and (n.endswith('_'+sample+')') or n.endswith('('+sample+')'))]
    if len(rows)!=1 or rows[0]['status'] not in ('COMPLETED','CACHED') or str(rows[0]['exit'])!='0':
        raise ValueError(sample+': '+name+' not complete')
    return rows[0]


def task_work(run, task):
    options=list((Path(run)/'work').glob(task['hash']+'*'))
    if len(options)!=1: raise ValueError('Producing task directory unavailable: '+task['hash'])
    return options[0]


def completed_file(published, work, basename):
    sources=[work/basename,*work.glob('*/'+basename)]
    sources=[p for p in sources if p.is_file()]
    if not published.is_file() or published.stat().st_size==0 or len(sources)!=1:
        raise ValueError('Published producer output pending: '+str(published))
    if published.stat().st_size!=sources[0].stat().st_size:
        raise ValueError('Published producer output copy incomplete: '+str(published))
    return published


def sample_manifest(item, lab, config, folder):
    sid=safe(lab['sample_submitter_id']);case=item['plan']['case_id'];run=Path(item['plan']['run_directory'])
    event_file=Path(config['state_directory'])/'bam_events'/run.name/(sid+'.json')
    if config.get('event_driven'):
        if not event_file.is_file():raise ValueError(sid+': completed BAM event pending')
        event=read(event_file)
        if event['native_job']!=item['job']:raise ValueError(sid+': BAM event belongs to another native attempt')
        bam=Path(event['bam']);bai=Path(event['bai']);flags=Path(event['flags'])
        redux=Path(event['redux_work']);metrics=Path(event['bamtools_work'])
        if any(not p.is_file() or p.stat().st_size==0 for p in (bam,bai,flags)):
            raise ValueError(sid+': completed BAM event output unavailable')
    else:
        # Retained only for explicit historical metadata reconstruction/tests.
        tasks=task_data(run)
        redux=task_work(run,producing_task(tasks,'REDUX',sid))
        metrics=task_work(run,producing_task(tasks,'BAMTOOLS',sid))
        bam=run/'results'/case/'alignments'/sid/(sid+'.redux.bam');bai=Path(str(bam)+'.bai')
        completed_file(bam,redux,bam.name);completed_file(bai,redux,bai.name)
        flags=run/'results'/case/'bamtools'/sid/(sid+'.bam_metric.flag_counts.tsv')
        completed_file(flags,metrics,flags.name)
    command=(metrics/'.command.sh').read_text().replace('\\\n',' ')
    tokens=shlex.split(command);reference=(metrics/tokens[tokens.index('-ref_genome')+1]).resolve()
    profiles=[p for p in config['reference_profiles'] if Path(p['reference']).resolve()==reference]
    if len(profiles)!=1: raise ValueError(sid+': exact producing reference profile is not registered')
    if any(not Path(profiles[0][k]).is_file() for k in ('reference','fai','dict')):
        raise ValueError(sid+': producing reference assets are unavailable')
    versions={};version_sources=[]
    for kind,work in [('REDUX',redux),('BAMTOOLS',metrics)]:
        version_file=work/'versions.yml'
        native=yaml.safe_load(version_file.read_text())
        values=[v for k,v in native.items() if k.split(':')[-1]==kind]
        if len(values)!=1: raise ValueError(sid+': task-emitted version metadata is ambiguous')
        versions[kind]=values[0];version_sources.append(str(version_file))
    version_path=folder/(sid+'.source_versions.json');save(version_path,versions)
    fq=[];observed_lanes=[]
    with (Path(item['case']['package'])/'samplesheet.csv').open() as f:
        for row in csv.DictReader(f):
            if row['sample_id']!=sid or row['sequence_type']!='dna':continue
            if row['filetype']!='fastq':raise ValueError(sid+': original FASTQ binding is unavailable')
            paths=[Path(x) for x in row['filepath'].split(';')]
            if len(paths)!=2 or any(not p.is_file() or p.stat().st_size==0 for p in paths):
                raise ValueError(sid+': complete FASTQ pair pending')
            info=dict(x.split(':',1) for x in row['info'].split(';') if ':' in x)
            if not info.get('lane'):raise ValueError(sid+': explicit lane metadata missing')
            lane=int(info['lane']);observed_lanes.append(lane)
            match=re.search(r'_(S\d+_L\d{3})_R1_',paths[0].name)
            label=match[1] if match else 'L'+str(lane).zfill(3)
            if any(x['label']==label for x in fq):
                raise ValueError(sid+': multiple chunks require an explicit released QC binding')
            fq.append(dict(label=label,r1=str(paths[0]),r2=str(paths[1]),
                           input_gib=sum(p.stat().st_size for p in paths)/2**30))
    if not fq or set(observed_lanes)!={int(x) for x in lab['lanes'].split(';')}:
        raise ValueError(sid+': not all laboratory lanes are bound')
    # A combined multi-run BAM cannot be assigned to a single sequencing run.
    objects=item['plan'].get('objects',[])
    seqruns={x['run_id'] for x in objects if x.get('sample_id')==sid}
    if seqruns and seqruns!={lab['run_id']}:raise ValueError(sid+': multi-run BAM needs per-run metric binding')
    summary=flags.with_name(sid+'.bam_metric.summary.tsv')
    return dict(run_id=lab['run_id'],sample_id=sid,bam=str(bam),bai=str(bai),bam_gib=bam.stat().st_size/2**30,
                bamtools_summary=str(summary) if summary.is_file() else None,
                bamtools_flags=str(flags),source_versions=str(version_path),fastqs=fq,
                source_case_package=item['case']['package'],source_native_job=item['job'],
                source_version_files=version_sources),profiles[0]


def build_manifest(run_id, labs, items, config, folder):
    expected={}
    for row in labs:
        if row['run_id']!=run_id:continue
        sid=row['sample_submitter_id'];phase=PHASES[row['phase']]
        if sid in expected and (expected[sid]['lanes']!=row['lanes'] or PHASES[expected[sid]['phase']]!=phase):
            raise ValueError('Conflicting laboratory scope: '+sid)
        expected[sid]=row
    if not expected:raise ValueError('Complete laboratory run membership unavailable')
    phases={PHASES[r['phase']] for r in expected.values()}
    if len(phases)!=1:raise ValueError('Sequencing-run phase is ambiguous')
    samples=[];profiles=[]
    for sid,lab in sorted(expected.items()):
        matches=[i for i in items if dict(sample_id=sid,run_id=run_id) in libraries(i)]
        if len(matches)!=1:raise ValueError(sid+': native BAM binding pending or ambiguous')
        sample,profile=sample_manifest(matches[0],lab,config,folder)
        samples.append(sample);profiles.append(profile)
    if len({p['reference'] for p in profiles})!=1:raise ValueError('Run contains different alignment references')
    def values(key):return '; '.join(sorted({r.get(key,'') for r in expected.values()}-{''}))
    metadata=dict(project_id=values('project_id'),batch_id=values('batch_id'),flowcell_id=values('flowcell_id'),
                  run_date=values('run_date'),tissue_type=values('tissue_preservation_method'),
                  report_date=datetime.now(ZoneInfo('America/Toronto')).date().isoformat())
    return dict(reference=profiles[0]['reference'],fai=profiles[0]['fai'],dict=profiles[0]['dict'],
                logo=config['logo'],css=config['css'],samples=samples,
                expected_samples={run_id:sorted(expected)},run_metadata={run_id:metadata},
                coverage_producer='parabricks.bammetrics',demultiplexing={run_id:config.get('demultiplexing',{}).get(run_id,[])},laboratory_sources=sorted({r['source'] for r in expected.values()})),phases.pop()



def deliver(config, state, result_dir=None):
    from_root=Path(result_dir) if result_dir else Path(state['execution_root'])/'results'/state['run_id']
    success=from_root/'_SUCCESS'
    if not success.is_file():return False
    receipt=read(success);report=from_root/(state['run_id']+'.multiqc_report.html')
    if receipt['run_id']!=state['run_id'] or hashlib.sha256(report.read_bytes()).hexdigest()!=receipt['report_sha256']:
        raise ValueError('Produced run report differs from its success receipt')
    target=Path(config['archive_root'])/state['phase']/'qc_reports'/state['run_id']
    evidence=target/'evidence'/state['analysis_id'];evidence.mkdir(parents=True,exist_ok=True)
    archive=evidence/'qc_sources.tar.gz'
    if not archive.exists():
        tmp=archive.with_suffix('.tmp')
        with tarfile.open(tmp,'w:gz') as tar:
            for name in ['qc','provenance','_SUCCESS',state['run_id']+'.multiqc_report_data.zip']:
                tar.add(from_root/name,arcname=name)
        tmp.replace(archive)
    manifest=Path(state['execution_root'])/'delivery.json'
    save(manifest,dict(reports=[dict(phase=state['phase'],run_id=state['run_id'],analysis_id=state['analysis_id'],
        source=str(report),filename=state['run_id']+'.laboratory.multiqc_report.html')]))
    result=subprocess.run([config['python'],config['stager'],'--manifest',str(manifest),
                           '--local-root',config['archive_root']],check=True,capture_output=True,text=True)
    save(Path(state['execution_root'])/'delivery_receipt.json',json.loads(result.stdout))
    return True


def job_states(ids):
    if not ids:return {}
    result=subprocess.run(['sacct','-X','-n','-j',','.join(sorted(ids)),
        '--format=JobIDRaw,State,ExitCode','-P'],check=True,capture_output=True,text=True)
    return {v[0]:v[1:3] for line in result.stdout.splitlines() if len(v:=line.split('|'))>=3}


def tick(config, runner=subprocess.run, status_reader=None, only_run=None):
    root=Path(config['state_directory']);root.mkdir(parents=True,exist_ok=True)
    items=candidates(config);labs=laboratory_rows(config,items)
    eligible={l['run_id'] for i in items if i['job'] not in config['baseline_job_ids'] for l in libraries(i)}
    states={read(p)['run_id']:read(p) for p in root.glob('*.run.json')}
    active={s['job_id'] for s in states.values() if s.get('job_id') and not s.get('external_managed') and s.get('status')!='delivered'}
    jobs=status_reader(active) if status_reader else {}
    report=[]
    for run_id in sorted(eligible|set(states)):
        if only_run is not None and run_id!=only_run:continue
        safe(run_id);path=root/(run_id+'.run.json');state=states.get(run_id,dict(run_id=run_id))
        if state.get('external_managed'):
            report.append(dict(run_id=run_id,status='external_managed',job_id=state['job_id']));continue
        if state.get('job_id'):
            if state.get('status')!='delivered' and deliver(config,state):
                state.update(status='delivered',delivered_at=now());save(path,state)
            elif state.get('status')!='delivered' and state['job_id'] in jobs:
                jobstate,exitcode=jobs[state['job_id']];state.update(scheduler_state=jobstate,scheduler_exit=exitcode)
                if jobstate.split()[0] in ('FAILED','TIMEOUT','CANCELLED','OUT_OF_MEMORY','NODE_FAIL','PREEMPTED'):
                    state['status']='execution_failed'
                elif jobstate=='COMPLETED':state['status']='completion_receipt_missing'
                save(path,state)
            report.append(dict(run_id=run_id,status=state['status'],job_id=state['job_id']));continue
        if state.get('submission_intent'):
            report.append(dict(run_id=run_id,status='submission_unconfirmed'));continue
        folder=Path(config['run_root'])/('lab-qc-'+run_id+'-'+config['workflow_version'].replace('/','-'))
        folder.mkdir(parents=True,exist_ok=True)
        try:
            manifest,phase=build_manifest(run_id,labs,items,config,folder/'sources')
        except (ValueError,KeyError,FileNotFoundError) as exc:
            state.update(status='waiting_for_bams',reason=str(exc),observed_at=now());save(path,state)
            report.append(dict(run_id=run_id,status=state['status'],reason=state['reason']));continue
        manifest_path=folder/'manifest.json';save(manifest_path,manifest);(folder/'logs').mkdir(exist_ok=True)
        command=['sbatch','--parsable','--account='+config['account'],
                 '--output='+str(folder/'logs/lab_head_%j.log'),
                 '--export=ALL,WIGITS_WORKFLOW_DIR='+config['workflow_directory']+
                 ',WIGITS_LAB_QC_CONFIG='+config.get('config_path',''),
                 str(Path(config['workflow_directory'])/'run_laboratory_qc.sbatch'),str(manifest_path),str(folder)]
        state.update(status='submission_unconfirmed',phase=phase,submission_intent=True,command=command,
                     execution_root=str(folder),analysis_id=folder.name,submitted_at=now(),
                     workflow_version=config['workflow_version'],workflow_commit=config['workflow_commit'])
        save(path,state)
        result=runner(command,check=True,capture_output=True,text=True,
                      env={k:v for k,v in os.environ.items() if k!='PYTHONPATH'})
        job=result.stdout.strip().split(';')[0]
        if not job.isdigit():raise ValueError('Unconfirmed sbatch response; saved intent prevents duplicate submission')
        state.update(status='submitted',job_id=job,submission_intent=False);save(path,state)
        report.append(dict(run_id=run_id,status='submitted',job_id=job))
    save(root/'last_tick.json',dict(version=VERSION,observed_at=now(),runs=report))
    return report


def bam_ready(config, metrics_dir, runner=subprocess.run):
    """One native completed-output notification; persist it before considering QC."""
    metrics_dir=Path(metrics_dir).resolve();metrics=metrics_dir.parent
    if (metrics/'.exitcode').read_text().strip()!='0':
        raise ValueError('BAMTOOLS event is not a successful completed task')
    flags=list(metrics_dir.glob('*.bam_metric.flag_counts.tsv'))
    if len(flags)!=1:raise ValueError('BAMTOOLS event sample is ambiguous')
    sid=flags[0].name.removesuffix('.bam_metric.flag_counts.tsv');safe(sid)
    matches=[i for i in candidates(config) if metrics.is_relative_to(Path(i['plan']['run_directory'])/'work')]
    if len(matches)!=1:raise ValueError('Completed BAM event is not bound to one registered native run')
    item=matches[0];run=Path(item['plan']['run_directory'])
    bound=[l for l in libraries(item) if l['sample_id']==sid]
    if len(bound)!=1:raise ValueError('Completed BAM event has no unique laboratory library binding')
    tokens=shlex.split((metrics/'.command.sh').read_text().replace(chr(92)+chr(10),' '))
    bam=(metrics/tokens[tokens.index('-bam_file')+1]).resolve();bai=Path(str(bam)+'.bai')
    redux=bam.parent.parent
    if not (redux/'versions.yml').is_file():raise ValueError('Native REDUX version source unavailable')
    event=dict(sample_id=sid,run_id=bound[0]['run_id'],native_job=item['job'],received_at=now(),
               bam=str(bam),bai=str(bai),flags=str(flags[0]),redux_work=str(redux),bamtools_work=str(metrics))
    save(Path(config['state_directory'])/'bam_events'/run.name/(sid+'.json'),event)
    return tick(dict(config,event_driven=True),runner,only_run=event['run_id'])


def finish(config, run_id, result_dir=None):
    path=Path(config['state_directory'])/(safe(run_id)+'.run.json');state=read(path)
    if state.get('external_managed') or state.get('status')=='delivered':return
    if result_dir and not Path(result_dir).resolve().is_relative_to(Path(state['execution_root']).resolve()):
        raise ValueError('Report output is outside the registered QC execution')
    if not deliver(config,state,result_dir):raise ValueError('Laboratory QC completion receipt unavailable')
    state.update(status='delivered',delivered_at=now(),delivered_from=str(result_dir) if result_dir else 'published outputs');save(path,state)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--config',required=True)
    events=parser.add_mutually_exclusive_group(required=True)
    events.add_argument('--bam-ready',help='Completed native BAMTOOLS output directory')
    events.add_argument('--deliver-run',help='Successful LAB_RUN_QC completion event')
    events.add_argument('--report-ready',help='Completed run report task success receipt')
    args=parser.parse_args();config=read(args.config);config['config_path']=str(Path(args.config).resolve())
    root=Path(config['state_directory']);root.mkdir(parents=True,exist_ok=True)
    with (root/'dispatch.lock').open('a') as lock:
        # Events wait briefly for one another; a nonblocking lock could lose the last BAM event.
        fcntl.flock(lock,fcntl.LOCK_EX)
        if args.bam_ready:bam_ready(config,args.bam_ready)
        elif args.report_ready:
            success=Path(args.report_ready).resolve();finish(config,read(success)['run_id'],success.parent)
        else:finish(config,args.deliver_run)


if __name__=='__main__':main()
