#!/usr/bin/env python3
"""Read-only Slurm/Nextflow resource report. No analysis or scheduler mutations.

One CSV row per native job ID, including unsuccessful attempts and controllers.
Completed CPU comes from the Slurm parent; memory is a labelled step proxy.
Live CPU is provisional batch-step accounting, never a zero from sacct.
"""
import argparse
import collections
import csv
from datetime import datetime
import html
import io
import json
from pathlib import Path
import re
import subprocess
from zoneinfo import ZoneInfo

VERSION = '0.1.0'
FIELDS = ['cluster','job_id','kind','step','task','attempt','state','exit_code',
          'allocated_cpus','allocated_memory_gib','elapsed_seconds','cpu_seconds',
          'average_busy_cpus','cpu_use_percent','cpu_measurement','allocated_core_hours',
          'observed_maxrss_gib','rss_fraction_percent','memory_source','queue_seconds',
          'allocated_gpus','gpu_card_hours','nextflow_cpu_percent','nextflow_peak_rss',
          'source','notes']
ACTIVE = {'RUNNING','COMPLETING','SUSPENDED'}


def seconds(value):
    if not value or value in {'Unknown','N/A','-'}:
        return None
    days, sep, clock = value.partition('-')
    parts = (clock if sep else value).split(':')
    return (float(days)*86400 if sep else 0) + sum(float(p)*60**i for i,p in enumerate(reversed(parts)))


def gib(value):
    match = re.fullmatch(r'([0-9.]+)([KMGTPE]?)(?:i?B)?', value or '', re.I)
    if not match:
        return None
    power = 'KMGTPE'.index(match[2].upper())+1 if match[2] else 0
    return float(match[1])*1024**power/1024**3


def records(text, delimiter='|'):
    return list(csv.DictReader(io.StringIO(text), delimiter=delimiter)) if text.strip() else []


def discover(root):
    jobs = {}
    def add(job, **details):
        job = str(job or '')
        if re.fullmatch(r'\d+(?:_\d+)?', job):
            target = jobs.setdefault(job, {})
            target.update({k:v for k,v in details.items() if v and not target.get(k)})
    # Bounded metadata locations only: never walk BAM/FASTQ payloads or work/.
    metadata = list((root/'state').glob('*execution*.json'))
    metadata += list((root/'state').glob('*/previous_execution.json'))
    metadata += list((root/'inputs').glob('*/dispatch.json'))
    def walk(value, source, label=''):
        if isinstance(value, dict):
            if 'job_id' in value:
                add(value['job_id'], source=source, task=value.get('sample', label))
            for key, child in value.items():
                walk(child, source, key)
        elif isinstance(value, list):
            for child in value:
                walk(child, source, label)
    for path in sorted(metadata):
        try:
            walk(json.loads(path.read_text()), str(path))
        except (ValueError, OSError):
            pass  # A concurrently written metadata file is optional discovery.
    for path in sorted((root/'logs').glob('*.log')):
        match = re.search(r'\.(\d+)\.log$', path.name)
        if match:
            add(match[1], source=str(path))
    for path in sorted((root/'results/pipeline_info').glob('*trace*.txt')):
        for row in records(path.read_text(), '\t'):
            add(row.get('native_id'), task=row.get('name'), attempt=row.get('attempt'),
                nextflow_cpu_percent=row.get('%cpu'), nextflow_peak_rss=row.get('peak_rss'),
                source=str(path))
    logs = list(root.glob('.nextflow.log*')) + list((root/'logs').glob('nextflow_head*.log'))
    for path in sorted(logs):
        for line in path.read_text(errors='replace').splitlines():
            match = re.search(r'\[SLURM\] submitted process (.+?) > jobId: (\d+);', line)
            if match:
                job = match[2]
                add(job, source=str(path))
                jobs[job]['task'] = match[1]
    return jobs


def query(command, warnings):
    try:
        result = subprocess.run(command, capture_output=True, text=True, timeout=30)
        if result.returncode:
            warnings.append(' '.join(command[:2])+': '+result.stderr.strip())
            return ''
        return result.stdout
    except (OSError, subprocess.TimeoutExpired) as error:
        warnings.append(str(error))
        return ''


def usage_rows(jobs, accounting, live):
    parents = {r['JobIDRaw']:r for r in accounting if '.' not in r['JobIDRaw']}
    steps = collections.defaultdict(list)
    for row in accounting:
        if '.' in row['JobIDRaw']:
            steps[row['JobIDRaw'].split('.')[0]].append(row)
    batches = {r['JobID']:r for r in live}
    result = []
    for job, details in sorted(jobs.items(), key=lambda x:x[0]):
        parent = parents.get(job, {})
        state = parent.get('State','UNAVAILABLE').split()[0].rstrip('+')
        tres = dict(x.split('=',1) for x in parent.get('AllocTRES','').split(',') if '=' in x)
        cpus = int(parent.get('AllocCPUS') or 0)
        elapsed = float(parent.get('ElapsedRaw') or 0)
        memory = gib(tres.get('mem',''))
        notes = []
        cpu, measurement = None, 'unavailable'
        if state in ACTIVE:
            batch = batches.get(job+'.batch', {})
            # Batch accounting is provisional; separate srun steps are not folded in.
            if batch.get('NTasks') and seconds(batch.get('AveCPU')) is not None:
                cpu = int(batch['NTasks'])*seconds(batch['AveCPU'])
                measurement = 'live_batch_estimate'
                notes.append('Provisional batch CPU; separate srun steps may be absent')
        elif state not in {'PENDING','UNAVAILABLE'}:
            cpu = seconds(parent.get('TotalCPU'))
            measurement = 'final_parent' if cpu is not None else 'unavailable'
            if state != 'COMPLETED':
                notes.append('Unsuccessful attempt; signal termination may undercount child CPU')
                if cpu == 0:
                    cpu, measurement = None, 'unavailable'
        observed, memory_source = None, ''
        batch = next((r for r in steps[job] if r['JobIDRaw']==job+'.batch'), {})
        observed = gib(batch.get('MaxRSS'))
        if observed is not None:
            memory_source = 'sacct.batch.MaxRSS'
        elif state in ACTIVE:
            observed = gib(batches.get(job+'.batch',{}).get('MaxRSS'))
            if observed is not None:
                memory_source = 'sstat.batch.MaxRSS (live)'
        if observed is None:
            values = [gib(r.get('MaxRSS')) for r in steps[job]]
            values = [v for v in values if v is not None]
            if values:
                observed, memory_source = max(values), 'maximum reported step MaxRSS'
        if observed is not None:
            notes.append('Memory is a step observation, not a guaranteed simultaneous whole-job peak')
        denominator = cpus*elapsed
        efficiency = 100*cpu/denominator if cpu is not None and denominator else None
        if efficiency is not None and efficiency > 100:
            notes.append('CPU percentage above100%; inspect accounting')
        gpu = float(tres['gres/gpu']) if 'gres/gpu' in tres else sum(
            float(v) for k,v in tres.items() if k.startswith('gres/gpu:'))
        if gpu:
            notes.append('GPU allocation only; GPU activity was not measured')
        task = details.get('task') or parent.get('JobName','')
        stage = task.split(' (')[0]
        if ':' in stage:
            stage = stage.rsplit(':',1)[-1]
        kind = 'workflow_task' if task.startswith('NFCORE_') or parent.get('JobName','').startswith('nf-') else 'support'
        queue = None
        try:
            queue = (datetime.fromisoformat(parent['Start'])-datetime.fromisoformat(parent['Eligible'])).total_seconds()
        except (KeyError,ValueError):
            pass
        values = dict(cluster=parent.get('Cluster',''),job_id=job,kind=kind,step=stage,
            task=task,attempt=details.get('attempt',''),state=state,exit_code=parent.get('ExitCode',''),
            allocated_cpus=cpus or None,allocated_memory_gib=memory,elapsed_seconds=elapsed,
            cpu_seconds=cpu,average_busy_cpus=cpu/elapsed if cpu is not None and elapsed else None,
            cpu_use_percent=efficiency,cpu_measurement=measurement,allocated_core_hours=denominator/3600,
            observed_maxrss_gib=observed,rss_fraction_percent=100*observed/memory if observed is not None and memory else None,
            memory_source=memory_source,queue_seconds=queue,allocated_gpus=gpu,gpu_card_hours=gpu*elapsed/3600,
            nextflow_cpu_percent=details.get('nextflow_cpu_percent',''),nextflow_peak_rss=details.get('nextflow_peak_rss',''),
            source=details.get('source',''),notes='; '.join(notes))
        result.append({k:('' if v is None else round(v,4) if isinstance(v,float) else v) for k,v in values.items()})
    return result


def render(root, rows, warnings, generated):
    esc = lambda v:html.escape(str(v), quote=True)
    completed = [r for r in rows if r['state']=='COMPLETED' and r['cpu_measurement']=='final_parent' and r['elapsed_seconds'] and r['allocated_cpus']]
    allocated = sum(r['allocated_cpus']*r['elapsed_seconds'] for r in completed)
    weighted = 100*sum(r['cpu_seconds'] for r in completed)/allocated if allocated else None
    fmt = lambda v:f'{v:,.2f}' if isinstance(v,(float,int)) else esc(v or 'Unavailable')
    table_columns = [('job_id','Job'),('step','Step'),('state','State'),('allocated_cpus','CPUs'),
        ('average_busy_cpus','Avg busy CPUs'),('cpu_use_percent','CPU use %'),('cpu_measurement','CPU source'),
        ('allocated_memory_gib','RAM GiB'),('observed_maxrss_gib','Observed RSS GiB'),('rss_fraction_percent','RSS / RAM %'),
        ('elapsed_seconds','Runtime'),('allocated_core_hours','Core-hours'),('notes','Notes')]
    table = '<tr>'+''.join('<th>'+label+'</th>' for _,label in table_columns)+'</tr>'
    for row in rows:
        runtime = int(row['elapsed_seconds'])
        duration = f'{runtime//3600}:{runtime%3600//60:02d}:{runtime%60:02d}'
        table += '<tr>'+''.join('<td>'+(duration if key=='elapsed_seconds' else fmt(row[key]))+'</td>' for key,_ in table_columns)+'</tr>'
    stage_rows = collections.defaultdict(list)
    for row in rows:
        stage_rows[row['kind'],row['step']].append(row)
    stages = '<tr><th>Kind / step</th><th>Jobs</th><th>Completed CPU use %</th><th>Metric coverage</th><th>Core-hours, all attempts</th></tr>'
    for (kind,stage), group in sorted(stage_rows.items()):
        measured = [r for r in group if r in completed]
        denom = sum(r['allocated_cpus']*r['elapsed_seconds'] for r in measured)
        pct = 100*sum(r['cpu_seconds'] for r in measured)/denom if denom else ''
        stages += '<tr>'+''.join('<td>'+fmt(v)+'</td>' for v in [kind+' / '+stage,len(group),pct,f'{len(measured)}/{len(group)}',sum(r['allocated_core_hours'] for r in group)])+'</tr>'
    warning_text = ''.join('<li>'+esc(v)+'</li>' for v in warnings)
    workload = ''
    plan_path = root/'state/input_resources.json'
    if plan_path.exists():
        try:
            plan = json.loads(plan_path.read_text())
            workload = '<div class="box"><h2>Recorded allocation workload</h2><p>'+esc(plan_path)+' · '+esc(plan.get('policy',''))+' · '+esc(plan.get('status',''))+'</p><table><tr><th>Sample</th><th>Role</th><th>Input GiB</th><th>Read pairs</th><th>Basis</th></tr>'
            for library in plan.get('libraries',[]):
                size = library.get('input_bytes',library.get('fastq_bytes'))
                workload += '<tr>'+''.join('<td>'+fmt(v)+'</td>' for v in [library.get('sample_id',''),library.get('role',''),size/1024**3 if size is not None else '',library.get('read_pairs'),library.get('workload_basis','')])+'</tr>'
            workload += '</table></div>'
        except (OSError,ValueError):
            workload = ''
    return f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>Resource use — {esc(root.name)}</title>
<style>body{{font:15px system-ui,sans-serif;margin:28px;color:#172b40;background:#f6f8fb}}h1{{font-size:24px}}p{{max-width:1100px;line-height:1.5}}.box{{background:white;padding:16px;margin:15px 0;border:1px solid #dce3ed;border-radius:8px}}.scroll{{overflow:auto}}table{{border-collapse:collapse;width:100%;font-size:13px}}td,th{{padding:9px;border-bottom:1px solid #e0e6ee;text-align:left;vertical-align:top}}th{{background:#eaf0f8;white-space:nowrap}}td{{max-width:340px;overflow-wrap:anywhere}}.stats{{display:flex;gap:24px;flex-wrap:wrap}}.stats b{{display:block;font-size:23px}}</style>
<h1>Job resource use: {esc(root.name)}</h1><p>Generated {esc(generated)} · Reporter v{VERSION} · Read-only snapshot. Running-job figures are provisional; refresh after completion.</p>
<div class="box stats"><div><b>{len(rows)}</b>Distinct jobs</div><div><b>{sum(r['state']=='COMPLETED' for r in rows)}</b>Completed</div><div><b>{sum(r['state'] in ACTIVE for r in rows)}</b>Active</div><div><b>{fmt(sum(r['allocated_core_hours'] for r in rows))}</b>Core-hours, all attempts</div><div><b>{fmt(weighted)}</b>Completed CPU use % ({len(completed)}/{len(rows)} jobs)</div></div>
<p>CPU use = actual CPU seconds / (allocated CPUs × elapsed seconds). Stage percentages are weighted by allocated CPU time. Controllers/support jobs are shown separately from workflow steps. Failed and cancelled attempts contribute to occupancy costs. Each job ID is counted once; parent totals are not added to their step totals. Native Nextflow CPU percentages count busy cores and are retained separately in the CSV.</p>
<p>Memory observations are labelled step-level RSS proxies, not guaranteed whole-job peaks. No automatic resource changes or pass/fail thresholds are applied. GPU card-hours describe allocation only; unavailable GPU activity is not zero activity.</p>
<p>Whole-run costs can include transfers, controllers and earlier failed or cancelled attempts. Compare equivalent steps and input scope before drawing conclusions about the BAM and FASTQ routes.</p>
{'<div class="box"><strong>Collection notes</strong><ul>'+warning_text+'</ul></div>' if warnings else ''}
{workload}
<div class="box scroll"><h2>By workflow step</h2><table>{stages}</table></div>
<div class="box scroll"><h2>Every job attempt</h2><table>{table}</table></div>
<p>Files beside this report: <a href="resource_usage.csv">resource_usage.csv</a>, <a href="sacct.psv">raw Slurm accounting</a>, <a href="sstat.psv">live accounting snapshot</a>. Sources: <a href="https://slurm.schedmd.com/sacct.html">Slurm accounting definitions</a>. No input BAM or FASTQ content was read.</p></html>'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir',type=Path,required=True)
    parser.add_argument('--out-dir',type=Path)
    parser.add_argument('--job-id',action='append',default=[],help='Additional ID or comma-separated IDs')
    parser.add_argument('--sacct-file',type=Path,help='Use saved accounting instead of querying Slurm')
    parser.add_argument('--sstat-file',type=Path,help='Use saved live accounting instead of querying Slurm')
    args = parser.parse_args()
    root = args.run_dir.resolve()
    output = args.out_dir or root/'results/pipeline_info/resource_usage'
    jobs = discover(root)
    for group in args.job_id:
        for job in group.split(','):
            if not re.fullmatch(r'\d+(?:_\d+)?',job):
                parser.error('Job IDs must be numeric Slurm IDs')
            jobs.setdefault(job,dict(source='explicit --job-id'))
    warnings = []
    if not jobs:
        warnings.append('No job IDs discovered yet; run again after task submission')
    fields = 'Cluster,JobIDRaw,JobName%180,State,ExitCode,AllocCPUS,AllocTRES%180,ElapsedRaw,TotalCPU,MaxRSS,TimelimitRaw,Submit,Eligible,Start,End,NodeList,NTasks'
    accounting = args.sacct_file.read_text() if args.sacct_file else query(
        ['sacct','--jobs',','.join(jobs),'--starttime=1970-01-01','--parsable2','--format='+fields],warnings) if jobs else ''
    parsed = records(accounting)
    active_ids = [r['JobIDRaw'] for r in parsed if '.' not in r['JobIDRaw'] and r['State'].split()[0] in ACTIVE]
    # Saved input mode is offline unless an explicit live file is supplied.
    live = args.sstat_file.read_text() if args.sstat_file else query(
        ['sstat','--allsteps','--jobs',','.join(active_ids),'--parsable2','--format=JobID,NTasks,AveCPU,MaxRSS'],warnings) if active_ids and not args.sacct_file else ''
    rows = usage_rows(jobs,parsed,records(live))
    missing = [r['job_id'] for r in rows if r['state']=='UNAVAILABLE']
    if missing:
        warnings.append('No accounting returned for: '+', '.join(missing))
    output.mkdir(parents=True,exist_ok=True)
    (output/'sacct.psv').write_text(accounting)
    (output/'sstat.psv').write_text(live)
    with (output/'resource_usage.csv').open('w',newline='') as handle:
        writer = csv.DictWriter(handle,fieldnames=FIELDS,lineterminator='\n')
        writer.writeheader();writer.writerows(rows)
    generated = datetime.now(ZoneInfo('America/New_York')).isoformat(timespec='seconds')
    (output/'resource_usage.html').write_text(render(root,rows,warnings,generated))
    print(json.dumps(dict(jobs=len(rows),out_dir=str(output),notes=warnings)))


if __name__=='__main__':
    main()
