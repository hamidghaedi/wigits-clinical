#!/usr/bin/env python3
"""Prepare native FASTQ inputs and complete-run laboratory QC scope.

Mappings are explicit: laboratory case/sample IDs bind libraries; analyte_type
and tumor_descriptor bind roles; fastq_parent_dir supplies the acquisition run;
lanes supplies expected lanes. Canonical staged files supply all actual chunks
and both mates. RunInfo or one FASTQ header supplies the native flowcell ID.
No read counting, payload hashing, transfer, analysis or per-case code generation.
The case samplesheet is selected; laboratory QC always retains the full run.
"""
import argparse
import csv
from datetime import datetime
import fcntl
import gzip
import glob
import hashlib
import json
from pathlib import Path
import re
import shlex
import sys
import xml.etree.ElementTree as ET

VERSION = '1.1.1'
FIELDS = 'group_id,subject_id,sample_id,sample_type,sequence_type,filetype,info,filepath'.split(',')
PHASES = {'Phase_1', 'phase_1.5', 'Phase_2'}
PHASE_ALIASES = {'Phase 1':'Phase_1','Phase 1.5':'phase_1.5','Phase 2':'Phase_2'}
REQUIRED = {'case_submitter_id', 'sample_submitter_id', 'analyte_type',
            'tumor_descriptor', 'specimen_type', 'tissue_preservation_method',
            'lanes', 'fastq_parent_dir'}
CASE_TOKEN = re.compile(r'(?i)(?<![A-Z0-9])KQ[-_ ]?0*([0-9]{1,4})(?![0-9])')


def case_key(value):
    """Same case-token rule as the tracker lab_manifest/seqera adapters."""
    matches={f'KQ{int(n):04d}' for n in CASE_TOKEN.findall(value or '')}
    if len(matches)>1: raise ValueError('Ambiguous patient identifier: '+value)
    return next(iter(matches)) if matches else value


def sample_key(value):
    """Format-only matching, preserving GDNA/RNA and numbered DNA aliquots.

    Phase1 lab standardizer passes2/3 handle separators. Tracker's case-token
    rule handles prefix zero padding. Never adopt the old broad pass4 which
    collapses GDNA/DNA or 1DNA/DNA. Non-cosmetic differences need sourced aliases.
    """
    value=re.sub(r'(?i)^KQ[-_ ]?0*([0-9]{1,4})(?=[-_ ]|$)',
                 lambda m:f'KQ{int(m[1]):04d}',value.strip())
    return re.sub(r'[-_\s]','',value).upper()


def source_paths(config, explicit=()):
    if explicit: return list(explicit),{}
    phases={};paths=set()
    for source in config.get('laboratory_sources',[]):
        for name in sorted(glob.glob(source['glob'],recursive=True)):
            path=Path(name).resolve()
            if path.suffix.lower()!='.csv':continue
            paths.add(path)
            phases[str(path)]=source.get('phase_by_file',{}).get(path.name,source.get('phase'))
    if not paths:raise ValueError('No original laboratory CSVs found in the configured laboratory_sources')
    return sorted(paths),phases


def load_aliases(config):
    result={};sources={}
    for name in config.get('sample_alias_files',[]):
        path=Path(name);sources[str(path.resolve())]=hashlib.sha256(path.read_bytes()).hexdigest()
        with path.open(encoding='utf-8-sig',newline='') as stream:
            for row in csv.DictReader(stream):
                old,new,evidence=(row.get(k,'').strip() for k in ('manifest_sample_id','storage_sample_id','evidence'))
                if not all((old,new,evidence)):raise ValueError('Alias requires original ID, storage ID and source evidence')
                if old in result and result[old]['storage_sample_id']!=new:raise ValueError('Conflicting sourced sample alias: '+old)
                if case_key(old)!=case_key(new):raise ValueError('Sample alias cannot change patient identity')
                result[old]=dict(storage_sample_id=new,evidence=evidence,source=str(path))
    return result,sources


def identifier(value):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', value):
        raise ValueError('Expected a literal identifier: '+value)
    return value


def save(path, value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value,indent=2)+'\n');tmp.replace(path)


def write_csv(path, rows):
    fields=list(dict.fromkeys(k for r in rows for k in r if not k.startswith('_')))
    with Path(path).open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,extrasaction='ignore',lineterminator='\n')
        writer.writeheader();writer.writerows(rows)


def laboratory(paths, phase, phase_map=None):
    result={};sources={}
    for path in paths:
        path=Path(path).resolve()
        sources[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
        with path.open(newline='',encoding='utf-8-sig') as stream:
            reader=csv.DictReader(stream)
            if not REQUIRED <= set(reader.fieldnames or []):
                raise ValueError('Missing laboratory columns in '+str(path)+': '+', '.join(sorted(REQUIRED-set(reader.fieldnames or []))))
            for number,original in enumerate(reader,2):
                row={k:(v or '').strip() for k,v in original.items() if k}
                if not row['sample_submitter_id']: continue
                parent=re.sub(r'/fastq\d*$','',row['fastq_parent_dir'].replace('\\','/').rstrip('/'))
                run=Path(parent).name if parent else ''
                row['_run_id']=run if run and run not in {'data','nfs','.'} else ''
                identifier(row['sample_submitter_id'])
                row['phase']=PHASE_ALIASES.get(row.get('phase'),row.get('phase')) or (phase_map or {}).get(str(path)) or phase
                row['_source']=str(path);row['_source_row']=number
                row['_original']=dict(original)
                key=row['_run_id'] or str(path),row['sample_submitter_id']
                if key in result and any(result[key].get(k)!=row.get(k) for k in REQUIRED|{'phase'}):
                    raise ValueError('Conflicting laboratory rows for '+row['sample_submitter_id'])
                result[key]=row
    if not result: raise ValueError('No laboratory libraries found')
    return list(result.values()),sources


def role(row):
    descriptor=row['tumor_descriptor'].casefold()
    if descriptor=='normal': return 'normal'
    if descriptor in {'primary','metastatic','recurrence','recurrent'}: return 'tumor'
    raise ValueError('Unresolved laboratory role: '+row['sample_submitter_id'])


def select(rows, case, tumour=None, normal=None, rna=None, dna_only=False, fastq_root=None):
    rows=[r for r in rows if case_key(r['case_submitter_id'])==case_key(case)]
    for row in rows:
        if CASE_TOKEN.search(row['case_submitter_id']) and case_key(row['sample_submitter_id'])!=case_key(case):
            raise ValueError('Laboratory case and sample IDs identify different patients: '+row['sample_submitter_id'])
    bound={}
    for key, analyte, kind, requested in [('tumor','DNA','tumor',tumour),('normal','DNA','normal',normal),('rna','RNA','tumor',rna)]:
        if key=='rna' and dna_only: continue
        choices=[r for r in rows if r['analyte_type'].upper()==analyte and role(r)==kind
                 and (not requested or r['sample_submitter_id']==requested)]
        if key=='rna' and not choices and not requested: continue
        if len(choices)>1 and fastq_root and not requested:
            staged=[r for r in choices if library_ready(r,fastq_root)]
            if len(staged)==1:choices=staged
            elif key=='rna':continue
        if len(choices)!=1:
            raise ValueError(f'Expected one {key} library for {case}; found {len(choices)}. Supply the matching laboratory manifests / explicit sample selector.')
        bound[key]=choices[0]
    return bound


def locate_fastqs(lab, fastq_root, allow_missing=False):
    """Match within the lab-declared run, with original names/paths retained."""
    run=lab['_run_id'];sample=lab['sample_submitter_id']
    if not run:
        if allow_missing:return []
        raise ValueError('Original laboratory run path is unresolved: '+sample)
    root=Path(fastq_root)/run
    key=sample_key(lab.get('_storage_sample_id',sample))
    # Inspect only this acquisition run; never search the entire filesystem.
    folders=[p for p in root.iterdir() if p.is_dir() and sample_key(p.name)==key] if root.is_dir() else []
    candidates=[]
    for folder in folders:
        for path in sorted(folder.glob('*.fastq.gz')):
            m=re.fullmatch(r'(.+)_S(\d+)_L(\d{3})_R([12])_(\d{3})\.fastq\.gz',path.name)
            if not m or sample_key(m[1])!=key:raise ValueError('Unbound or conflicting FASTQ in sample folder: '+str(path))
            candidates.append((path,m.groups()[1:]))
    if len(folders)>1:raise ValueError('Multiple storage folders match the same laboratory library: '+sample)
    if not candidates and not allow_missing:raise ValueError('No staged FASTQs for declared laboratory library: '+sample)
    return candidates


def library_ready(lab, fastq_root):
    candidates=locate_fastqs(lab,fastq_root,allow_missing=True)
    expected={int(x) for x in re.split(r'[;, ]+',lab['lanes']) if x}
    pairs={}
    for path,(number,lane,mate,chunk) in candidates:
        if not path.stat().st_size:return False
        pairs.setdefault((number,lane,chunk),set()).add(int(mate))
    return bool(pairs) and {int(k[1]) for k in pairs}==expected and all(p=={1,2} for p in pairs.values())


def header_flowcell(path, lane, mate):
    # One header only; no counting, whole-file integrity pass or BAM scan.
    with gzip.open(path,'rt') as stream: header=stream.readline(4096).strip().split()
    fields=header[0].split(':') if header else []
    if len(fields)!=7 or not fields[0].startswith('@') or fields[3]!=str(int(lane)):
        raise ValueError('No usable Illumina header; supply this run\'s RunInfo.xml: '+str(path))
    if len(header)>1 and header[1].split(':')[0]!=str(mate):
        raise ValueError('FASTQ mate header differs from filename: '+str(path))
    return fields[2]


def fastq_rows(case, bound, fastq_root, run_infos):
    rows=[];objects=[]
    for key,lab in bound.items():
        sample=lab['sample_submitter_id'];run=lab['_run_id']
        native_sample=sample
        expected={int(x) for x in re.split(r'[;, ]+',lab['lanes']) if x}
        if not expected: raise ValueError('Missing laboratory lanes: '+sample)
        pairs={}
        for path,parts in locate_fastqs(lab,fastq_root):
            sample_number,lane,mate,chunk=parts
            if path.stat().st_size==0: raise ValueError('Empty staged FASTQ: '+str(path))
            pairs.setdefault((sample_number,lane,chunk),{})[int(mate)]=path
        observed={int(k[1]) for k in pairs}
        if observed!=expected: raise ValueError(f'{sample}: expected lanes {sorted(expected)}, found {sorted(observed)}')
        for (number,lane,chunk),mates in sorted(pairs.items()):
            if set(mates)!={1,2}: raise ValueError(f'{sample}: missing mate in lane {lane}, chunk {chunk}')
            flowcells=[run_infos[run]] if run in run_infos else [header_flowcell(mates[m],lane,m) for m in (1,2)]
            if len(set(flowcells))!=1: raise ValueError('FASTQ mate flowcells differ: '+sample)
            flowcell=identifier(flowcells[0])
            rows.append(dict(zip(FIELDS,[case,case,native_sample,'tumor' if key=='rna' else key,
                'rna' if key=='rna' else 'dna','fastq',
                f'library_id:{native_sample};lane:{lane};flowcell:{flowcell}',';'.join(str(mates[m].resolve()) for m in (1,2))])))
            objects.extend(dict(sample_id=sample,run_id=run,lane=lane,mate=str(m),chunk=chunk,
                storage_sample_id=mates[m].parent.name,
                match_method='source_confirmed_alias' if lab.get('_alias_evidence') else
                    'exact' if mates[m].name.startswith(sample+'_S') and mates[m].parent.name==sample else 'phase1_separator_or_patient_padding',
                destination=str(mates[m].resolve()),expected_bytes=mates[m].stat().st_size) for m in (1,2))
    return rows,objects


def register_labs(config, rows, sources):
    """Persist entire DNA-run scopes, independently of the selected case."""
    project=Path(config['project_root']);root=Path(config.get('qc_registry_root',project/'file_manifests/lab_qc'))
    root.mkdir(parents=True,exist_ok=True)
    with (root/'.registry.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        qc_path=root/'config.json'
        qc=json.loads(qc_path.read_text()) if qc_path.exists() else dict(config['lab_qc'])
        for run in sorted({r['_run_id'] for r in rows}):
            full=[r for r in rows if r['_run_id']==run and r['analyte_type'].upper()=='DNA']
            if not full: continue
            for row in full: row['analyte_type']='DNA'
            digest=hashlib.sha256(json.dumps(full,sort_keys=True).encode()).hexdigest()
            directory=root/'laboratory'/run/digest[:12];directory.mkdir(parents=True,exist_ok=True)
            path=directory/'manifest.csv'
            # A revised whole-run source supersedes its prior active QC membership.
            old_paths=[p for p in qc.get('laboratories',[]) if Path(p).parent.parent.name==run]
            incoming={r['sample_submitter_id']:r for r in full}
            for old_path in old_paths:
                with Path(old_path).open() as stream: old=list(csv.DictReader(stream))
                for old_row in old:
                    sample=old_row['sample_submitter_id']
                    if sample not in incoming or any(incoming[sample].get(k)!=old_row.get(k) for k in ('lanes','phase')):
                        raise ValueError('Cannot shrink or alter registered full-run QC scope: '+run)
            previous=[p for p in qc.get('laboratories',[]) if p not in old_paths]
            if not path.exists(): write_csv(path,full);save(directory/'provenance.json',dict(sources=sources))
            qc['laboratories']=sorted(set(previous+[str(path)]))
        qc['package_roots']=sorted(set(qc.get('package_roots',[])+[str(config.get('_manifest_root',config.get('manifest_root',project/'file_manifests')))]))
        qc['config_path']=str(qc_path)
        save(qc_path,qc)
        return qc_path


def prepare(args):
    config=json.loads(args.config.read_text())
    requested=args.case_id or args.sample_id
    if not requested:raise ValueError('Provide the patient/sample ID')
    case=identifier(case_key(requested))
    required={'project_root','scratch_root','workflow_version','workflow_directory','workflow_commit',
              'native_source','native_source_version','native_source_commit','reference_root','extension_resources','lab_qc',
              'resource_version','driver_panel_sha256'}
    if not required<=config.keys(): raise ValueError('Incomplete shared site configuration: '+', '.join(sorted(required-config.keys())))
    paths,phases=source_paths(config,args.lab_manifest or [])
    rows,sources=laboratory(paths,args.phase,phases)
    aliases,alias_sources=load_aliases(config);sources.update(alias_sources)
    for row in rows:
        if row['sample_submitter_id'] in aliases:
            row['_storage_sample_id']=aliases[row['sample_submitter_id']]['storage_sample_id']
            row['_alias_evidence']=aliases[row['sample_submitter_id']]
    declared=[r for r in rows if case_key(r['case_submitter_id'])==case]
    if requested!=case:
        exact=[r for r in declared if sample_key(r['sample_submitter_id'])==sample_key(requested)]
        if len(exact)==1:
            selector='rna_sample' if exact[0]['analyte_type'].upper()=='RNA' else 'normal_sample' if role(exact[0])=='normal' else 'tumour_sample'
            if not getattr(args,selector):setattr(args,selector,exact[0]['sample_submitter_id'])
    fastq_root=Path(config.get('fastq_root',Path(config['scratch_root'])/'fastq'))
    bound=select(rows,case,args.tumour_sample,args.normal_sample,args.rna_sample,args.dna_only,fastq_root)
    declared_bound=dict(bound)
    rna_declared=[r for r in declared if r['analyte_type'].upper()=='RNA']
    rna_status='NOT_DECLARED' if not rna_declared else 'NOT_SELECTED' if args.dna_only else 'AMBIGUOUS'
    if 'rna' in bound:
        rna_status='READY' if library_ready(bound['rna'],fastq_root) else 'NOT_STAGED'
        if rna_status!='READY':bound.pop('rna')
    if any(r['phase'] not in PHASES for r in bound.values()):
        raise ValueError('Selected original source needs its phase configured once in laboratory_sources')
    if bound['tumor']['tissue_preservation_method'].casefold() not in {'frozen','fresh frozen','fresh-frozen'}:
        raise ValueError('Routine native route requires recorded fresh-frozen tumour; FFPE needs a separately authorized route')
    if args.blood_ch and bound['normal']['specimen_type'].casefold()!='blood':
        raise ValueError('Blood CH requires a laboratory-confirmed blood normal')
    if args.rna_neoantigens and 'rna' not in bound:
        raise ValueError('RNA neoantigens require selected tumour RNA')
    run_infos={}
    info_root=config.get('run_info_root')
    paths=list(args.run_info)
    if info_root:
        paths.extend(Path(info_root)/(r['_run_id']+'.RunInfo.xml') for r in bound.values())
    for path in paths:
        if not path.is_file(): continue
        run=ET.parse(path).getroot().find('Run')
        if run is None or not run.get('Id') or not run.findtext('Flowcell'):
            raise ValueError('Unusable RunInfo.xml: '+str(path))
        run_infos[run.get('Id')]=run.findtext('Flowcell')
        sources[str(path.resolve())]=hashlib.sha256(path.read_bytes()).hexdigest()
    samplesheet,objects=fastq_rows(case,bound,fastq_root,run_infos)
    aid=identifier(args.analysis_id or f'{case}-native-{datetime.now().strftime("%Y%m%d")}-r1')
    config['_manifest_root']=str(config.get('manifest_root',Path(config['project_root'])/'file_manifests')).format(phase=bound['tumor']['phase'])
    out=Path(config['_manifest_root'])/aid
    if out.exists(): raise ValueError('Prepared inputs already exist; use them to run/resume, or choose another --analysis-id: '+str(out))
    # Enroll only acquisition runs bound to this case, retaining every library.
    runs={r['_run_id'] for r in bound.values()}
    qc=register_labs(config,[r for r in rows if r['_run_id'] in runs],sources)
    out.mkdir(parents=True)
    run_root=Path(config.get('run_root',Path(config['scratch_root'])/'analysis/WiGiTS'/config['workflow_version']))
    workflow=Path(config['workflow_directory']);reference=Path(config['reference_root'])
    with (out/'samplesheet.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=FIELDS,lineterminator='\n');writer.writeheader();writer.writerows(samplesheet)
    laboratory_path=out/'laboratory.json'
    save(laboratory_path,dict(sources=sources,selected_libraries=bound,full_run_qc_config=str(qc)))
    # This is the authoritative generated case input. CSV/options/site files are
    # projections, not independent places to edit patient bindings.
    case_manifest=dict(schema_version='wigits.case-input.v1',case_id=case,phase=bound['tumor']['phase'],
        input_selector=requested,analysis_id=aid,run_directory=str(run_root/aid),
        assay_mode='WGTS' if 'rna' in bound else 'WGS',route='native_hartwig_fresh_frozen',
        options=dict(blood_ch=args.blood_ch,rna_neoantigens=args.rna_neoantigens),
        native_source=config['native_source'],native_source_commit=config['native_source_commit'],
        resource_version=config['resource_version'],driver_panel_sha256=config['driver_panel_sha256'],
        selected_libraries=declared_bound,declared_libraries=declared,
        runnable_libraries=bound,rna=dict(status=rna_status,declared_sample_ids=[r['sample_submitter_id'] for r in rna_declared]),
        files=objects,source_files=sources,workflow_version=config['workflow_version'],
        workflow_commit=config['workflow_commit'],full_run_qc_config=str(qc),status='READY',
        matching_policy='exact_then_unique_separator_and_case_prefix_normalization; sourced_noncosmetic_aliases',
        original_laboratory_values_preserved=True)
    save(out/'case_manifest.json',case_manifest)
    # Reuse approved clinical/runtime sidecars from enrolled preparation sources.
    # They are supplementary context, never an alternative identity source.
    contexts=[]
    for pattern in config.get('context_sources',[]):
        for name in sorted(glob.glob(pattern)):
            path=Path(name);record=json.loads(path.read_text())
            if record.get('case_id')==case:contexts.append((path,record))
    runtime=[]
    if contexts:
        labels={r.get('native_cancer_type') for _,r in contexts if r.get('native_cancer_type')}
        if len(labels)>1:raise ValueError('Conflicting approved native clinical context')
        save(out/'clinical_context.json',dict(sources=[str(p) for p,_ in contexts],records=[r for _,r in contexts]))
        for path,_ in contexts:
            extra=path.with_name('native_runtime.config')
            if extra.is_file() and extra.read_text() not in [p.read_text() for p in runtime]:runtime.append(extra)
        if len(runtime)>1:raise ValueError('Conflicting approved native runtime sidecars')
        if runtime:(out/'native_runtime.config').write_bytes(runtime[0].read_bytes())
        case_manifest['approved_context_sources']=[str(p) for p,_ in contexts]
        case_manifest['approved_runtime_sources']=[str(p) for p in runtime]
        save(out/'case_manifest.json',case_manifest)
    resources=dict(config['extension_resources'])
    optional=dict(schema_version='wigits.optional-extensions.v1',resources=resources,cases=[dict(case_id=case,
        results_case_dir=str(run_root/aid/'results'/case),native_genome='GRCh38_hmf',
        tumour_sample=bound['tumor']['sample_submitter_id'],normal_sample=bound['normal']['sample_submitter_id'],
        rna_sample=bound.get('rna',{}).get('sample_submitter_id',''),
        laboratory=dict(source=str(laboratory_path),normal_specimen_type=bound['normal']['specimen_type'],tumour_rna_available='rna' in bound),
        options=dict(blood_ch=args.blood_ch,rna_neoantigens=args.rna_neoantigens))])
    save(out/'optional_extensions.json',optional)
    values=dict(WIGITS_PROJECT_ROOT=config['project_root'],WIGITS_SCRATCH_ROOT=config['scratch_root'],
        WIGITS_RUN_ROOT=str(run_root),WIGITS_RUN_ID=aid,WIGITS_WORKFLOW_DIR=str(workflow),
        WIGITS_EXTENSIONS_CONFIG=str(out/'optional_extensions.json'),WIGITS_LAB_QC_CONFIG=str(qc),
        ONCO_SOURCE=config['native_source'],REF_CONFIG=str(reference/'reference_data.config'))
    if runtime:values['WIGITS_NATIVE_EXTRA_CONFIG']=str(out/'native_runtime.config')
    for key,cfgkey in [('WIGITS_SLURM_ACCOUNT','account'),('WIGITS_NEXTFLOW_MODULE','nextflow_module'),('WIGITS_APPTAINER_MODULE','apptainer_module')]:
        if config.get(cfgkey): values[key]=config[cfgkey]
    site=config.get('site','fir')
    if site not in {'fir','frontenac'}:raise ValueError('Unknown configured site')
    (out/'site.env').write_text('\n'.join('export '+k+'='+shlex.quote(str(v)) for k,v in values.items())+
                              '\nsource "$WIGITS_WORKFLOW_DIR/conf/'+site+'_paths.env"\n')
    save(out/'input_plan.json',dict(case_id=case,run_id=aid,run_directory=str(run_root/aid),
        libraries=bound,objects=objects))
    save(out/'release.json',dict(preparer_version=VERSION,workflow_version=config['workflow_version'],
        workflow_commit=config['workflow_commit'],native_source_version=config['native_source_version'],
        native_source=config['native_source'],native_source_commit=config['native_source_commit'],sources=sources,
        resource_version=config['resource_version'],driver_panel_sha256=config['driver_panel_sha256'],
        site_config_source=str(args.config.resolve()),site_config_sha256=hashlib.sha256(args.config.read_bytes()).hexdigest(),
        laboratory_qc_config=str(qc),route='native_hartwig_fresh_frozen',input_files=len(objects),
        blood_ch=args.blood_ch,rna_neoantigens=args.rna_neoantigens))
    return dict(prepared=str(out),case_id=case,lane_rows=len(samplesheet),fastqs=len(objects),
        run_directory=str(run_root/aid),lab_qc_config=str(qc),analysis_submitted=False,
        library_ids={key:lab['sample_submitter_id'] for key,lab in declared_bound.items()},rna_status=rna_status)


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('sample_id',nargs='?',help='Patient ID or full laboratory sample ID; original lab files are discovered automatically')
    p.add_argument('--lab-manifest',type=Path,action='append',help='Optional explicit original CSV override; normally scan laboratory_sources')
    p.add_argument('--case-id',help='Backward-compatible explicit patient selector')
    p.add_argument('--phase',choices=sorted(PHASES),help='Optional phase override where absent in the source configuration')
    p.add_argument('--config',type=Path,default=Path(__file__).resolve().parents[1]/'config/native_fir.json',help='Shared site JSON, automatically found in the project config directory')
    p.add_argument('--analysis-id',help='Optional explicit identifier; default <case>-native-<date>-r1')
    for kind in ('tumour','normal','rna'):p.add_argument('--'+kind+'-sample')
    p.add_argument('--dna-only',action='store_true',help='Explicitly omit laboratory RNA libraries')
    p.add_argument('--blood-ch',action='store_true')
    p.add_argument('--rna-neoantigens',action='store_true')
    p.add_argument('--run-info',type=Path,action='append',default=[],help='Optional original RunInfo.xml; otherwise read one FASTQ header')
    p.add_argument('--json',action='store_true',help='Machine-readable output instead of the short human summary')
    return p


if __name__=='__main__':
    try:
        args=parser().parse_args();result=prepare(args)
        if args.json:print(json.dumps(result,indent=2))
        else:
            print('Prepared: '+result['prepared'])
            for key,sample in result['library_ids'].items():print(key+': '+sample)
            print('RNA: '+result['rna_status'])
            profile=json.loads(args.config.read_text());folder=Path(result['prepared'])
            launcher=profile.get('launcher',str(Path(profile['project_root'])/'scripts/run_native_fir.sh'))
            print('Run:\n'+shlex.join(['bash',launcher,
                '--manifest',str(folder/'samplesheet.csv'),'--config',str(folder/'site.env')]))
    except (ValueError,KeyError,FileNotFoundError) as error:
        print(str(error),file=sys.stderr);sys.exit(2)
