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
import hashlib
import json
from pathlib import Path
import re
import shlex
import sys
import xml.etree.ElementTree as ET

VERSION = '1.0.0'
FIELDS = 'group_id,subject_id,sample_id,sample_type,sequence_type,filetype,info,filepath'.split(',')
PHASES = {'Phase_1', 'phase_1.5', 'Phase_2'}
PHASE_ALIASES = {'Phase 1':'Phase_1','Phase 1.5':'phase_1.5','Phase 2':'Phase_2'}
REQUIRED = {'case_submitter_id', 'sample_submitter_id', 'analyte_type',
            'tumor_descriptor', 'specimen_type', 'tissue_preservation_method',
            'lanes', 'fastq_parent_dir'}


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


def laboratory(paths, phase):
    result={};sources={}
    for path in paths:
        path=Path(path).resolve()
        sources[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
        with path.open(newline='',encoding='utf-8-sig') as stream:
            reader=csv.DictReader(stream)
            if not REQUIRED <= set(reader.fieldnames or []):
                raise ValueError('Missing laboratory columns in '+str(path)+': '+', '.join(sorted(REQUIRED-set(reader.fieldnames or []))))
            for original in reader:
                row={k:(v or '').strip() for k,v in original.items() if k}
                if not row['sample_submitter_id']: continue
                row['_run_id']=identifier(Path(row['fastq_parent_dir']).parent.name)
                identifier(row['sample_submitter_id'])
                row['phase']=PHASE_ALIASES.get(row.get('phase'),row.get('phase')) or phase
                if row['phase'] not in PHASES: raise ValueError('Unknown laboratory phase')
                key=row['_run_id'],row['sample_submitter_id']
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


def select(rows, case, tumour=None, normal=None, rna=None, dna_only=False):
    rows=[r for r in rows if r['case_submitter_id']==case or r['case_submitter_id'].split('-',1)[0]==case]
    bound={}
    for key, analyte, kind, requested in [('tumor','DNA','tumor',tumour),('normal','DNA','normal',normal),('rna','RNA','tumor',rna)]:
        if key=='rna' and dna_only: continue
        choices=[r for r in rows if r['analyte_type'].upper()==analyte and role(r)==kind
                 and (not requested or r['sample_submitter_id']==requested)]
        if key=='rna' and not choices and not requested: continue
        if len(choices)!=1:
            raise ValueError(f'Expected one {key} library for {case}; found {len(choices)}. Supply the matching laboratory manifests / explicit sample selector.')
        bound[key]=choices[0]
    return bound


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
        sample=lab['sample_submitter_id'];run=lab['_run_id'];folder=Path(fastq_root)/run/sample
        expected={int(x) for x in re.split(r'[;, ]+',lab['lanes']) if x}
        if not expected: raise ValueError('Missing laboratory lanes: '+sample)
        pattern=re.compile(re.escape(sample)+r'_S(\d+)_L(\d{3})_R([12])_(\d{3})\.fastq\.gz$')
        pairs={}
        for path in sorted(folder.glob('*.fastq.gz')):
            match=pattern.fullmatch(path.name)
            if not match: raise ValueError('Staged FASTQ name does not match the laboratory library: '+str(path))
            sample_number,lane,mate,chunk=match.groups()
            if path.stat().st_size==0: raise ValueError('Empty staged FASTQ: '+str(path))
            pairs.setdefault((sample_number,lane,chunk),{})[int(mate)]=path
        observed={int(k[1]) for k in pairs}
        if observed!=expected: raise ValueError(f'{sample}: expected lanes {sorted(expected)}, found {sorted(observed)}')
        for (number,lane,chunk),mates in sorted(pairs.items()):
            if set(mates)!={1,2}: raise ValueError(f'{sample}: missing mate in lane {lane}, chunk {chunk}')
            flowcells=[run_infos[run]] if run in run_infos else [header_flowcell(mates[m],lane,m) for m in (1,2)]
            if len(set(flowcells))!=1: raise ValueError('FASTQ mate flowcells differ: '+sample)
            flowcell=identifier(flowcells[0])
            rows.append(dict(zip(FIELDS,[case,case,sample,'tumor' if key=='rna' else key,
                'rna' if key=='rna' else 'dna','fastq',
                f'library_id:{sample};lane:{lane};flowcell:{flowcell}',';'.join(str(mates[m].resolve()) for m in (1,2))])))
            objects.extend(dict(sample_id=sample,run_id=run,lane=lane,mate=str(m),chunk=chunk,
                destination=str(mates[m].resolve()),expected_bytes=mates[m].stat().st_size) for m in (1,2))
    return rows,objects


def register_labs(config, rows, sources):
    """Persist entire DNA-run scopes, independently of the selected case."""
    project=Path(config['project_root']);root=project/'file_manifests/lab_qc'
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
        qc['package_roots']=[str(project/'file_manifests')]
        qc['config_path']=str(qc_path)
        save(qc_path,qc)
        return qc_path


def prepare(args):
    config=json.loads(args.config.read_text());case=identifier(args.case_id)
    required={'project_root','scratch_root','workflow_version','workflow_directory','workflow_commit',
              'native_source','native_source_version','native_source_commit','reference_root','extension_resources','lab_qc',
              'resource_version','driver_panel_sha256'}
    if not required<=config.keys(): raise ValueError('Incomplete shared site configuration: '+', '.join(sorted(required-config.keys())))
    rows,sources=laboratory(args.lab_manifest,args.phase)
    bound=select(rows,case,args.tumour_sample,args.normal_sample,args.rna_sample,args.dna_only)
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
    samplesheet,objects=fastq_rows(case,bound,Path(config['scratch_root'])/'fastq',run_infos)
    aid=identifier(args.analysis_id or f'{case}-native-{datetime.now().strftime("%Y%m%d")}-r1')
    out=Path(config['project_root'])/'file_manifests'/aid
    if out.exists(): raise ValueError('Prepared inputs already exist; use them to run/resume, or choose another --analysis-id: '+str(out))
    # Enroll only acquisition runs bound to this case, retaining every library.
    runs={r['_run_id'] for r in bound.values()}
    qc=register_labs(config,[r for r in rows if r['_run_id'] in runs],sources)
    out.mkdir()
    run_root=Path(config['scratch_root'])/'analysis/WiGiTS'/config['workflow_version']
    workflow=Path(config['workflow_directory']);reference=Path(config['reference_root'])
    with (out/'samplesheet.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=FIELDS,lineterminator='\n');writer.writeheader();writer.writerows(samplesheet)
    laboratory_path=out/'laboratory.json'
    save(laboratory_path,dict(sources=sources,selected_libraries=bound,full_run_qc_config=str(qc)))
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
    for key,cfgkey in [('WIGITS_SLURM_ACCOUNT','account'),('WIGITS_NEXTFLOW_MODULE','nextflow_module'),('WIGITS_APPTAINER_MODULE','apptainer_module')]:
        if config.get(cfgkey): values[key]=config[cfgkey]
    (out/'site.env').write_text('\n'.join('export '+k+'='+shlex.quote(str(v)) for k,v in values.items())+
                              '\nsource "$WIGITS_WORKFLOW_DIR/conf/fir_paths.env"\n')
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
        run_directory=str(run_root/aid),lab_qc_config=str(qc),analysis_submitted=False)


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--lab-manifest',type=Path,action='append',required=True,help='Original laboratory CSV; repeat for tumour/normal/RNA runs')
    p.add_argument('--case-id',required=True)
    p.add_argument('--phase',choices=sorted(PHASES),required=True,help='Used only where the laboratory CSV has no phase')
    p.add_argument('--config',type=Path,required=True,help='Shared site JSON, reused across samples')
    p.add_argument('--analysis-id',help='Optional explicit identifier; default <case>-native-<date>-r1')
    for kind in ('tumour','normal','rna'):p.add_argument('--'+kind+'-sample')
    p.add_argument('--dna-only',action='store_true',help='Explicitly omit laboratory RNA libraries')
    p.add_argument('--blood-ch',action='store_true')
    p.add_argument('--rna-neoantigens',action='store_true')
    p.add_argument('--run-info',type=Path,action='append',default=[],help='Optional original RunInfo.xml; otherwise read one FASTQ header')
    return p


if __name__=='__main__':
    try: print(json.dumps(prepare(parser().parse_args()),indent=2))
    except (ValueError,KeyError,FileNotFoundError) as error:
        print(str(error),file=sys.stderr);sys.exit(2)
