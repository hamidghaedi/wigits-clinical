#!/usr/bin/env python3
"""Native BAM intake v0.2.0: full-stream verification and established preparation.

SAM decoding/checksums: samtools1.22.1. Mate preparation: GATK4.6.1.0 Picard
FixMateInformation. Source files never change. Full checksums include every
record, all sequence/quality data and non-mate tags; mapped coordinates are
checked separately because Picard may position unmapped reads beside their mate.
CRC checksums are standard samtools order-independent conservation checks, not
cryptographic file identity. Only DNA Illumina BAM input is released here.
"""
import argparse
import csv
import fcntl
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import subprocess
import time

VERSION = '0.2.0'


def identity(path):
    s = path.stat()
    return dict(path=str(path.resolve()), bytes=s.st_size, mtime_ns=s.st_mtime_ns)


def write_json(path, value):
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(value, indent=2) + '\n')
    os.replace(tmp, path)


def header(text):
    sq, rg, pg = {}, {}, []
    for line in text.splitlines():
        f = line.split('\t')
        tags = dict(v.split(':', 1) for v in f[1:] if ':' in v)
        if f[0] == '@SQ':
            sq[tags['SN']] = tags
        elif f[0] == '@RG':
            if tags.get('ID') in rg:
                raise ValueError('Duplicate read-group ID')
            rg[tags.get('ID')] = tags.get('SM')
        elif f[0] == '@PG':
            pg.append(tags)
    return sq, rg, pg


def check_header(text, dictionary, sample, reference):
    sq, rg, pg = header(text)
    expected, _, _ = header(dictionary.read_text())
    if {k:v.get('LN') for k,v in sq.items()} != {k:v.get('LN') for k,v in expected.items()}:
        raise ValueError('BAM/reference sequence dictionary mismatch')
    if not rg or set(rg.values()) != {sample}:
        raise ValueError('Read-group sample binding is missing or mixed')
    for line in text.splitlines():
        if line.startswith('@RG\t'):
            tags=dict(x.split(':',1) for x in line.split('\t')[1:] if ':' in x)
            if tags.get('PL','').upper()!='ILLUMINA':
                raise ValueError('Released intake requires explicit ILLUMINA read-group platform')
    md5 = [k for k in sq if sq[k].get('M5')]
    if any(sq[k]['M5'] != expected[k].get('M5') for k in md5):
        raise ValueError('BAM/reference M5 mismatch')
    if len(md5) == len(sq):
        return 'all header sequence M5 values match configured reference dictionary'
    # Older Parabricks BAMs lack M5. Preserve the weaker but explicit origin
    # declaration; do not claim that matching sequence lengths proves identity.
    declared = set()
    for record in pg:
        for token in shlex.split(record.get('CL', '')):
            if token.endswith(('.fa', '.fasta')):
                declared.add(Path(token).name)
    if declared != {reference.name}:
        raise ValueError('Reference identity unresolved: incomplete M5 and no unambiguous source FASTA declaration')
    return 'source PG FASTA declaration plus sequence dictionary; BAM lacks complete M5'


def checksum_rows(path):
    return sorted(' '.join(x.split()) for x in path.read_text().splitlines() if x.strip() and not x.startswith('#'))


def scan(bam, sample, prefix, audit_binary, log, commands, limit=None):
    full_path, mapped_path = Path(str(prefix)+'.content.chk'), Path(str(prefix)+'.mapped.chk')
    audit_path = Path(str(prefix)+'.audit.json')
    common = ['samtools', 'checksum', '-F', '0', '-b', '0xfd5', '-c', '-C', '-t', '*,MC,MQ', '-T']
    processes = []
    with full_path.open('w') as full, mapped_path.open('w') as mapped:
        try:
            c1 = subprocess.Popen(common, stdin=subprocess.PIPE, stdout=full, stderr=log)
            c2 = subprocess.Popen(common + ['-P'], stdin=subprocess.PIPE, stdout=mapped, stderr=log)
            decode = (['samtools','head','-n',str(limit),str(bam)] if limit else
                      ['samtools', 'view', '--no-PG', '-h', '-@', '4', str(bam)])
            commands.extend([common, common+['-P'], decode])
            d = subprocess.Popen(decode, stdout=subprocess.PIPE, stderr=log)
            processes = [c1, c2, d]
            cmd = [str(audit_binary), sample, str(audit_path), str(c2.stdin.fileno()), '1']
            commands.append(cmd)
            a = subprocess.Popen(cmd, stdin=d.stdout, stdout=c1.stdin,
                                 pass_fds=(c2.stdin.fileno(),), stderr=log)
            processes.append(a)
            d.stdout.close(); c1.stdin.close(); c2.stdin.close()
            codes = [p.wait() for p in [a, d, c1, c2]]
            if any(codes):
                raise RuntimeError(f'BAM full traversal/checksum failed: {codes}')
        finally:
            for p in processes:
                if p.poll() is None:
                    p.terminate()
    return dict(scope=f'first {limit} records only; not a full-file PASS' if limit else 'all records',
                audit=json.loads(audit_path.read_text()), audit_path=str(audit_path),
                content_checksum=str(full_path), mapped_checksum=str(mapped_path))


def check_audit(audit, prepared=False):
    hard = ['missing_or_invalid_rg', 'hard_clipped_supplementary', 'invalid_mc',
            'missing_primary_sequence', 'missing_primary_quality']
    if prepared:
        hard += ['missing_mc', 'coordinate_errors']
    failures = {k:audit[k] for k in hard if audit[k]}
    if not audit['primary']:
        failures['no_primary_records'] = True
    if failures:
        raise ValueError(f'Incompatible BAM records: {failures}')


def prepare(args):
    args.outdir.mkdir(parents=True, exist_ok=False)
    state = dict(status='INSPECTING', version=VERSION, source=identity(args.bam),
                 sample=args.sample, started_utc=datetime.now(timezone.utc).isoformat(), commands=[])
    receipt = args.outdir/'receipt.json'
    write_json(receipt, state)
    with (args.outdir/'commands.log').open('w') as log:
        def run(cmd, **kwargs):
            state['commands'].append(cmd)
            print(json.dumps(cmd), file=log, flush=True)
            return subprocess.run(cmd, check=True, stderr=log, **kwargs)
        try:
            state['samtools_version'] = run(['samtools','--version'], stdout=subprocess.PIPE).stdout.decode().splitlines()[0]
            if state['samtools_version'] != 'samtools 1.22.1':
                raise ValueError('Pinned samtools1.22.1 required')
            state['gatk_version']=run(shlex.split(args.gatk_command)+['--version'],stdout=subprocess.PIPE).stdout.decode()
            if '4.6.1.0' not in state['gatk_version']:
                raise ValueError('Pinned GATK4.6.1.0 required')
            digest = hashlib.sha256()
            with args.reference.open('rb') as h:
                for block in iter(lambda:h.read(8*1024**2), b''):digest.update(block)
            if digest.hexdigest() != args.reference_sha256:
                raise ValueError('Configured reference SHA-256 mismatch')
            state['reference'] = dict(**identity(args.reference), sha256=digest.hexdigest())
            derived_dictionary=args.outdir/'reference.dict'
            run(['samtools','dict','-o',str(derived_dictionary),str(args.reference)])
            actual_ref,_,_=header(derived_dictionary.read_text())
            supplied_ref,_,_=header(args.dictionary.read_text())
            if {k:(v.get('LN'),v.get('M5')) for k,v in actual_ref.items()}!={k:(v.get('LN'),v.get('M5')) for k,v in supplied_ref.items()}:
                raise ValueError('Supplied dictionary does not describe the verified FASTA')
            run(['samtools','quickcheck','-v',str(args.bam)])
            h = run(['samtools','view','-H',str(args.bam)],stdout=subprocess.PIPE).stdout.decode()
            state['reference_binding'] = check_header(h,derived_dictionary,args.sample,args.reference)
            # Header/sample incompatibility is rejected before reading the large BAM.
            audit_binary = args.outdir/'bam_intake_audit'
            run(['g++','-O3','-std=c++17',str(Path(__file__).with_name('bam_intake_audit.cpp')),'-o',str(audit_binary)])
            state['probe']=scan(args.bam,args.sample,args.outdir/'probe',audit_binary,log,state['commands'],limit=20000)
            check_audit(state['probe']['audit'])
            state['source_index']='present; fresh prepared-copy index required' if Path(str(args.bam)+'.bai').exists() else 'missing; prepared-copy index required'
            print(f'{args.sample}: full source inspection and conservation checksums',flush=True)
            state['before'] = scan(args.bam,args.sample,args.outdir/'source',audit_binary,log,state['commands'])
            check_audit(state['before']['audit'])
            state['status']='PREPARING_MATE_INFORMATION_AND_INDEX';write_json(receipt,state)
            work=args.outdir/'sorting';work.mkdir()
            out=args.outdir/f'{args.sample}.oncoanalyser.bam'
            threads=max(1,args.threads-8)
            name_sort=['samtools','sort','-N','-@',str(threads),'-m','1G','-l','0',
                       '-T',str(work/'name'),'-O','BAM',str(args.bam)]
            fixed=work/'mate_fixed.bam'
            fix=shlex.split(args.gatk_command)+['--java-options','-Xmx16g -XX:ParallelGCThreads=4',
                 'FixMateInformation','--INPUT','/dev/stdin','--OUTPUT',str(fixed),
                 '--ASSUME_SORTED','true','--SORT_ORDER','unsorted','--ADD_MATE_CIGAR','true',
                 '--IGNORE_MISSING_MATES','false','--COMPRESSION_LEVEL','1','--TMP_DIR',str(work)]
            coordinate_sort=['samtools','sort','-@',str(args.threads-1),'-m','1G','-T',str(work/'coordinate'),'-o',str(out),str(fixed)]
            processes=[]
            print(f'{args.sample}: preparing a new BAM with Picard FixMateInformation',flush=True)
            try:
                state['commands'].extend([name_sort,fix])
                p1=subprocess.Popen(name_sort,stdout=subprocess.PIPE,stderr=log);processes.append(p1)
                # GATK emits a textual completion result on stdout. Use its normal
                # BAM file output rather than mixing that text with binary data.
                p2=subprocess.Popen(fix,stdin=p1.stdout,stdout=log,stderr=log);processes.append(p2);p1.stdout.close()
                codes=[p.wait() for p in reversed(processes)]
                if any(codes):raise RuntimeError(f'Mate-information preparation failed: {codes}')
            finally:
                for p in processes:
                    if p.poll() is None:p.terminate()
            run(coordinate_sort)
            run(['samtools','index','-@',str(args.threads-1),str(out)])
            run(['samtools','quickcheck','-v',str(out)])
            run(['samtools','idxstats',str(out)],stdout=subprocess.DEVNULL)
            after_header=run(['samtools','view','-H',str(out)],stdout=subprocess.PIPE).stdout.decode()
            check_header(after_header,derived_dictionary,args.sample,args.reference)
            state['status']='VERIFYING_PREPARED_COPY';write_json(receipt,state)
            print(f'{args.sample}: verifying prepared records and preservation',flush=True)
            state['after']=scan(out,args.sample,args.outdir/'prepared',audit_binary,log,state['commands'])
            check_audit(state['after']['audit'],prepared=True)
            for key in ['records','primary','secondary','supplementary']:
                if state['before']['audit'][key]!=state['after']['audit'][key]:
                    raise ValueError('Record conservation failed: '+key)
            for key in ['content_checksum','mapped_checksum']:
                if checksum_rows(Path(state['before'][key]))!=checksum_rows(Path(state['after'][key])):
                    raise ValueError('Sequence/quality/non-mate conservation check failed: '+key)
            if identity(args.bam)!=state['source']:
                raise ValueError('Source changed during preparation')
            state.update(status='PASS',output=identity(out),index=identity(Path(str(out)+'.bai')),
                         completed_utc=datetime.now(timezone.utc).isoformat())
            fixed.unlink()  # Exact generated intermediate; source BAM is untouched.
            write_json(receipt,state)
            print(f'{args.sample}: PASS {out}',flush=True)
        except Exception as error:
            state.update(status='FAILED',error=str(error));write_json(receipt,state);raise


def publish(sheet, directory, output):
    with sheet.open() as handle:reader=csv.DictReader(handle);fieldnames=reader.fieldnames;rows=list(reader)
    for row in rows:
        if row['filetype']!='bam' or row['sequence_type']!='dna':
            raise ValueError('Released intake supports Illumina DNA BAM rows only')
        record=json.loads((directory/row['sample_id']/'receipt.json').read_text())
        if record['status']!='PASS' or record['sample']!=row['sample_id'] or record['source']!=identity(Path(row['filepath'])):
            raise ValueError('Intake not passed or source changed')
        for key in ['output','index']:
            if record[key]!=identity(Path(record[key]['path'])):
                raise ValueError('Prepared input changed after verification')
        row['filepath']=record['output']['path']
    tmp=output.with_suffix('.tmp')
    with tmp.open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fieldnames);writer.writeheader();writer.writerows(rows)
    os.replace(tmp,output)


def verify_launch(sheet, directory, reference=None, reference_sha256=None):
    with sheet.open() as handle:rows=list(csv.DictReader(handle))
    alignment_types={'bam','cram','bam_redux','cram_redux','redux_dir'}
    alignment_rows=[r for r in rows if r['filetype'] in alignment_types]
    if not alignment_rows:return
    if directory is None:raise ValueError('BAM intake required before analytical launch')
    if reference is None or reference_sha256 is None:
        raise ValueError('Verified BAM reference binding required')
    for row in alignment_rows:
        if row['filetype']!='bam' or row['sequence_type']!='dna':raise ValueError('Unreleased alignment input type')
        data=json.loads((directory/row['sample_id']/'receipt.json').read_text())
        if data['status']!='PASS' or data['version']!=VERSION or data['sample']!=row['sample_id']:
            raise ValueError('Matching passed intake required')
        ref=data['reference']
        if ref['sha256']!=reference_sha256 or {k:ref[k] for k in ['path','bytes','mtime_ns']}!=identity(reference):
            raise ValueError('Reference changed or differs from intake reference')
        if data['output']!=identity(Path(row['filepath'])):
            raise ValueError('Launch path is not the verified prepared BAM')
        if data['source']!=identity(Path(data['source']['path'])) or data['index']!=identity(Path(data['index']['path'])):
            raise ValueError('Source or prepared index changed')
        check_audit(data['after']['audit'],prepared=True)
        if any(data[x]['scope']!='all records' for x in ['before','after']):
            raise ValueError('A bounded probe cannot release BAM input')
        if any(data['before']['audit'][k]!=data['after']['audit'][k] for k in ['records','primary','secondary','supplementary']):
            raise ValueError('Record conservation evidence mismatch')
        for key in ['content_checksum','mapped_checksum']:
            if checksum_rows(Path(data['before'][key]))!=checksum_rows(Path(data['after'][key])):
                raise ValueError('Intake conservation evidence mismatch')


def batch(args):
    args.outdir.mkdir(parents=True,exist_ok=True)
    plan_path=args.outdir/'dispatch.json'
    state=json.loads(plan_path.read_text()) if plan_path.exists() else dict(jobs={},status='INTAKE_RUNNING')
    if state['status']=='PASS':
        publish(args.samplesheet,args.outdir,args.prepared_samplesheet)
        verify_launch(args.prepared_samplesheet,args.outdir,args.reference,args.reference_sha256)
        return
    with args.samplesheet.open() as handle:rows=list(csv.DictReader(handle))
    if len(rows)!=2 or {r['sample_type'] for r in rows}!={'tumor','normal'} or len({r['sample_id'] for r in rows})!=2:
        raise ValueError('Release requires one explicitly bound tumour/normal DNA pair')
    for row in rows:
        if row['filetype']!='bam' or row['sequence_type']!='dna':raise ValueError('Illumina DNA BAM input required')
        name=row['sample_id']
        if name not in state['jobs']:
            size=Path(row['filepath']).stat().st_size/(1024**3)
            hours=math.ceil(max(24,size/4*1.5)/4)*4
            if hours>240:raise ValueError('Intake workload exceeds released resource envelope')
            cmd=['sbatch','--parsable','--account='+args.account,'--cpus-per-task=64','--mem=192G',
                 f'--time={hours}:00:00','--output='+str(args.outdir/f'{name}.%j.log'),
                 str(Path(__file__).with_name('bam_intake_worker.sbatch')),str(Path(__file__).resolve()),
                 row['filepath'],name,str(args.outdir/name),str(args.reference),str(args.dictionary),args.reference_sha256,args.gatk_command]
            job=subprocess.check_output(cmd,text=True).strip().split(';')[0]
            state['jobs'][name]=dict(job_id=job,command=cmd,cpus=64,memory_gb=192,time_hours=hours)
            write_json(plan_path,state)
    while True:
        statuses={}
        for sample,j in state['jobs'].items():
            text=subprocess.check_output(['sacct','-n','-X','-j',j['job_id'],'--format=State,ExitCode','-P'],text=True).strip()
            statuses[sample]=text
            if any(x in text for x in ['FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL']):
                state.update(status='FAILED',accounting=statuses);write_json(plan_path,state)
                raise RuntimeError('Intake failed; no analytical launch: '+sample+' '+text)
        if all(x.startswith('COMPLETED|0:0') for x in statuses.values()):break
        time.sleep(30)
    publish(args.samplesheet,args.outdir,args.prepared_samplesheet)
    verify_launch(args.prepared_samplesheet,args.outdir,args.reference,args.reference_sha256)
    state.update(status='PASS',prepared_samplesheet=str(args.prepared_samplesheet));write_json(plan_path,state)


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='mode',required=True)
    for mode in ['prepare','batch']:
        q=sub.add_parser(mode)
        q.add_argument('--outdir',type=Path,required=True)
        q.add_argument('--reference',type=Path,required=True)
        q.add_argument('--dictionary',type=Path,required=True)
        q.add_argument('--reference-sha256',required=True)
        q.add_argument('--gatk-command',required=True)
        if mode=='prepare':
            q.add_argument('--bam',type=Path,required=True);q.add_argument('--sample',required=True);q.add_argument('--threads',type=int,required=True)
        else:
            q.add_argument('--samplesheet',type=Path,required=True);q.add_argument('--prepared-samplesheet',type=Path,required=True);q.add_argument('--account',required=True)
    q=sub.add_parser('verify-launch');q.add_argument('--samplesheet',type=Path,required=True);q.add_argument('--receipts',type=Path)
    q.add_argument('--reference',type=Path);q.add_argument('--reference-sha256')
    args=p.parse_args()
    if args.mode=='prepare':prepare(args)
    elif args.mode=='batch':
        args.outdir.mkdir(parents=True,exist_ok=True)
        with (args.outdir/'dispatcher.lock').open('w') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            batch(args)
    else:verify_launch(args.samplesheet,args.receipts,args.reference,args.reference_sha256)


if __name__=='__main__':main()
