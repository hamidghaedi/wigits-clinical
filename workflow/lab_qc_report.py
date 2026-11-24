#!/usr/bin/env python3
"""Run-level laboratory QC from native tool outputs; invoked by LAB_RUN_QC.

Each manifest row is one actual laboratory library, grouped by sequencing run
(not paired analysis). Sum fastp before_filtering reads/bases/Q30 bases across
all recorded lanes; Q30 = summed Q30 bases / summed bases. fastp is QC-only,
with all trimming/filtering disabled and no FASTQ outputs. No post-trimming
yield or raw instrument PF/undetermined read count is inferred. Picard WGS
metrics remain their native filtered measurements; BAMTOOLS duplicate flag
counts remain their actual producer's observations, not Picard duplicate metrics.
Copy the native flag text verbatim into a separately labelled samtools-format
parser input. Original source paths and small-file hashes remain in provenance.
There is no specimen clearance, caller/filter change or BAM rewrite here.
"""
import argparse
import base64
import csv
import hashlib
import html
import json
from pathlib import Path
import re
import shutil
import zipfile


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def metric_rows(path):
    lines = Path(path).read_text().splitlines()
    for i, line in enumerate(lines):
        if line.startswith('## METRICS CLASS'):
            end = next((j for j in range(i+2, len(lines)) if not lines[j] or lines[j].startswith('#')), len(lines))
            return list(csv.DictReader(lines[i+1:end], delimiter='\t'))
    raise ValueError('No Picard metric table: ' + str(path))


def raw_summary(paths):
    result = dict(total_reads=0, total_bases=0, q30_bases=0)
    for path in paths:
        data = json.loads(Path(path).read_text())['summary']
        before, after = data['before_filtering'], data['after_filtering']
        if before['total_reads'] != after['total_reads'] or before['total_bases'] != after['total_bases']:
            raise ValueError('QC-only fastp unexpectedly changed read/base yield: ' + str(path))
        for key in result:
            result[key] += before[key]
    result['q30_pct'] = 100 * result['q30_bases'] / result['total_bases'] if result['total_bases'] else None
    return result


def duplicate_summary(path):
    counts = {}
    for line in Path(path).read_text().splitlines():
        match = re.match(r'^(\d+) \+ (\d+) (primary(?: duplicates)?)(?:\s|$)', line)
        if match:
            counts[match[3]] = int(match[1]) + int(match[2])
    if not counts.get('primary') or 'primary duplicates' not in counts:
        raise ValueError('Missing primary duplicate/primary counts: ' + str(path))
    return dict(primary_reads=counts['primary'], primary_duplicates=counts['primary duplicates'],
                bam_duplicate_pct=100 * counts['primary duplicates'] / counts['primary'])


def load_samples(manifest, run):
    samples = [s for s in manifest['samples'] if s['run_id'] == run]
    if not samples or len({s['sample_id'] for s in samples}) != len(samples):
        raise ValueError('Expected nonempty unique run membership')
    for s in samples:
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9._-]*', s['sample_id']):
            raise ValueError('Unsafe sample ID')
        if not s['fastqs'] or len({x['label'] for x in s['fastqs']}) != len(s['fastqs']):
            raise ValueError('Missing or duplicate lane pairs')
    return samples


def prepare(args):
    manifest = json.loads(Path(args.manifest).read_text())
    samples = load_samples(manifest, args.run)
    provenance, data, lane_data, source_hashes = {}, {}, {}, {}
    root = Path(args.metrics)
    Path('qc').mkdir(); Path('provenance').mkdir()
    shutil.copy2(args.logo, 'logo.jpg'); shutil.copy2(args.css, 'report.css')
    for s in samples:
        sid = s['sample_id']; paths = [root / (sid+'_'+q['label']+'.fastp.json') for q in s['fastqs']]
        raw = raw_summary(paths)
        wgs = metric_rows(root / (sid+'.wgs_metrics.txt'))[0]
        insert = metric_rows(root / (sid+'.insert_size_metrics.txt'))
        fr = next((x for x in insert if x.get('PAIR_ORIENTATION') == 'FR'), None)
        if fr is None:
            raise ValueError('No measured FR insert-size metrics: '+sid)
        dup = duplicate_summary(s['bamtools_flags'])
        data[sid] = dict(raw_reads=raw['total_reads'], raw_gb=raw['total_bases']/1e9,
                         raw_q30_pct=raw['q30_pct'], bam_duplicate_pct=dup['bam_duplicate_pct'],
                         mean_insert_size=float(fr['MEAN_INSERT_SIZE']))
        provenance[sid] = dict(source=s, raw=raw, bamtools=dup, picard_wgs=wgs, picard_insert=insert)
        for q,path in zip(s['fastqs'],paths):
            lane = raw_summary([path])
            lane_data[sid+' '+q['label']] = dict(reads=lane['total_reads'], yield_gb=lane['total_bases']/1e9, q30_pct=lane['q30_pct'])
        for path in root.glob(sid+'.*'):
            shutil.copy2(path, Path('qc')/path.name); source_hashes[str(path)] = digest(path)
        for q in s['fastqs']:
            for path in root.glob(sid+'_'+q['label']+'.*'):
                shutil.copy2(path, Path('qc')/path.name); source_hashes[str(path)] = digest(path)
        # BAMTOOLS emits this standard flag-count format; preserve bytes and name its real producer.
        dest=Path('qc')/(sid+'.flagstat.txt');shutil.copy2(s['bamtools_flags'],dest)
        source_hashes[s['bamtools_flags']] = digest(s['bamtools_flags'])
    demux = manifest.get('demultiplexing', {}).get(args.run, [])
    for item in demux:
        path=Path(item);dest=Path('qc/demultiplexing')/path.name;dest.parent.mkdir(exist_ok=True)
        shutil.copy2(path,dest); source_hashes[str(path)]=digest(path)
    description = ('All recorded lane pairs; fastp QC-only (no trimming/filtering or FASTQ rewriting). '
                   'Coverage: native Picard CollectWgsMetrics on all-lane REDUX BAMs, MAPQ20/baseQ20, '
                   'cap250, overlapping mates counted once, primary chr1-22/X/Y territory on the exact HMF reference. '
                   'Duplication/flag plots: Hartwig BAMTOOLS, not Picard MarkDuplicates. '
                   'Fastp duplication is an independent estimate; its disabled filtering does not constitute a laboratory QC pass. ')
    if not demux:
        description += 'Original demultiplexing/PF/undetermined-read statistics were not supplied; FASTQ yield includes assigned reads only. '
    description += manifest.get('run_notes', {}).get(args.run, '')
    headers = {
        'raw_reads': dict(title='Raw reads', description='Individual mate reads, summed across all lane pairs', format='{:,.0f}'),
        'raw_gb': dict(title='Raw Gb', description='FASTQ bases before QC processing / 1e9', format='{:,.2f}'),
        'raw_q30_pct': dict(title='Q30 bases', description='100 * summed raw Q30 bases / summed raw bases', suffix='%', max=100, format='{:,.2f}'),
        'bam_duplicate_pct': dict(title='BAM duplicates', description='BAMTOOLS: 100 * primary duplicate reads / primary reads; PCR and optical combined', suffix='%', max=100, format='{:,.2f}'),
        'mean_insert_size': dict(title='Mean insert', description='Picard CollectInsertSizeMetrics FR mean; bp', suffix=' bp', format='{:,.2f}')}
    config=dict(title='MOHCCN sequencing-run QC', subtitle=args.run,
        report_comment=description, custom_logo='logo.jpg', custom_logo_width=50,
        custom_logo_title="Queen's University Department of Pathology and Molecular Medicine",
        custom_css_files=['report.css'],
        report_header_info=[{'Run ID':args.run},{'Pipeline version':args.release_tag}],
        module_order=[{'custom_content':{}},{'fastp':{'name':'fastp (raw FASTQ QC only)','generalstats':False}},
                      {'samtools':{'name':'Alignment flags (Hartwig BAMTOOLS; samtools format)'}},
                      {'picard':{'name':'Picard (deduplicated coverage and insert size)'}}],
        sample_names_replace_regex=True,
        sample_names_replace={r'\.wgs_metrics$':'',r'\.insert_size_metrics$':'',r'\.fastp$':''},
        custom_data={'lab_summary':dict(plot_type='generalstats',headers=[{k:v} for k,v in headers.items()],data=data),
                     'lab_lane_yield':dict(section_name='Assigned FASTQ yield by sample and lane',
                         description='Individual reads, including both mates; not PF/undetermined demultiplexing counts.',
                         plot_type='table',pconfig={'id':'lab_lane_yield','title':'FASTQ yield by lane'},data=lane_data)},
        table_sample_merge={f'Lane {i:03}':{'type':'regex','pattern':f'_S[0-9]+_L{i:03}$'} for i in range(1,9)})
    # Preserve upstream task-recorded versions, and the versions emitted by these QC tasks.
    versions={}
    for s in samples:
        v=Path(s['source_versions'])
        if v.exists():
            import yaml
            native=yaml.safe_load(v.read_text())
            for process in ['BAMTOOLS','REDUX','BWAMEM2_ALIGN']:
                if process in native: versions[process]=native[process]
            source_hashes[str(v)]=digest(v)
    for p in root.glob('*.version.txt'):
        version=p.read_text().strip()
        if '.fastp.version' in p.name: versions['Raw FASTQ QC']={'fastp':version.removeprefix('fastp ')}
        elif '.picard.version' in p.name: versions['BAM QC']={'Picard':version.removeprefix('Version:')}
    versions['Report']={'MultiQC':Path('multiqc.version.txt').read_text().strip().removeprefix('multiqc, version ')}
    config['software_versions']=versions
    Path('multiqc_config.json').write_text(json.dumps(config,indent=2)+'\n')
    Path('provenance/summary.json').write_text(json.dumps(provenance,indent=2)+'\n')
    Path('provenance/sources.json').write_text(json.dumps(source_hashes,indent=2)+'\n')
    shutil.copy2(args.manifest,'provenance/input_manifest.json')
    Path('provenance/release.json').write_text(json.dumps(dict(tag=args.release_tag,commit=args.release_commit,
        logo_sha256=digest(args.logo),demultiplexing_supplied=bool(demux)),indent=2)+'\n')
    with Path('provenance/sample_summary.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['sample_id',*next(iter(data.values()))]);w.writeheader()
        for sid,row in data.items(): w.writerow(dict(sample_id=sid,**row))


def finalize(args):
    manifest=json.loads(Path(args.manifest).read_text());samples=load_samples(manifest,args.run)
    report=Path(args.run+'.multiqc_report.html');text=report.read_text()
    if args.release_tag not in text: raise ValueError('Missing release label')
    for s in samples:
        if s['sample_id'] not in text: raise ValueError('Missing sample: '+s['sample_id'])
    archives=list(Path('.').glob('*_data.zip'))
    if len(archives)!=1: raise ValueError('Expected one MultiQC data ZIP')
    with zipfile.ZipFile(archives[0]) as z:
        name=next(n for n in z.namelist() if n.endswith('multiqc_general_stats.txt'))
        rows=list(csv.DictReader(z.read(name).decode().splitlines(),delimiter='\t'))
        ids={r['Sample'] for r in rows}
    expected={s['sample_id'] for s in samples}
    if ids!=expected: raise ValueError(f'General Statistics membership differs: {ids} vs {expected}')
    required=['MEAN_COVERAGE','MEDIAN_COVERAGE','picard_target_bases_30X',
              'raw_reads','raw_q30_pct','bam_duplicate_pct','mean_insert_size']
    for suffix in required:
        col=next((k for k in rows[0] if k.endswith('-'+suffix)),None)
        if col is None or any(r.get(col) in ('',None) for r in rows):
            raise ValueError('Missing populated summary metric: '+suffix)
    release=json.loads(Path('provenance/release.json').read_text())
    images=re.findall(r'data:image/[^;]+;base64,([A-Za-z0-9+/=]+)',text)
    if not any(hashlib.sha256(base64.b64decode(x)).hexdigest()==release['logo_sha256'] for x in images):
        raise ValueError('Expected branding image not embedded')
    if not re.search(r'#page_title\s+\.custom_logo\s*\{[^}]*width:\s*50px',text):
        raise ValueError('Missing 50px branding CSS')
    Path('_SUCCESS').write_text(json.dumps(dict(run_id=args.run,samples=sorted(expected),
        release_tag=args.release_tag,release_commit=args.release_commit,report_sha256=digest(report),
        data_sha256=digest(archives[0])),indent=2)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['prepare','finalize'])
    p.add_argument('--manifest',required=True);p.add_argument('--run',required=True)
    p.add_argument('--release-tag',required=True);p.add_argument('--release-commit',required=True)
    p.add_argument('--metrics');p.add_argument('--logo');p.add_argument('--css')
    a=p.parse_args();(prepare if a.mode=='prepare' else finalize)(a)


if __name__=='__main__': main()
