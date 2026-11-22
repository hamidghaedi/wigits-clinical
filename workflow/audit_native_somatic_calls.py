#!/usr/bin/env python3
"""Read-only native SAGE/PURPLE/ORANGE QC; not a filter or validation of accuracy.

Inputs are explicit case/result paths in a JSON list. Stream VCFs on Slurm;
count VCF records, never equate these to independent mutations or PURPLE TMB.
Only FILTER=PASS counts as passing. Preserve native AF separately from AD/DP.
Diagnostic counters (not rejection rules): ALT <5; normal ALT >=3 with AF >=2%;
normal DP <20; >=6 ALT with FSB/RSB exactly 0/1; AFBQ <20; AFMQ <30; MUC <3.
Missing metrics remain missing. C>T/G>A alone never establishes FFPE artefact.
TCW C>T/C>G fractions describe context; they are not a signature refit.
ORANGE's somaticVariants are catalog entries, not all necessarily reported:
retain the reported flag, transcripts and exact allele key. No BAM is read,
no caller runs, and no source report/result is modified. Internal PASS-site
tables support cross-case recurrence checks, not new reporting selection.
"""
import argparse
import collections
import csv
import datetime
import gzip
import json
from pathlib import Path
import statistics

VERSION = '0.1.0'
COMPLEMENT = str.maketrans('ACGT', 'TGCA')


def number(value):
    try:
        return float(value) if value not in (None, '.', '') else None
    except ValueError:
        return None


def first(value):
    return number(value.split(',')[0]) if value else None


def support(sample):
    ad = sample.get('AD', '').split(',')
    dp = first(sample.get('DP'))
    alt = number(ad[1]) if len(ad) == 2 else None
    result = {'ref': number(ad[0]) if len(ad) == 2 else None,
              'alt': alt, 'dp': dp, 'native_af': first(sample.get('AF')),
              'alt_over_dp': alt / dp if alt is not None and dp else None}
    for field in ('AFBQ', 'AFMQ', 'MUC'):
        result[field] = first(sample.get(field))
    for field in ('FSB', 'RSB', 'AED', 'RC_CNT'):
        result[field] = sample.get(field)
    return result


def diagnostics(tumor, normal):
    flags = []
    if tumor['alt'] is not None and tumor['alt'] < 5:
        flags.append('tumor_alt_lt5')
    if normal['dp'] is not None and normal['dp'] < 20:
        flags.append('normal_dp_lt20')
    if (normal['alt'] is not None and normal['alt'] >= 3 and
            normal['native_af'] is not None and normal['native_af'] >= .02):
        flags.append('normal_alt_ge3_af_ge2pct')
    for field in ('FSB', 'RSB'):
        values = (tumor[field] or '').split(',')
        if len(values) == 2 and tumor['alt'] is not None and tumor['alt'] >= 6:
            if number(values[1]) in (0, 1):
                flags.append(field.lower() + '_alt_exclusive')
    for field, threshold in [('AFBQ', 20), ('AFMQ', 30), ('MUC', 3)]:
        if tumor[field] is not None and tumor[field] < threshold:
            flags.append(field.lower() + '_below_' + str(threshold))
    return flags


def substitution(ref, alt, context):
    if len(ref) != 1 or len(alt) != 1 or ref not in 'ACGT' or alt not in 'ACGT':
        return None, None
    if ref in 'AG':
        ref, alt = ref.translate(COMPLEMENT), alt.translate(COMPLEMENT)
        context = context.translate(COMPLEMENT)[::-1]
    label = ref + '>' + alt
    tcw = (len(context) == 3 and context[1] == ref and context[0] == 'T'
           and context[2] in 'AT' and label in ('C>T', 'C>G'))
    return label, tcw


def allele(chrom, pos, ref, alt):
    return ':'.join([chrom, str(pos), ref, alt])


def scan_vcf(path, tumor_id, normal_id, selected, sites=None):
    counts = collections.Counter()
    filters = collections.Counter()
    diags = collections.Counter()
    changes = collections.Counter()
    low_changes = collections.Counter()
    hist = collections.Counter()
    missing = collections.Counter()
    found = {}
    afs = []
    header = []
    columns = None
    with gzip.open(path, 'rt') as stream:
        for line in stream:
            if line.startswith('##'):
                if any(x in line.lower() for x in ('tinc', 'version', 'source=')):
                    header.append(line.rstrip())
                continue
            if line.startswith('#CHROM'):
                columns = line.rstrip().split('\t')
                ti, ni = columns.index(tumor_id), columns.index(normal_id)
                continue
            row = line.rstrip().split('\t')
            counts['all_records'] += 1
            filters[row[6]] += 1
            key = allele(row[0], row[1], row[3], row[4])
            passed = row[6] == 'PASS'
            if not passed and key not in selected:
                continue
            fmt = row[8].split(':')
            tumor = support(dict(zip(fmt, row[ti].split(':'))))
            normal = support(dict(zip(fmt, row[ni].split(':'))))
            flags = diagnostics(tumor, normal)
            info = dict(x.split('=', 1) if '=' in x else (x, True)
                        for x in row[7].split(';'))
            record = {'allele': key, 'filter': row[6], 'qual': number(row[5]),
                      'tier': info.get('TIER'), 'gene': info.get('GENE'),
                      'trinucleotide': info.get('TNC'), 'tumor': tumor,
                      'normal': normal, 'diagnostic_flags': flags,
                      'native_hotspot': info.get('HOTSPOT'),
                      'purple_cn': info.get('PURPLE_CN'),
                      'purple_macn': info.get('PURPLE_MACN')}
            if key in selected:
                found[key] = record
            if not passed:
                continue
            counts['pass_records'] += 1
            if ',' in row[4]:
                counts['multiallelic_records'] += 1
            label, tcw = substitution(row[3], row[4], info.get('TNC', ''))
            if label:
                counts['pass_snv'] += 1
                changes[label] += 1
                if tcw:
                    counts['tcw_ct_or_cg'] += 1
            else:
                counts['pass_non_snv'] += 1
            af = tumor['native_af']
            if af is None:
                hist['missing'] += 1
            else:
                afs.append(af)
                band = ('<0.02' if af < .02 else '0.02-<0.05' if af < .05 else
                        '0.05-<0.10' if af < .10 else '0.10-<0.20' if af < .20 else
                        '0.20-<0.50' if af < .50 else '>=0.50')
                hist[band] += 1
                if label and af < .10:
                    low_changes[label] += 1
            diags.update(flags)
            for field, value in tumor.items():
                if value is None:
                    missing['tumor_' + field] += 1
            for field, value in normal.items():
                if value is None:
                    missing['normal_' + field] += 1
            if normal['alt'] is not None and normal['alt'] > 0:
                counts['pass_with_any_normal_alt'] += 1
            if sites:
                sites.writerow([key, af, tumor['alt'], tumor['dp'],
                                normal['alt'], normal['dp'], normal['native_af'],
                                ','.join(flags)])
    return {'counts': dict(counts), 'filter_counts': dict(filters),
            'vaf_bins': dict(hist), 'median_native_vaf': statistics.median(afs) if afs else None,
            'snv_substitutions': dict(changes), 'snv_substitutions_vaf_lt10pct': dict(low_changes),
            'diagnostic_counts': dict(diags), 'missing_metrics': dict(missing),
            'selected_records': found, 'source_header': header, 'samples': columns[9:]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--index', type=int, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    case = json.loads(args.manifest.read_text())[args.index]
    args.out.mkdir(parents=True, exist_ok=True)
    orange = json.loads(Path(case['orange']).read_text())
    selected = {allele(v['chromosome'], v['position'], v['ref'], v['alt']): v
                for v in orange['purple']['somaticVariants']}
    result = {'version': VERSION, 'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'case': case, 'orange_fit': orange['purple']['fit'],
              'orange_characteristics': orange['purple']['characteristics'],
              'orange_catalog_rows': len(orange['purple']['somaticVariants']),
              'orange_reported_rows': sum(v.get('reported') is True for v in orange['purple']['somaticVariants']),
              'orange_catalog': list(selected.values()), 'signatures': orange.get('sigAllocations'),
              'configured_primary_tumor': orange.get('configuredPrimaryTumor')}
    for source in ('sage', 'purple'):
        if source == 'purple':
            with gzip.open(args.out / (case['case_id'] + '.pass_sites.tsv.gz'), 'wt') as fh:
                writer = csv.writer(fh, delimiter='\t')
                writer.writerow(['allele', 'native_af', 'tumor_alt', 'tumor_dp',
                                 'normal_alt', 'normal_dp', 'normal_native_af', 'diagnostic_flags'])
                result[source] = scan_vcf(case[source], orange['sampleId'], orange['referenceId'], selected, writer)
        else:
            result[source] = scan_vcf(case[source], orange['sampleId'], orange['referenceId'], selected)
    (args.out / (case['case_id'] + '.json')).write_text(json.dumps(result, indent=2) + '\n')
    print(case['case_id'], {s: result[s]['counts'] for s in ('sage', 'purple')}, flush=True)


if __name__ == '__main__':
    main()
