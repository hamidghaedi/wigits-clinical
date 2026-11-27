"""Translate supported tumour SAGE local phase sets for pVACtools 5.4.

Only the somatic caller's tumour namespace establishes heterozygous phase.
Independent germline LPS integers, distance, VAF and normal genotypes cannot.
This writes derivative context/evidence; original calls and selected candidates
are never changed. pVACtools requires the selected anchors in its context VCF.
"""
from bisect import bisect_left
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path

import neo_interface as ni

VERSION = "sage-lps-to-pvac-hp-v1"


def key(fields):
    return tuple(fields[i] for i in (0, 1, 3, 4))


def records(path, tumour):
    samples = ni.vcf_samples(str(path))
    index = samples.index(tumour) if tumour in samples else None
    with ni._open(str(path)) as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            values = ni.sample_values(fields[8], fields[9 + index]) if index is not None else {}
            yield fields, values


def phase_fields(fields, values):
    info = ni.parse_info(fields[7])
    return dict(lps=ni.int_list(info.get('LPS'), default=[]) or [],
                read_counts=ni.int_list(info.get('LPS_RC'), default=[]) or [],
                ad=ni.int_list(values.get('AD'), default=[]) or [],
                gt=values.get('GT', '.'), filter=fields[6])


def heterozygous(gt):
    return gt.replace('|', '/') in ('0/1', '1/0')


def export_context(somatic, germline, phase_source, tumour, anchors, output, evidence, window=100):
    """Return counts; HP is absent (not '.') on unresolved records.

    One source allele with multiple LPS makes each of those labels ambiguous
    for this adapter. We do not flatten overlapping sets into a haplotype.
    Homozygous-alt handling remains pVACtools' own genotype-based behaviour.
    """
    positions = defaultdict(list)
    for chrom, pos, _, _ in anchors:
        positions[chrom].append(int(pos))
    for chrom in positions:
        positions[chrom].sort()

    def near(fields):
        nearby = positions[fields[0]]
        index = bisect_left(nearby, int(fields[1]) - window)
        return index < len(nearby) and nearby[index] <= int(fields[1]) + window

    context = {}
    omitted = []
    anchor_positions = {(a[0], a[1]) for a in anchors}
    missing_tumour = 0
    for origin, path in [('somatic', somatic), ('germline', germline)]:
        if not path:
            continue
        for fields, values in records(path, tumour):
            if fields[6] != 'PASS' or not near(fields):
                continue
            allele = key(fields)
            if allele in context:
                continue  # identical allele retains somatic anchor/context evidence
            if origin == 'germline' and (fields[0], fields[1]) in anchor_positions:
                # The pVAC converter tests main-VCF position, not exact allele,
                # when labelling context origin. A different germline allele
                # here would be incorrectly labelled somatic.
                omitted.append(dict(allele=list(allele), source=str(path),
                    status='AMBIGUOUS_SHARED_POSITION_CONTEXT_NOT_EXPORTED', original=phase_fields(fields, values)))
                continue
            if origin == 'somatic' and allele not in anchors:
                # pVACtools determines somatic/germline context by presence in
                # the main selected VCF. Supplying a non-selected somatic allele
                # would change its WT peptide as if that allele were germline.
                omitted.append(dict(allele=list(allele), source=str(path),
                    status='NONSELECTED_SOMATIC_CONTEXT_NOT_EXPORTED', original=phase_fields(fields, values)))
                continue
            if not values:
                missing_tumour += 1
            context[allele] = (fields, values, origin, str(path))

    source = {}
    duplicates, ambiguous = set(), set()
    for fields, values in records(phase_source, tumour):
        phase = phase_fields(fields, values)
        if len(phase['lps']) > 1:
            ambiguous.update((fields[0], label) for label in phase['lps'])
        allele = key(fields)
        if allele in context:
            if allele in source:
                duplicates.add(allele)
            source[allele] = phase

    headers = []
    # Union definitions: germline FORMAT/INFO declarations may differ from somatic.
    seen = set()
    for path in (somatic, germline):
        if not path:
            continue
        with ni._open(str(path)) as stream:
            for line in stream:
                if line.startswith('#CHROM'):
                    break
                if line.startswith('##FORMAT=<ID=HP,'):
                    continue
                token = line.split(',')[0] if line.startswith(('##INFO=<', '##FORMAT=<', '##FILTER=<', '##contig=<')) else line
                if token not in seen:
                    headers.append(line)
                    seen.add(token)
    headers += ['##FORMAT=<ID=HP,Number=1,Type=String,Description="Opaque same-alt tumour SAGE phase-group token; not parental phase">\n',
                '##FORMAT=<ID=SW_ORIGINAL_HP,Number=.,Type=String,Description="Unmodified source HP; not used across caller namespaces">\n',
                f'##neoPhasingAdapter={VERSION}\n',
                '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + tumour + '\n']
    ledger, counts = [], Counter()
    with open(output, 'w') as stream:
        stream.writelines(headers)
        for allele, (fields, values, origin, path) in context.items():
            native = phase_fields(fields, values)
            phase = source.get(allele)
            hp = None
            if not values:
                status = 'NO_TUMOUR_GENOTYPE'
            elif ',' in fields[4] or not heterozygous(native['gt']):
                status = 'NATIVE_GENOTYPE_ONLY'
            elif allele in duplicates:
                status = 'DUPLICATE_PHASE_SOURCE'
            elif not phase:
                status = 'NO_SOMATIC_PHASE_RECORD'
            elif not (heterozygous(phase['gt']) or phase['gt'] in ('.','./.','.|.')):
                status = 'PHASE_SOURCE_GENOTYPE_UNRESOLVED'
            elif len(phase['ad']) != 2 or phase['ad'][1] <= 0:
                status = 'NO_TUMOUR_ALT_READ_SUPPORT'
            elif len(phase['lps']) != 1 or (fields[0], phase['lps'][0]) in ambiguous:
                status = 'AMBIGUOUS_OR_ABSENT_LPS'
            elif len(phase['read_counts']) != 1 or phase['read_counts'][0] <= 0:
                status = 'NO_POSITIVE_PHASE_READ_SUPPORT'
            else:
                status = 'SUPPORTED_SOMATIC_LPS'
                hp = f"SAGE_{fields[0]}_{phase['lps'][0]}"
            # Crucial: missing HP must be omitted from FORMAT, otherwise pVAC
            # can compare two missing values as equal and infer a false link.
            fmt = [name for name in fields[8].split(':') if name != 'HP']
            written = dict(values)
            if 'HP' in values:
                fmt.append('SW_ORIGINAL_HP')
                written['SW_ORIGINAL_HP'] = values['HP']
            if hp:
                fmt.append('HP')
                written['HP'] = hp
            stream.write('\t'.join(fields[:8] + [':'.join(fmt), ':'.join(written.get(name, '.') for name in fmt)]) + '\n')
            counts[status] += 1
            counts['anchors' if allele in anchors else 'proximal_' + origin] += 1
            ledger.append(dict(allele=list(allele), anchor=allele in anchors, origin=origin,
                               source=path, original=native, original_hp=values.get('HP'),
                               tumour_phase_source=str(phase_source), phase=phase,
                               status=status, hp=hp))
    summary = dict(adapter=VERSION, tumour_sample=tumour, phase_source=str(phase_source),
                   source_kind='native_sage_somatic' if phase_source != somatic else 'purple_somatic_fallback',
                   window_bp=window, counts=dict(counts), missing_tumour_genotypes=missing_tumour,
                   omitted_context=omitted,
                   limitations=['PASS context only; filtered calls are not promoted',
                                'Multiple/overlapping LPS unresolved; different/absent LPS does not prove trans',
                                'Independent germline LPS cannot link to somatic LPS',
                                'Non-selected somatic context is recorded but not mislabelled as germline by pVACtools',
                                'pVACtools 5.4 applies only same-transcript missense proximal consequences',
                                'HP export is not proof of peptide uptake; inspect pVAC proximal tables'],
                   records=ledger)
    Path(evidence).write_text(json.dumps(summary, indent=2) + '\n')
    return dict(counts)


def summarise_usage(predictions, evidence, output):
    """Record actual pVAC converter uptake, separately from available phase.

    A converter row proves context acceptance. Peptide-window application is
    checked separately by the pinned-tool integration test, not inferred here.
    """
    source = json.loads(Path(evidence).read_text())
    tables = []
    for path in sorted(Path(predictions).rglob('*.proximal_variants.tsv')):
        with path.open() as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))
        tables.append(dict(path=str(path), rows=len(rows),
                           main_somatic_variants=sorted({r['main_somatic_variant'] for r in rows}),
                           types=dict(Counter(r['type'] for r in rows))))
    result = dict(adapter=source['adapter'], evidence=str(evidence),
                  source_phase_counts=source['counts'], converter_tables=tables,
                  status='RECORDED' if tables else 'NO_CONVERTER_TABLES',
                  interpretation='Converter uptake; not proof that every context falls inside every predicted peptide')
    Path(output).write_text(json.dumps(result,indent=2)+'\n')
    return result
