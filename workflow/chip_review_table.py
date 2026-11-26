#!/usr/bin/env python3
"""Write the review table for the clonal haematopoiesis calls of one case.

The variant table from chip_candidates.py says what was called and how it was classified. This step
hands over the evidence a reviewer needs to check each call without opening the VCF: depth and
supporting reads in blood and tumour, the read-context breakdown, base and mapping quality, strand and
fragment orientation, position within the read, repeat context, and the outcome of each published
filtering criterion taken one at a time.

One row per variant, one value per cell. Multi-valued VCF fields are split into named columns.

Outputs into --outdir:
  <case>.chip_review.tsv            one row per variant, every variant
  <case>.chip_review_columns.tsv    one row per column: where the value comes from and what it means

The criteria columns apply the thresholds of Vlasschaert et al., Blood 2023;141(18):2214-2223, to SAGE
calls. Their pipeline calls with Mutect2, so these columns are our application of their criteria and
are not a reproduction of their output.
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chip_candidates import _open, parse_info, parse_impact  # noqa: E402
from chip_curation import CuratedDefinition  # noqa: E402

# Published sequence criteria.
MIN_DEPTH = 20
MIN_ALT_READS = 5
MIN_ALLELE_FRACTION = 0.02
HOMOPOLYMER_MIN_LENGTH = 3

RC_PARTS = ("full", "partial_core", "core", "realigned", "reference", "total")
UMI_PARTS = ("total_none", "total_single", "total_dual_strand", "alt_none", "alt_single",
             "alt_dual_strand")


def split(value, n, cast=float):
    """Split a comma-separated FORMAT value into n cells; blanks where the value is absent."""
    parts = (value or "").split(",")
    out = []
    for i in range(n):
        try:
            v = cast(parts[i])
            out.append(int(v) if cast is int else v)
        except (IndexError, ValueError, TypeError):
            out.append("")
    return out


def wilson(k, n, z=1.959964):
    if not n:
        return "", ""
    p = k / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return max(0.0, centre - half), min(1.0, centre + half)


def p_if_heterozygote(k, n):
    """Two-sided exact binomial probability of a count this far from n/2 when the true fraction is 0.5."""
    if not n:
        return ""
    obs = math.comb(n, k)
    total = sum(math.comb(n, i) for i in range(n + 1) if math.comb(n, i) <= obs)
    return total / (2 ** n)


def variant_type(ref, alt):
    if len(ref) == len(alt):
        return "SNV" if len(ref) == 1 else "MNV"
    return "insertion" if len(alt) > len(ref) else "deletion"


def sample_cells(prefix, g):
    """Flatten one sample's FORMAT fields into named cells."""
    ad = split(g.get("AD"), 2, int)
    c = {f"{prefix}_depth": split(g.get("DP"), 1, int)[0],
         f"{prefix}_ref_reads": ad[0],
         f"{prefix}_alt_reads": ad[1],
         f"{prefix}_allele_fraction": split(g.get("AF"), 1)[0]}
    for name, v in zip(RC_PARTS, split(g.get("RC_CNT"), 6, int)):
        c[f"{prefix}_context_reads_{name}"] = v
    for name, v in zip(RC_PARTS, split(g.get("RC_QUAL"), 6, int)):
        c[f"{prefix}_context_quality_{name}"] = v
    c[f"{prefix}_base_quality_recalibrated_all"], c[f"{prefix}_base_quality_recalibrated_alt"] = \
        split(g.get("ARCBQ"), 2, int)
    c[f"{prefix}_base_quality_sequencer_alt"] = split(g.get("ASBQ"), 1, int)[0]
    c[f"{prefix}_base_quality_final"] = split(g.get("AFBQ"), 1, int)[0]
    c[f"{prefix}_map_quality_all"], c[f"{prefix}_map_quality_alt"] = split(g.get("ARMQ"), 2, int)
    c[f"{prefix}_map_quality_alt_final"] = split(g.get("AFMQ"), 1, int)[0]
    c[f"{prefix}_read_forward_fraction_ref"], c[f"{prefix}_read_forward_fraction_alt"] = \
        split(g.get("RSB"), 2)
    c[f"{prefix}_fragment_forward_fraction_ref"], c[f"{prefix}_fragment_forward_fraction_alt"] = \
        split(g.get("FSB"), 2)
    c[f"{prefix}_edge_distance_alt"], c[f"{prefix}_edge_distance_all"] = split(g.get("AED"), 2)
    c[f"{prefix}_unique_fragment_positions_alt"] = split(g.get("MUC"), 1, int)[0]
    c[f"{prefix}_improper_pairs"] = split(g.get("RC_IPC"), 1, int)[0]
    j = split(g.get("RC_JIT"), 3, int)
    c[f"{prefix}_jitter_shortened"], c[f"{prefix}_jitter_lengthened"] = j[0], j[1]
    c[f"{prefix}_jitter_quality_penalty"] = j[2]
    c[f"{prefix}_simple_alt_matches"] = split(g.get("SAC"), 1, int)[0]
    for name, v in zip(UMI_PARTS, split(g.get("UMI_CNT"), 6, int)):
        c[f"{prefix}_fragments_{name}"] = v
    return c


def both(fraction):
    if fraction == "":
        return ""
    return "yes" if 0.0 < fraction < 1.0 else "no"


def sample_columns(prefix, label):
    s = f"SAGE FORMAT, {label} sample"
    d = [(f"{prefix}_depth", f"{s}: DP", f"Read depth in {label} at this position, duplicates removed."),
         (f"{prefix}_ref_reads", f"{s}: AD", f"Reads in {label} supporting the reference allele."),
         (f"{prefix}_alt_reads", f"{s}: AD", f"Reads in {label} supporting the variant allele."),
         (f"{prefix}_allele_fraction", f"{s}: AF",
          f"Variant allele fraction in {label}, from read-context support over coverage.")]
    meaning = {"full": "matching the full read context, core and both flanks",
               "partial_core": "matching the core and part of the flanks",
               "core": "matching the core only",
               "realigned": "matching after local realignment",
               "reference": "matching the reference context",
               "total": "counted at this position"}
    for p in RC_PARTS:
        d.append((f"{prefix}_context_reads_{p}", f"{s}: RC_CNT", f"Reads in {label} {meaning[p]}."))
    for p in RC_PARTS:
        d.append((f"{prefix}_context_quality_{p}", f"{s}: RC_QUAL",
                  f"Summed quality of the reads in {label} {meaning[p]}."))
    d += [
        (f"{prefix}_base_quality_recalibrated_all", f"{s}: ARCBQ",
         f"Mean recalibrated base quality at this position in {label}, all reads."),
        (f"{prefix}_base_quality_recalibrated_alt", f"{s}: ARCBQ",
         f"Mean recalibrated base quality at this position in {label}, variant reads."),
        (f"{prefix}_base_quality_sequencer_alt", f"{s}: ASBQ",
         f"Mean base quality as reported by the sequencer, variant reads in {label}."),
        (f"{prefix}_base_quality_final", f"{s}: AFBQ",
         f"Mean base quality in {label} after all SAGE adjustments."),
        (f"{prefix}_map_quality_all", f"{s}: ARMQ", f"Mean mapping quality in {label}, all reads."),
        (f"{prefix}_map_quality_alt", f"{s}: ARMQ", f"Mean mapping quality in {label}, variant reads."),
        (f"{prefix}_map_quality_alt_final", f"{s}: AFMQ",
         f"Mean mapping quality of variant reads in {label} after all SAGE adjustments."),
        (f"{prefix}_read_forward_fraction_ref", f"{s}: RSB",
         f"Fraction of reference reads in {label} on the forward strand."),
        (f"{prefix}_read_forward_fraction_alt", f"{s}: RSB",
         f"Fraction of variant reads in {label} on the forward strand. Values at 0 or 1 mean every "
         f"variant read is on one strand."),
        (f"{prefix}_fragment_forward_fraction_ref", f"{s}: FSB",
         f"Fraction of reference fragments in {label} in forward orientation."),
        (f"{prefix}_fragment_forward_fraction_alt", f"{s}: FSB",
         f"Fraction of variant fragments in {label} in forward orientation. Values at 0 or 1 mean "
         f"every variant fragment has one orientation."),
        (f"{prefix}_edge_distance_alt", f"{s}: AED",
         f"Mean distance of the variant from the read end in {label}, as a fraction of read length, "
         f"variant reads. Low values mean the variant sits at read ends."),
        (f"{prefix}_edge_distance_all", f"{s}: AED",
         f"Mean distance of this position from the read end in {label}, all reads."),
        (f"{prefix}_unique_fragment_positions_alt", f"{s}: MUC",
         f"Minimum count of distinct fragment start or end coordinates among variant reads in "
         f"{label}. A low count with several reads suggests copies of one fragment."),
        (f"{prefix}_improper_pairs", f"{s}: RC_IPC", f"Variant reads in {label} in improper pairs."),
        (f"{prefix}_jitter_shortened", f"{s}: RC_JIT",
         f"Reads in {label} carrying a shortened repeat at this site."),
        (f"{prefix}_jitter_lengthened", f"{s}: RC_JIT",
         f"Reads in {label} carrying a lengthened repeat at this site."),
        (f"{prefix}_jitter_quality_penalty", f"{s}: RC_JIT",
         f"Quality penalty applied in {label} for repeat jitter."),
        (f"{prefix}_simple_alt_matches", f"{s}: SAC",
         f"Reads in {label} carrying the variant base without the full read context."),
    ]
    um = {"total_none": "all fragments seen once", "total_single": "all fragments built from copies",
          "total_dual_strand": "all fragments with both strands observed",
          "alt_none": "variant fragments seen once", "alt_single": "variant fragments built from copies",
          "alt_dual_strand": "variant fragments with both strands observed"}
    for p in UMI_PARTS:
        d.append((f"{prefix}_fragments_{p}", f"{s}: UMI_CNT", f"Count in {label} of {um[p]}."))
    return d


HEAD = [
    ("case", "run", "Case identifier."),
    ("gene", "PAVE IMPACT", "Gene symbol."),
    ("chrom", "VCF", "Chromosome, GRCh38."),
    ("pos", "VCF", "Position, GRCh38, 1-based."),
    ("ref", "VCF", "Reference allele."),
    ("alt", "VCF", "Variant allele."),
    ("variant_type", "derived", "SNV, MNV, insertion or deletion."),
    ("hgvs_coding", "PAVE IMPACT", "Coding change on the canonical transcript."),
    ("hgvs_protein", "PAVE IMPACT", "Protein change on the canonical transcript."),
    ("consequence", "PAVE IMPACT", "Consequence on the canonical transcript."),
    ("coding_effect", "PAVE IMPACT", "Coding effect class."),
    ("transcript", "PAVE IMPACT", "Canonical Ensembl transcript."),
    ("refseq_id", "PAVE PAVE_TI", "RefSeq identifier of the canonical transcript."),
    ("exon", "PAVE PAVE_TI", "Exon number on the canonical transcript."),
    ("codon", "PAVE PAVE_TI", "Codon number on the canonical transcript."),
    ("curated_accession", "curated definition", "RefSeq accession the curated rule for this gene is written against."),
    ("transcript_matches_curated_accession", "derived",
     "Whether the canonical transcript is the one the curated rule is written against. When it is "
     "not, residue numbering can differ and a curated position may be missed."),
    ("classification", "chip_candidates", "Clone candidate or constitutional."),
    ("confidence_band", "chip_candidates", "Confident, supported or provisional."),
    ("in_panel", "chip_candidates", "Whether the gene is in the CHIP panel."),
    ("curated_gene_status", "curated definition",
     "retained: on the reference study's list. candidate_not_retained: a candidate gene the study "
     "did not keep. not_curated: outside their candidate list."),
    ("gene_age_association_ukb", "curated definition",
     "Association of mutations in this gene with age in UK Biobank at five or more supporting "
     "reads, odds ratio per year with significance."),
    ("in_curated_definition", "chip_curation",
     "Whether the gene is retained and the change is one the curated rule recognises in that gene."),
    ("curation_note", "chip_curation", "Reason for the value in in_curated_definition."),
    ("known_artefact", "chip_curation", "Matches a recurrent artefact or a binomial-test failure in the reference cohorts."),
    ("cohort_validated", "chip_curation", "Associated with age and the TERT risk allele in the reference cohorts."),
    ("meets_depth", "derived", f"Blood depth of {MIN_DEPTH} or more."),
    ("meets_alt_reads", "derived", f"{MIN_ALT_READS} or more supporting reads in blood."),
    ("meets_allele_fraction", "derived", f"Blood allele fraction of {MIN_ALLELE_FRACTION:.0%} or more."),
    ("meets_orientation", "derived",
     "Variant fragments in blood seen in both orientations, taken from SAGE fragment strand bias "
     "as the counterpart of the F1R2 and F2R1 requirement."),
    ("clears_homopolymer_rule", "derived",
     f"No for an insertion or deletion in a single-base repeat of {HOMOPOLYMER_MIN_LENGTH} or more."),
    ("meets_all_sequence_criteria", "derived", "All five sequence criteria met."),
    ("published_criteria_outcome", "derived",
     "retained when the call meets every sequence criterion, is inside the curated definition, is "
     "not classified as constitutional and matches no known artefact. Otherwise the criteria it "
     "misses, in full."),
    ("blood_af_ci95_low", "derived", "Lower bound of the 95 percent Wilson interval on the blood allele fraction."),
    ("blood_af_ci95_high", "derived", "Upper bound of the 95 percent Wilson interval on the blood allele fraction."),
    ("detection_floor_at_this_depth", "chip_candidates",
     "Allele fraction at which five reads are seen with 95 percent power at this depth."),
    ("p_if_inherited_heterozygote", "derived",
     "Two-sided exact binomial probability of the blood read counts if the variant were inherited "
     "and heterozygous. Small values argue against an inherited heterozygous variant."),
    ("blood_alt_on_both_strands", "derived", "Variant reads in blood seen on both strands."),
    ("present_in_tumour", "chip_candidates", "Whether the variant has support in the tumour sample."),
    ("tumour_to_blood_af_ratio", "derived",
     "Tumour allele fraction over blood allele fraction. Near or above 1 at a low blood fraction is "
     "what tumour material in the blood draw or an infiltrating clone would give."),
]
TAIL = [
    ("sage_tier", "SAGE INFO: TIER", "HOTSPOT, PANEL, HIGH_CONFIDENCE or LOW_CONFIDENCE. Thresholds are most permissive at hotspots."),
    ("sage_quality", "VCF QUAL", "SAGE variant quality."),
    ("sage_filter", "SAGE FILTER", "Filters SAGE applied, kept as evidence. PASS when none."),
    ("sage_map_quality_factor", "SAGE INFO: MQF", "Mapping-quality heuristic used in the SAGE quality filter."),
    ("sage_quality_probability", "SAGE INFO: TQP", "Probability used in the SAGE quality filter for the sample being called, blood here."),
    ("max_edge_distance", "SAGE INFO: MED", "Largest distance of the variant from a read end among supporting reads, as a fraction of read length."),
    ("min_events_in_read", "SAGE INFO: RC_NM", "Fewest mismatch or indel events in any supporting read."),
    ("repeat_count", "SAGE INFO: RC_REPC", "Repeat count in the read context."),
    ("repeat_sequence", "SAGE INFO: RC_REPS", "Repeat unit in the read context."),
    ("reference_repeat_count", "SAGE INFO: REP_C", "Repeat count in the reference at this site."),
    ("reference_repeat_sequence", "SAGE INFO: REP_S", "Repeat unit in the reference at this site."),
    ("microhomology", "SAGE INFO: MH or RC_MH", "Microhomology at an insertion or deletion."),
    ("nearby_indel", "SAGE INFO: NEARBY_INDEL", "An insertion or deletion overlaps the core of this variant."),
    ("local_phase_set", "SAGE INFO: LPS", "Identifier shared by variants observed on the same reads."),
    ("local_phase_set_reads", "SAGE INFO: LPS_RC", "Reads supporting the phase set."),
    ("trinucleotide", "SAGE INFO: TNC", "Reference trinucleotide context."),
    ("read_context_left_flank", "SAGE INFO: RC_INFO", "Left flank of the read context."),
    ("read_context_core", "SAGE INFO: RC_INFO", "Core of the read context."),
    ("read_context_right_flank", "SAGE INFO: RC_INFO", "Right flank of the read context."),
    ("mappability", "PAVE INFO: MAPPABILITY", "GEM mappability in a 150-base window; 1 is uniquely mappable."),
    ("gnomad_frequency", "PAVE INFO: GND_FREQ", "Population allele frequency in gnomAD."),
    ("clinvar", "PAVE INFO: CLNSIG", "ClinVar significance."),
    ("clinvar_conflict", "PAVE INFO: CLNSIGCONF", "ClinVar conflicting interpretations."),
    ("warnings", "chip_candidates", "Every warning raised for this call."),
]
COLUMNS = (HEAD[:39] + sample_columns("blood", "blood")
           + HEAD[39:] + sample_columns("tumour", "tumour") + TAIL)


def transcript_details(info, transcript):
    for entry in str(info.get("PAVE_TI", "")).split(","):
        f = entry.split("|")
        if len(f) >= 10 and f[2] == transcript:
            return f[7], f[8], f[9]
    return "", "", ""


def build(vcf, variants_tsv, outdir, case, blood_sample, tumour_sample, curated_genes=None):
    os.makedirs(outdir, exist_ok=True)
    with open(variants_tsv) as fh:
        called = list(csv.DictReader(fh, delimiter="\t"))
    order = {(r["chrom"], r["pos"], r["ref"], r["alt"]): i for i, r in enumerate(called)}
    by_key = {(r["chrom"], r["pos"], r["ref"], r["alt"]): r for r in called}
    curation = CuratedDefinition(curated_genes) if curated_genes and os.path.isfile(curated_genes) else None

    rows, samples = [], []
    with _open(vcf) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            f = line.rstrip("\n").split("\t")
            key = (f[0], f[1], f[3], f[4])
            c = by_key.get(key)
            if c is None:
                raise SystemExit(f"{key} is in the VCF and missing from {variants_tsv}")
            info = parse_info(f[7])
            fmt = f[8].split(":")
            geno = {s: dict(zip(fmt, f[9 + i].split(":"))) for i, s in enumerate(samples)}
            imp = parse_impact(info.get("IMPACT", ""))
            row = {k: c.get(k, "") for k in (
                "case", "gene", "chrom", "pos", "ref", "alt", "hgvs_coding", "hgvs_protein",
                "consequence", "coding_effect", "transcript", "classification", "confidence_band",
                "in_panel", "curated_gene_status", "gene_age_association_ukb",
                "in_curated_definition", "curation_note", "known_artefact", "cohort_validated",
                "detection_floor_at_this_depth", "present_in_tumour", "sage_filter", "warnings")}
            row["variant_type"] = variant_type(f[3], f[4])
            row["refseq_id"], row["exon"], row["codon"] = transcript_details(info, imp.get("transcript", ""))
            acc = curation.genes.get(row["gene"], {}).get("accession", "") if curation else ""
            row["curated_accession"] = acc
            row["transcript_matches_curated_accession"] = (
                "" if not acc or not row["refseq_id"] else
                "yes" if row["refseq_id"].split(".")[0] == acc.split(".")[0] else "no")

            row.update(sample_cells("blood", geno.get(blood_sample, {})))
            row.update(sample_cells("tumour", geno.get(tumour_sample, {})))

            dp, alt, af = row["blood_depth"], row["blood_alt_reads"], row["blood_allele_fraction"]
            lo, hi = wilson(alt, dp) if dp != "" and alt != "" else ("", "")
            row["blood_af_ci95_low"], row["blood_af_ci95_high"] = lo, hi
            row["p_if_inherited_heterozygote"] = (
                p_if_heterozygote(min(alt, dp), dp) if dp != "" and alt != "" else "")
            row["blood_alt_on_both_strands"] = both(row["blood_read_forward_fraction_alt"])
            taf = row["tumour_allele_fraction"]
            row["tumour_to_blood_af_ratio"] = (taf / af) if af not in ("", 0) and taf != "" else ""

            reps, repc = info.get("RC_REPS", ""), info.get("RC_REPC", "")
            is_indel = len(f[3]) != len(f[4])
            homopolymer = (is_indel and isinstance(reps, str) and len(reps) == 1
                           and str(repc).isdigit() and int(repc) >= HOMOPOLYMER_MIN_LENGTH)
            checks = [
                ("meets_depth", dp != "" and dp >= MIN_DEPTH, f"depth below {MIN_DEPTH}"),
                ("meets_alt_reads", alt != "" and alt >= MIN_ALT_READS,
                 f"fewer than {MIN_ALT_READS} supporting reads"),
                ("meets_allele_fraction", af != "" and af >= MIN_ALLELE_FRACTION,
                 f"allele fraction below {MIN_ALLELE_FRACTION:.0%}"),
                ("meets_orientation", both(row["blood_fragment_forward_fraction_alt"]) == "yes",
                 "variant fragments in one orientation only"),
                ("clears_homopolymer_rule", not homopolymer, "insertion or deletion in a homopolymer"),
            ]
            misses = []
            for name, ok, why in checks:
                row[name] = "yes" if ok else "no"
                if not ok:
                    misses.append(why)
            row["meets_all_sequence_criteria"] = "no" if misses else "yes"
            if row["in_curated_definition"] == "no":
                misses.append("outside the curated definition")
            if row["classification"] != "clone candidate":
                misses.append("classified as constitutional")
            if row["known_artefact"] == "yes":
                misses.append("matches a known artefact")
            row["published_criteria_outcome"] = "retained" if not misses else "; ".join(misses)

            row["sage_tier"] = info.get("TIER", "")
            row["sage_quality"] = f[5]
            row["sage_map_quality_factor"] = info.get("MQF", "")
            row["sage_quality_probability"] = info.get("TQP", "")
            row["max_edge_distance"] = info.get("MED", "")
            row["min_events_in_read"] = info.get("RC_NM", "")
            row["repeat_count"], row["repeat_sequence"] = repc, reps
            row["reference_repeat_count"] = info.get("REP_C", "")
            row["reference_repeat_sequence"] = info.get("REP_S", "")
            row["microhomology"] = info.get("MH", "") or info.get("RC_MH", "")
            row["nearby_indel"] = "yes" if info.get("NEARBY_INDEL") else "no"
            row["local_phase_set"] = info.get("LPS", "")
            row["local_phase_set_reads"] = info.get("LPS_RC", "")
            row["trinucleotide"] = info.get("TNC", "")
            rc = str(info.get("RC_INFO", "")).split("-")
            row["read_context_left_flank"] = rc[2] if len(rc) > 4 else ""
            row["read_context_core"] = rc[3] if len(rc) > 4 else ""
            row["read_context_right_flank"] = rc[4] if len(rc) > 4 else ""
            row["mappability"] = info.get("MAPPABILITY", "")
            row["gnomad_frequency"] = info.get("GND_FREQ", "")
            row["clinvar"] = info.get("CLNSIG", "")
            row["clinvar_conflict"] = info.get("CLNSIGCONF", "")
            rows.append((order[key], row))

    names = [c[0] for c in COLUMNS]
    if len(set(names)) != len(names):
        raise SystemExit("duplicate column names in the review table definition")
    rows.sort(key=lambda x: x[0])
    for _, row in rows:
        extra, missing = set(row) - set(names), set(names) - set(row)
        if extra or missing:
            raise SystemExit(f"review table definition is out of step: extra {sorted(extra)}, "
                             f"missing {sorted(missing)}")

    def cell(v):
        if isinstance(v, float):
            return f"{v:.3g}" if 0 < abs(v) < 1e-3 else f"{v:.4f}".rstrip("0").rstrip(".")
        return "" if v is True else str(v)

    table = os.path.join(outdir, f"{case}.chip_review.tsv")
    with open(table, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(names)
        for _, row in rows:
            w.writerow([cell(row[n]) for n in names])
    dictionary = os.path.join(outdir, f"{case}.chip_review_columns.tsv")
    with open(dictionary, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["column", "source", "description"])
        w.writerows(COLUMNS)
    if len(rows) != len(called):
        raise SystemExit(f"{len(called)} variants in the table and {len(rows)} in the VCF")
    return table, dictionary, [r for _, r in rows]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annotated-vcf", required=True)
    ap.add_argument("--variants-tsv", required=True, help="<case>.chip_variants.tsv from chip_candidates.py")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--case", required=True)
    ap.add_argument("--blood-sample", required=True)
    ap.add_argument("--tumour-sample", required=True)
    ap.add_argument("--curated-genes", default=None)
    a = ap.parse_args(argv)
    table, dictionary, rows = build(a.annotated_vcf, a.variants_tsv, a.outdir, a.case,
                                    a.blood_sample, a.tumour_sample, a.curated_genes)
    print(f"variants {len(rows)}  columns {len(COLUMNS)}")
    print(f"meet all sequence criteria {sum(r['meets_all_sequence_criteria'] == 'yes' for r in rows)}")
    print(f"retained under the published criteria {sum(r['published_criteria_outcome'] == 'retained' for r in rows)}")
    print(table)
    print(dictionary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
