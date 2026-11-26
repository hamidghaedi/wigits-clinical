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

Outputs also include <case>.chip_review_summary.json, the counts under each supporting-read threshold
and whether the sample should be inspected as a possible poor library.

The criteria columns apply the thresholds and gene-specific rules of Vlasschaert et al., Blood
2023;141(18):2214-2223, to SAGE calls. Their pipeline calls with Mutect2, so these columns are our
application of their criteria and are not a reproduction of their output.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chip_candidates import _open, parse_info, parse_impact  # noqa: E402
from chip_curation import (CuratedDefinition, short_protein, BINOMIAL_GENES, BINOMIAL_EXEMPT,  # noqa: E402
                           MISSENSE_RANGES, ASXL1_G646W_MIN_ALLELE_FRACTION)

# Sequence criteria of the paper. It sets the supporting-read threshold per cohort from the association
# with age and the TERT risk allele: 5 for the UK Biobank exomes and 3 for the All of Us genomes.
MIN_DEPTH = 20
MIN_ALT_READS_EXOME = 5
MIN_ALT_READS_GENOME = 3
MIN_ALLELE_FRACTION = 0.02
# Applied by the Rauh laboratory in addition; the paper keeps ASXL1 p.G646Wfs*12 despite it.
HOMOPOLYMER_MIN_LENGTH = 3
# Probability below which the read counts are taken to differ from an inherited heterozygote.
BINOMIAL_P = 0.01
# The paper advises inspecting a sample with this many retained variants as a possible poor library.
INSPECT_SAMPLE_AT = 4

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
    ("hgvs_coding", "PAVE IMPACT", "Coding change on the Ensembl canonical transcript."),
    ("hgvs_protein", "PAVE IMPACT", "Protein change on the Ensembl canonical transcript."),
    ("consequence", "PAVE IMPACT", "Consequence on the Ensembl canonical transcript."),
    ("coding_effect", "PAVE IMPACT", "Coding effect class."),
    ("transcript", "PAVE IMPACT", "Ensembl canonical transcript."),
    ("curated_accession", "authors' gene list", "RefSeq accession the curated rule for this gene is written against."),
    ("curated_transcript", "chip_curation", "Ensembl transcript on which the curated rule was evaluated."),
    ("transcript_basis", "chip_curation",
     "curated_accession when that transcript maps to the accession of the rule; canonical_fallback "
     "when no Ensembl transcript does and the canonical one was used."),
    ("curated_hgvs_coding", "PAVE PAVE_TI", "Coding change on the curated transcript."),
    ("curated_hgvs_protein", "PAVE PAVE_TI", "Protein change on the curated transcript. This is the change tested against the whitelist."),
    ("curated_effect", "PAVE PAVE_TI", "Consequence on the curated transcript."),
    ("curated_exon", "PAVE PAVE_TI", "Exon on the curated transcript, Ensembl numbering."),
    ("curated_codon", "PAVE PAVE_TI", "Codon on the curated transcript."),
    ("classification", "chip_candidates", "Clone candidate or constitutional."),
    ("confidence_band", "chip_candidates", "Confident, supported or provisional."),
    ("in_panel", "chip_candidates", "Whether the gene is in the CHIP panel."),
    ("curated_gene_status", "chip_gene_status.tsv",
     "retained: one of the 58 genes the paper keeps. removed_in_publication: one of the 16 it "
     "removes. added_after_publication: added to the authors' list since. not_curated: outside it."),
    ("gene_age_association_ukb", "supplementary table S6",
     "Association of mutations in this gene with age in UK Biobank at five or more supporting reads."),
    ("whitelist_status", "chip_curation",
     "whitelist, manual_review, not_whitelisted, not_on_curated_transcript or not_curated, following "
     "the authors' filter script."),
    ("whitelist_rule", "chip_curation", "The rule that produced the status."),
    ("in_curated_definition", "chip_curation",
     "yes when the change is whitelisted and the gene is retained or added since; review when the "
     "script sends it to manual review; no otherwise."),
    ("curation_note", "chip_curation", "Reason for the status, in words."),
    ("known_artefact", "supplementary table S2", "Matches a recurrent variant the paper removes as an artefact."),
    ("recurrent_inherited_variant", "supplementary table S4", "Matches a variant group the paper removes as inherited."),
    ("cohort_validated", "supplementary table S3", "Matches a variant associated with age or the TERT risk allele in the reference cohorts."),
    ("gene_specific_rule", "derived",
     "Outcome of a rule the paper sets for particular variants: ASXL1 p.G646Wfs*12 needs an allele "
     "fraction of 10 percent, ASXL1 p.G645Vfs*58 is removed, and TET2 or CBL domain missense is "
     "removed when its read counts fit an inherited heterozygote, three TET2 changes excepted."),
    ("meets_depth", "derived", f"Blood depth of {MIN_DEPTH} or more."),
    ("meets_alt_reads_5", "derived", f"{MIN_ALT_READS_EXOME} or more supporting reads in blood, the threshold the paper sets for UK Biobank exomes."),
    ("meets_alt_reads_3", "derived", f"{MIN_ALT_READS_GENOME} or more supporting reads in blood, the threshold the paper sets for All of Us genomes."),
    ("meets_allele_fraction", "derived", f"Blood allele fraction of {MIN_ALLELE_FRACTION:.0%} or more."),
    ("meets_orientation", "derived",
     "Variant fragments in blood seen in both orientations, from SAGE fragment strand bias, the "
     "counterpart of requiring F1R2 and F2R1 support."),
    ("clears_homopolymer_rule", "derived",
     f"No for an insertion or deletion in a single-base repeat of {HOMOPOLYMER_MIN_LENGTH} or more. "
     "ASXL1 p.G646Wfs*12 is exempt, as in the paper."),
    ("published_criteria_outcome_min_ad5", "derived",
     "retained when the call meets every criterion at five supporting reads, is inside the curated "
     "definition, is not constitutional and falls to no artefact, inherited-variant or gene-specific "
     "rule. Otherwise every criterion it misses."),
    ("published_criteria_outcome_min_ad3", "derived", "The same at three supporting reads."),
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
SPLIT = [c[0] for c in HEAD].index("tumour_to_blood_af_ratio")
COLUMNS = (HEAD[:SPLIT] + sample_columns("blood", "blood")
           + HEAD[SPLIT:] + sample_columns("tumour", "tumour") + TAIL)


def build(vcf, variants_tsv, outdir, case, blood_sample, tumour_sample, panel_dir=None):
    os.makedirs(outdir, exist_ok=True)
    with open(variants_tsv) as fh:
        called = list(csv.DictReader(fh, delimiter="\t"))
    order = {(r["chrom"], r["pos"], r["ref"], r["alt"]): i for i, r in enumerate(called)}
    by_key = {(r["chrom"], r["pos"], r["ref"], r["alt"]): r for r in called}
    curation = CuratedDefinition(panel_dir) if panel_dir and os.path.isdir(panel_dir) else None

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
                "in_panel", "curated_gene_status", "gene_age_association_ukb", "curated_transcript",
                "transcript_basis", "curated_hgvs_coding", "curated_hgvs_protein", "curated_effect",
                "curated_exon", "curated_codon", "whitelist_status", "whitelist_rule",
                "in_curated_definition", "curation_note", "known_artefact",
                "recurrent_inherited_variant", "cohort_validated",
                "detection_floor_at_this_depth", "present_in_tumour", "sage_filter", "warnings")}
            row["variant_type"] = variant_type(f[3], f[4])
            row["curated_accession"] = curation.genes.get(row["gene"], {}).get("accession", "") if curation else ""

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
            rule = row["whitelist_rule"]
            homopolymer = (is_indel and isinstance(reps, str) and len(reps) == 1
                           and str(repc).isdigit() and int(repc) >= HOMOPOLYMER_MIN_LENGTH
                           and rule != "asxl1_g646w")

            # rules the paper sets for particular variants
            specific, specific_miss = "", ""
            short = short_protein(row["curated_hgvs_protein"])
            if rule == "asxl1_g646w":
                ok = af != "" and af >= ASXL1_G646W_MIN_ALLELE_FRACTION
                specific = (f"ASXL1 p.G646Wfs*12 needs an allele fraction of "
                            f"{ASXL1_G646W_MIN_ALLELE_FRACTION:.0%}: {'met' if ok else 'not met'}")
                if not ok:
                    specific_miss = "ASXL1 p.G646Wfs*12 below 10% allele fraction"
            elif rule == "asxl1_g645v":
                specific = "ASXL1 p.G645Vfs*58 is removed as an artefact"
                specific_miss = specific
            elif (row["gene"] in BINOMIAL_GENES and short and row["whitelist_status"] == "whitelist"
                  and any(a <= int(short[1:-1]) <= b for a, b in MISSENSE_RANGES[row["gene"]])):
                p = row["p_if_inherited_heterozygote"]
                if (row["gene"], short) in BINOMIAL_EXEMPT:
                    specific = f"{row['gene']} {short} is exempt from the binomial test"
                elif p != "" and p >= BINOMIAL_P:
                    specific = (f"{row['gene']} domain missense with read counts that fit an inherited "
                                f"heterozygote (p = {p:.2g})")
                    specific_miss = "read counts fit an inherited heterozygote"
                else:
                    specific = f"{row['gene']} domain missense, read counts differ from an inherited heterozygote"
            row["gene_specific_rule"] = specific

            shared = [
                ("meets_depth", dp != "" and dp >= MIN_DEPTH, f"depth below {MIN_DEPTH}"),
                ("meets_allele_fraction", af != "" and af >= MIN_ALLELE_FRACTION,
                 f"allele fraction below {MIN_ALLELE_FRACTION:.0%}"),
                ("meets_orientation", both(row["blood_fragment_forward_fraction_alt"]) == "yes",
                 "variant fragments in one orientation only"),
                ("clears_homopolymer_rule", not homopolymer, "insertion or deletion in a homopolymer"),
            ]
            misses = []
            for name, ok, why in shared:
                row[name] = "yes" if ok else "no"
                if not ok:
                    misses.append(why)
            if row["in_curated_definition"] == "no":
                misses.append("outside the curated definition")
            elif row["in_curated_definition"] == "review":
                misses.append("awaiting manual review")
            if row["classification"] != "clone candidate":
                misses.append("classified as constitutional")
            if row["known_artefact"] == "yes":
                misses.append("matches a known artefact")
            if row["recurrent_inherited_variant"] == "yes":
                misses.append("matches a recurrent inherited variant")
            if specific_miss:
                misses.append(specific_miss)
            for n in (MIN_ALT_READS_EXOME, MIN_ALT_READS_GENOME):
                ok = alt != "" and alt >= n
                row[f"meets_alt_reads_{n}"] = "yes" if ok else "no"
                these = ([] if ok else [f"fewer than {n} supporting reads"]) + misses
                row[f"published_criteria_outcome_min_ad{n}"] = "retained" if not these else "; ".join(these)

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
    flat = [r for _, r in rows]
    kept = {n: [r for r in flat if r[f"published_criteria_outcome_min_ad{n}"] == "retained"]
            for n in (MIN_ALT_READS_EXOME, MIN_ALT_READS_GENOME)}
    inside = [r for r in flat if r["in_curated_definition"] == "yes"]
    summary = {
        "case": case,
        "variants": len(flat),
        "inside_the_curated_definition": len(inside),
        "awaiting_manual_review": sum(1 for r in flat if r["in_curated_definition"] == "review"),
        "retained_at_5_supporting_reads": len(kept[MIN_ALT_READS_EXOME]),
        "retained_at_3_supporting_reads": len(kept[MIN_ALT_READS_GENOME]),
        "retained_variants_at_3_supporting_reads": [
            f"{r['gene']} {r['curated_hgvs_protein'] or r['curated_hgvs_coding']}" for r in kept[MIN_ALT_READS_GENOME]],
        "inside_the_definition_and_below_threshold": [
            {"variant": f"{r['gene']} {r['curated_hgvs_protein'] if r['curated_hgvs_protein'] not in ('', 'p.?') else r['curated_hgvs_coding']}",
             "blood_alt_reads": r["blood_alt_reads"], "blood_depth": r["blood_depth"],
             "blood_allele_fraction": r["blood_allele_fraction"],
             "present_in_tumour": r["present_in_tumour"], "cohort_validated": r["cohort_validated"],
             "misses": r[f"published_criteria_outcome_min_ad{MIN_ALT_READS_GENOME}"]}
            for r in inside if r[f"published_criteria_outcome_min_ad{MIN_ALT_READS_GENOME}"] != "retained"],
        "inspect_sample_as_possible_poor_library": len(kept[MIN_ALT_READS_GENOME]) >= INSPECT_SAMPLE_AT,
        "inspection_basis": f"{INSPECT_SAMPLE_AT} or more retained variants in one sample, as the paper advises",
    }
    with open(os.path.join(outdir, f"{case}.chip_review_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1, default=str)
    return table, dictionary, flat, summary


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annotated-vcf", required=True)
    ap.add_argument("--variants-tsv", required=True, help="<case>.chip_variants.tsv from chip_candidates.py")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--case", required=True)
    ap.add_argument("--blood-sample", required=True)
    ap.add_argument("--tumour-sample", required=True)
    ap.add_argument("--panel-dir", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "panels"))
    a = ap.parse_args(argv)
    table, dictionary, rows, summary = build(a.annotated_vcf, a.variants_tsv, a.outdir, a.case,
                                             a.blood_sample, a.tumour_sample, a.panel_dir)
    print(f"variants {len(rows)}  columns {len(COLUMNS)}")
    print(json.dumps(summary, indent=1, default=str))
    print(table)
    print(dictionary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
