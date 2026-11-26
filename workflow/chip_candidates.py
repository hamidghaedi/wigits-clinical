#!/usr/bin/env python3
"""Build the clonal haematopoiesis variant table for one case.

Reads the PAVE-annotated calls made over the blood alignment and writes every variant within the panel,
with the evidence behind it and a confidence band. Nothing is removed for being small or weakly
supported. A variant below the detection floor is reported with a warning, and the floor itself is
computed from the depth achieved in this sample rather than fixed in advance.

The reporting rule is set out in docs/chip_reporting_proposal.md section 1.

Outputs into --outdir:
  <case>.chip_variants.tsv        one row per variant, every variant
  <case>.chip_summary.json        depth, detection floor, panel scope, counts per band and class
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
from datetime import datetime, timezone

try:
    from chip_curation import CuratedDefinition, parse_pave_ti
except ImportError:  # the module sits beside this one in workflow/
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from chip_curation import CuratedDefinition, parse_pave_ti

RULES_VERSION = "docs/chip_reporting_proposal.md@2026-11-26; whitelist Annovar_Whitelist_Filter_WDL@8af2f67"

# A variant is treated as constitutional rather than as a clone at or above this allele fraction.
CONSTITUTIONAL_AF = 0.35
# Population frequency above which an inherited variant is the better explanation.
COMMON_IN_POPULATION = 0.001
# Supporting-read counts that define the bands, applied together with the computed floor.
CONFIDENT_ALT_READS = 5
SUPPORTED_ALT_READS = 3


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_info(field):
    out = {}
    for part in field.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
        elif part:
            out[part] = True
    return out


def parse_impact(value):
    """PAVE IMPACT: gene,transcript,consequence,codingEffect,spliceRegion,hgvsCoding,hgvsProtein,..."""
    if not value:
        return {}
    f = value.split(",")
    keys = ["gene", "transcript", "consequence", "coding_effect", "splice_region",
            "hgvs_coding", "hgvs_protein", "other_reportable", "worst_coding_effect", "gene_count"]
    return {k: (f[i] if i < len(f) else "") for i, k in enumerate(keys)}


def binomial_floor(depth, min_alt, power=0.95):
    """Smallest allele fraction at which min_alt or more reads are seen with the stated power."""
    if not depth:
        return None
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = (lo + hi) / 2
        # P(X >= min_alt) for X ~ Binomial(depth, mid)
        p = 1.0 - sum(math.comb(depth, k) * mid ** k * (1 - mid) ** (depth - k)
                      for k in range(min_alt))
        if p >= power:
            hi = mid
        else:
            lo = mid
    return hi


def classify(af, alt_reads, gnomad_af, sage_filter, floor_confident):
    """Return (classification, band, warnings). Nothing here removes a variant from the report."""
    warnings = []
    filters = [f for f in (sage_filter or "").split(",") if f and f != "PASS"]

    if gnomad_af is not None and gnomad_af > COMMON_IN_POPULATION:
        classification = "constitutional"
        warnings.append(f"present in the population at {gnomad_af:.4f}, inherited variant likely")
    elif af is not None and af >= CONSTITUTIONAL_AF:
        classification = "constitutional"
        warnings.append(f"allele fraction {af:.2f} is consistent with an inherited variant")
    else:
        classification = "clone candidate"

    if af is None or alt_reads is None:
        band = "provisional"
        warnings.append("allele fraction or supporting read count unavailable")
    elif af >= 0.05 and alt_reads >= CONFIDENT_ALT_READS:
        band = "confident"
    elif af >= 0.02 or alt_reads >= SUPPORTED_ALT_READS:
        band = "supported"
        if floor_confident and af < floor_confident:
            warnings.append(f"allele fraction {af:.3f} is below the {floor_confident:.3f} this depth "
                            f"resolves with 95 percent power")
    else:
        band = "provisional"
        warnings.append(f"allele fraction {af:.3f} with {alt_reads} supporting reads approaches the "
                        f"error floor of the assay")

    if "readStrandBias" in filters:
        warnings.append("supporting reads are imbalanced across strands")
        band = "provisional"
    if any(f.startswith("minMQF") or f == "mapQualRefAltDiff" for f in filters):
        warnings.append("mapping quality at this position is poor")
    if "minAvgRawBaseQual" in filters:
        warnings.append("base qualities of the supporting reads are low")
    if "dedupSnvMnv" in filters:
        warnings.append("duplicate of a neighbouring call")
    return classification, band, warnings


def build(vcf, outdir, case, blood_sample, tumour_sample, panel_file=None, mean_depth=None,
          panel_dir=None):
    os.makedirs(outdir, exist_ok=True)
    panel = []
    if panel_file and os.path.isfile(panel_file):
        for line in open(panel_file):
            if line.startswith("#") or line.startswith("gene\t"):
                continue
            g = line.split("\t")[0].strip()
            if g:
                panel.append(g)

    curation = CuratedDefinition(panel_dir) if panel_dir and os.path.isdir(panel_dir) else None

    samples, rows = [], []
    depths = []
    with _open(vcf) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            f = line.rstrip("\n").split("\t")
            info = parse_info(f[7])
            fmt = f[8].split(":")
            geno = {s: dict(zip(fmt, f[9 + i].split(":"))) for i, s in enumerate(samples)}
            blood = geno.get(blood_sample, {})
            tumour = geno.get(tumour_sample, {})

            def num(d, key, cast=float):
                try:
                    return cast(d.get(key))
                except (TypeError, ValueError):
                    return None

            dp = num(blood, "DP", int)
            af = num(blood, "AF")
            ad = blood.get("AD", "")
            alt_reads = None
            if "," in ad:
                try:
                    alt_reads = int(ad.split(",")[1])
                except (ValueError, IndexError):
                    alt_reads = None
            if dp:
                depths.append(dp)
            imp = parse_impact(info.get("IMPACT", ""))
            gnomad = None
            if info.get("GND_FREQ") not in (None, True, ""):
                try:
                    gnomad = float(info["GND_FREQ"])
                except ValueError:
                    gnomad = None
            floor_conf = binomial_floor(dp, CONFIDENT_ALT_READS) if dp else None
            sage_filter = info.get("SAGE_FILTER", f[6])
            classification, band, warnings = classify(af, alt_reads, gnomad, sage_filter, floor_conf)

            gene = imp.get("gene", "")
            cur = {k: "" for k in ("curated_gene_status", "gene_age_association_ukb", "curated_transcript",
                                   "transcript_basis", "curated_effect", "curated_hgvs_coding",
                                   "curated_hgvs_protein", "curated_exon", "curated_codon",
                                   "whitelist_status", "whitelist_rule", "in_curated_definition",
                                   "curation_note", "known_artefact", "recurrent_inherited_variant",
                                   "cohort_validated")}
            if curation is not None:
                ann, basis = curation.annotation(gene, parse_pave_ti(info.get("PAVE_TI", "")),
                                                 imp.get("transcript", ""))
                ev = curation.evaluate(gene, ann, f[3], f[4])
                art, germ, val = curation.reference_cohorts(gene, ann, f[3], f[4])
                cur.update({
                    "curated_gene_status": curation.gene_status(gene),
                    "gene_age_association_ukb": curation.evidence.get(gene, ""),
                    "curated_transcript": ann["transcript"] if ann else "",
                    "transcript_basis": basis if gene in curation.genes else "",
                    "curated_effect": ann["effects"] if ann else "",
                    "curated_hgvs_coding": ann["hgvs_coding"] if ann else "",
                    "curated_hgvs_protein": ann["hgvs_protein"] if ann else "",
                    "curated_exon": ann["exon"] if ann else "",
                    "curated_codon": ann["codon"] if ann else "",
                    "whitelist_status": ev["status"], "whitelist_rule": ev["rule"],
                    "in_curated_definition": curation.in_definition(gene, ev),
                    "curation_note": ev["note"],
                    "known_artefact": "yes" if art else "no",
                    "recurrent_inherited_variant": "yes" if germ else "no",
                    "cohort_validated": "yes" if val else "no"})
                status = cur["curated_gene_status"]
                if status == "removed_in_publication":
                    cur["curation_note"] += "; the gene is one of the 16 the paper removes"
                elif status == "added_after_publication":
                    cur["curation_note"] += "; the gene was added to the authors' list after publication"
                if cur["in_curated_definition"] == "no":
                    warnings.append(f"outside the curated definition: {cur['curation_note']}")
                elif cur["in_curated_definition"] == "review":
                    warnings.append(f"needs manual review: {cur['curation_note']}")
                if art:
                    warnings.append("matches a recurrent artefact reported in the reference cohorts")
                if germ:
                    warnings.append("matches a variant the reference cohorts found to be inherited")
                if val:
                    warnings.append("associated with age or the TERT risk allele in the reference "
                                    "cohorts, which supports a genuine clone")

            if classification == "clone candidate" and (num(tumour, "AF") or 0) > 0:
                warnings.append("also called in the tumour, so tumour material in the blood draw and "
                                "a clone infiltrating the tumour are both open explanations")

            rows.append({
                "case": case,
                "gene": imp.get("gene", ""),
                "chrom": f[0], "pos": f[1], "ref": f[3], "alt": f[4],
                "consequence": imp.get("consequence", ""),
                "coding_effect": imp.get("coding_effect", ""),
                "hgvs_coding": imp.get("hgvs_coding", ""),
                "hgvs_protein": imp.get("hgvs_protein", ""),
                "transcript": imp.get("transcript", ""),
                "blood_depth": dp if dp is not None else "",
                "blood_alt_reads": alt_reads if alt_reads is not None else "",
                "blood_allele_fraction": f"{af:.4f}" if af is not None else "",
                "tumour_depth": tumour.get("DP", ""),
                "tumour_allele_fraction": tumour.get("AF", ""),
                "present_in_tumour": ("yes" if (num(tumour, "AF") or 0) > 0 else "no"),
                "gnomad_frequency": f"{gnomad:.6f}" if gnomad is not None else "",
                "clinvar": info.get("CLNSIG", ""),
                "mappability": info.get("MAPPABILITY", ""),
                "trinucleotide": info.get("TNC", ""),
                "repeat_count": info.get("RC_REPC", ""),
                "repeat_sequence": info.get("RC_REPS", ""),
                "sage_filter": sage_filter,
                "detection_floor_at_this_depth": f"{floor_conf:.4f}" if floor_conf else "",
                "classification": classification,
                "confidence_band": band,
                "in_panel": "yes" if imp.get("gene", "") in panel else "no",
                **cur,
                "warnings": "; ".join(warnings),
            })

    rows.sort(key=lambda r: (r["classification"] != "clone candidate",
                             {"confident": 0, "supported": 1, "provisional": 2}[r["confidence_band"]],
                             -(float(r["blood_allele_fraction"]) if r["blood_allele_fraction"] else 0)))

    tsv = os.path.join(outdir, f"{case}.chip_variants.tsv")
    if rows:
        with open(tsv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(rows)

    median_depth = sorted(depths)[len(depths) // 2] if depths else None
    depth_for_floor = int(mean_depth) if mean_depth else median_depth
    summary = {
        "case": case,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "rules_version": RULES_VERSION,
        "reporting_rule": ("Every variant called within the panel is reported. Variants below the "
                           "detection floor carry a warning and are not removed."),
        "blood_sample": blood_sample,
        "tumour_sample": tumour_sample,
        "panel_genes": len(panel),
        "panel_genes_with_a_call": len({r["gene"] for r in rows if r["in_panel"] == "yes"}),
        "panel_genes_without_a_call": sorted(set(panel) - {r["gene"] for r in rows}),
        "variants_reported": len(rows),
        "median_depth_at_called_positions": median_depth,
        "mean_depth_of_sample": mean_depth,
        "detection_floor": {
            "basis": f"binomial power at depth {depth_for_floor}",
            "allele_fraction_for_5_reads_at_95pc_power":
                round(binomial_floor(depth_for_floor, CONFIDENT_ALT_READS), 4) if depth_for_floor else None,
            "allele_fraction_for_3_reads_at_95pc_power":
                round(binomial_floor(depth_for_floor, SUPPORTED_ALT_READS), 4) if depth_for_floor else None,
        },
        "curated_definition": None if curation is None else {
            "source": "panels/vlasschaert, Annovar_Whitelist_Filter_WDL at 8af2f67",
            "genes": len(curation.genes),
            "retained_genes": sum(1 for s in curation.status.values() if s == "retained"),
            "calls_by_gene_status": {s: sum(1 for r in rows if r["curated_gene_status"] == s)
                                     for s in ("retained", "removed_in_publication",
                                               "added_after_publication", "not_curated")},
            "calls_by_whitelist_status": {s: sum(1 for r in rows if r["whitelist_status"] == s)
                                          for s in ("whitelist", "manual_review", "not_whitelisted",
                                                    "not_on_curated_transcript", "not_curated")},
            "calls_inside": sum(1 for r in rows if r["in_curated_definition"] == "yes"),
            "calls_for_manual_review": sum(1 for r in rows if r["in_curated_definition"] == "review"),
            "matching_a_known_artefact": sum(1 for r in rows if r["known_artefact"] == "yes"),
            "matching_a_recurrent_inherited_variant":
                sum(1 for r in rows if r["recurrent_inherited_variant"] == "yes"),
            "validated_in_the_reference_cohorts":
                sum(1 for r in rows if r["cohort_validated"] == "yes"),
        },
        "classification_counts": {k: sum(1 for r in rows if r["classification"] == k)
                                   for k in ("clone candidate", "constitutional")},
        "confidence_counts": {k: sum(1 for r in rows if r["confidence_band"] == k)
                               for k in ("confident", "supported", "provisional")},
        "clone_candidates_by_band": {
            k: sum(1 for r in rows if r["classification"] == "clone candidate"
                   and r["confidence_band"] == k) for k in ("confident", "supported", "provisional")},
        "inputs": {"annotated_vcf": os.path.abspath(vcf), "panel": panel_file},
        "outputs": {"variants": tsv},
    }
    path = os.path.join(outdir, f"{case}.chip_summary.json")
    with open(path, "w") as fh:
        json.dump(summary, fh, indent=1)
    return summary


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annotated-vcf", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--case", required=True)
    ap.add_argument("--blood-sample", required=True)
    ap.add_argument("--tumour-sample", required=True)
    ap.add_argument("--panel-file", default=None)
    ap.add_argument("--panel-dir", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "panels"),
                    help="directory holding vlasschaert/ and the chip_*.tsv tables")
    ap.add_argument("--mean-depth", type=float, default=None,
                    help="mean depth of the blood sample, used for the reported detection floor")
    a = ap.parse_args(argv)
    s = build(a.annotated_vcf, a.outdir, a.case, a.blood_sample, a.tumour_sample,
              a.panel_file, a.mean_depth, a.panel_dir)
    print(json.dumps({k: s[k] for k in ("variants_reported", "classification_counts",
                                        "confidence_counts", "clone_candidates_by_band",
                                        "curated_definition",
                                        "detection_floor", "panel_genes",
                                        "panel_genes_with_a_call")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
