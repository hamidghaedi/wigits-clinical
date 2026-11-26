#!/usr/bin/env python3
"""Unified neoantigen candidate table for one case.

Joins the pVACseq aggregated report with Hartwig NEO's own peptide scores, LILAC allele status,
SAGE_APPEND RNA support, PURPLE purity/copy number and the ORANGE immune-escape flags, then assigns an
interpretable tier. Field provenance is docs/neoantigen_interface.md; the tiering rules are stated in
docs/neoantigen_extension_proposal.md section 3 item 4 and repeated in the `--explain` output.

Eligibility is kept separate from ranking: a candidate whose variant or HLA basis is unsound is excluded
or held regardless of how good its binding prediction is.

Outputs (into --outdir):
  <case>.neo_candidates.tsv          one row per (variant, allele, peptide), tiered, with reasons
  <case>.neo_concordance.tsv         pVACseq vs Hartwig NEO on peptides both scored
  <case>.neo_candidates.summary.json counts per tier, per RNA state, concordance summary, provenance
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import neo_interface as ni  # noqa: E402

RULES_VERSION = "docs/neoantigen_extension_proposal.md@2026-11-24"

# --- thresholds (named, so a decision never rests on an invisible number) -----------------------
BINDING_STRONG_NM = 500.0        # pVACseq default strong-binder IC50
BINDING_PERCENTILE = 2.0         # and percentile rank
CLONAL_SUBCL_MAX = 0.25          # SUBCL above this is treated as subclonal
ALLELE_CN_MIN = 0.5              # LILAC tumour copy number below this = allele effectively lost
FFPE_TRANSITIONS = {("C", "T"), ("G", "A")}   # deamination signature
FFPE_MUC_MIN = 3                 # minimum unique fragment coordinates for a C>T/G>A call
# Minimum ratio of wild-type to mutant IC50 required before a candidate is prioritised. A ratio of 1.0
# means the two peptides bind equally well, so the mutation adds nothing a T cell can distinguish, and a
# ratio below 1.0 means the wild-type peptide is the better binder. Published work does not settle on a
# single cut-off, so 2.0 is used here and recorded with every candidate as fold_change.
MIN_FOLD_CHANGE = 2.0

TIERS = ("Prioritise", "Review", "Hold", "Exclude")


def _f(value, default=None):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def normalise_variant_key(text: str):
    """pVACseq ID 'chr13-27793619-27793620-G-T' -> ('chr13', '27793620', 'G', 'T').

    pVACseq writes a half-open start as well as the position; the second coordinate is the VCF POS.
    """
    parts = text.split("-")
    if len(parts) == 5:
        chrom, _start, pos, ref, alt = parts
        return (chrom, pos, ref, alt)
    if len(parts) == 4:
        chrom, pos, ref, alt = parts
        return (chrom, pos, ref, alt)
    return None


class CandidateBuilder:
    def __init__(self, results_case_dir, pvacseq_dir, outdir, case, selected_vcf=None):
        self.layout = ni.CaseLayout.from_results_dir(results_case_dir)
        self.pvacseq_dir = pvacseq_dir
        self.outdir = outdir
        self.case = case
        # RNA_STATUS / RDP / RAD / RAF live on the exported selected VCF, not on PURPLE's somatic VCF:
        # they are written by neo_export.py from SAGE_APPEND. Reading them from PURPLE yields None for
        # every record and silently sinks every candidate out of the expressed tiers.
        self.selected_vcf = selected_vcf or os.path.join(outdir, f"{case}.somatic.selected.vcf")
        os.makedirs(outdir, exist_ok=True)

    # -- inputs ---------------------------------------------------------------------------
    def aggregated_rows(self):
        hits = [f for f in os.listdir(self.pvacseq_dir) if f.endswith("aggregated.tsv")]
        if not hits:
            raise FileNotFoundError(f"no *aggregated.tsv in {self.pvacseq_dir}")
        path = os.path.join(self.pvacseq_dir, sorted(hits)[0])
        self.aggregated_path = path
        return ni.read_tsv(path)

    def all_epitope_index(self):
        """peptide -> best (lowest) IC50 per allele, from the full epitope table."""
        hits = [f for f in os.listdir(self.pvacseq_dir) if f.endswith("all_epitopes.tsv")]
        if not hits:
            return {}
        path = os.path.join(self.pvacseq_dir, sorted(hits)[0])
        self.all_epitopes_path = path
        index = {}
        with open(path) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                pep = row.get("MT Epitope Seq") or row.get("Epitope Seq")
                allele = row.get("HLA Allele")
                if not pep or not allele:
                    continue
                key = (pep, ni.normalise_allele(allele))
                ic50 = _f(row.get("Median MT IC50 Score") or row.get("Best MT IC50 Score"))
                if ic50 is None:
                    continue
                if key not in index or ic50 < index[key]:
                    index[key] = ic50
        return index

    def neo_peptides(self):
        """(peptide, allele) -> Hartwig NEO scores; plus variant-level NEO facts."""
        path = self.layout.files.get("neo_peptide_scores_tsv")
        peptides, variants = {}, {}
        if path:
            for row in ni.read_tsv(path):
                allele = ni.normalise_allele(row["Allele"])
                peptides[(row["Peptide"], allele)] = {
                    "neo_likelihood": row.get("Likelihood"),
                    "neo_likelihood_rank": row.get("LikelihoodRank"),
                    "neo_exp_likelihood": row.get("ExpLikelihood"),
                    "neo_effective_tpm": row.get("EffectiveTpm"),
                    "neo_recog_sim": row.get("RecogSim"),
                    "neo_other_allele_recog_sim": row.get("OtherAlleleRecogSim"),
                    "neo_allele_cn": row.get("AlleleCN"),
                    "neo_allele_disrupted": row.get("AlleleDisrupted"),
                }
        ndata = self.layout.files.get("neo_data_tsv")
        if ndata:
            for row in ni.read_tsv(ndata):
                info = row.get("VariantInfo", "")
                if ":" in info and ";" not in info:
                    chrom, pos, ref, alt = (info.split(":") + ["", "", "", ""])[:4]
                    variants[(chrom, pos, ref, alt)] = {
                        "neo_variant_type": row.get("VariantType"),
                        "neo_variant_copy_number": row.get("VariantCopyNumber"),
                        "neo_subclonal_likelihood": row.get("SubclonalLikelihood"),
                        "neo_nmd_min": row.get("NmdMin"),
                    }
        return peptides, variants

    def rna_index(self):
        """variant key -> RNA support written by neo_export.py (three-state, never a bare zero)."""
        out = {}
        if not os.path.isfile(self.selected_vcf):
            return out
        samples = ni.vcf_samples(self.selected_vcf)
        t_idx = samples.index(self.layout.tumour_dna)
        with ni._open(self.selected_vcf) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                info = ni.parse_info(f[7])
                vals = ni.sample_values(f[8], f[9 + t_idx])
                rad = ni.int_list(vals.get("RAD"), default=[0, 0])
                out[(f[0], f[1], f[3], f[4])] = {
                    "rna_status": info.get("RNA_STATUS", "NOT_COVERED"),
                    "rna_depth_sage": vals.get("RDP", ""),
                    "rna_alt_reads": rad[1] if len(rad) > 1 else 0,
                    "rna_vaf_sage": vals.get("RAF", ""),
                    "gene_tpm": (vals.get("GX", "") or "").split("|")[-1] if vals.get("GX", ".") != "." else "",
                    "transcript_tpm": (vals.get("TX", "") or "").split("|")[-1] if vals.get("TX", ".") != "." else "",
                }
        return out

    def somatic_index(self):
        """variant key -> SAGE/PURPLE facts needed for eligibility and FFPE screening."""
        out = {}
        path = self.layout.files["purple_somatic_vcf"]
        samples = ni.vcf_samples(path)
        t_idx = samples.index(self.layout.tumour_dna)
        n_idx = samples.index(self.layout.normal_dna)
        with ni._open(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if f[6] != "PASS":
                    continue
                info = ni.parse_info(f[7])
                tum = ni.sample_values(f[8], f[9 + t_idx])
                norm = ni.sample_values(f[8], f[9 + n_idx])
                n_ad = ni.int_list(norm.get("AD"), default=[0, 0])
                out[(f[0], f[1], f[3], f[4])] = {
                    "tier_sage": info.get("TIER"),
                    "hotspot": "HOTSPOT" in info,
                    "purple_vcn": info.get("PURPLE_VCN"),
                    "purple_cn": info.get("PURPLE_CN"),
                    "purple_af": info.get("PURPLE_AF"),
                    "subclonal_likelihood": info.get("SUBCL", "0"),
                    "biallelic": "BIALLELIC" in info,
                    "mappability": info.get("MAPPABILITY"),
                    "trinucleotide": info.get("TNC"),
                    "unique_fragments": tum.get("MUC"),
                    "normal_alt_reads": n_ad[1] if len(n_ad) > 1 else 0,
                    "ref": f[3],
                    "alt": f[4],
                }
        return out

    def lilac_index(self):
        out = {}
        path = self.layout.files.get("lilac_tsv")
        if not path:
            return out
        for row in ni.read_tsv(path):
            allele = ni.normalise_allele(row["Allele"])
            somatic_hits = sum(_f(row.get(k), 0.0) or 0.0 for k in
                               ("SomaticMissense", "SomaticNonsenseOrFrameshift", "SomaticSplice", "SomaticInframeIndel"))
            prev = out.get(allele)
            cn = _f(row.get("TumorCopyNumber"))
            if prev is None or (cn is not None and cn < (prev["lilac_tumour_cn"] or cn)):
                out[allele] = {"lilac_tumour_cn": cn, "lilac_somatic_hits": somatic_hits,
                               "lilac_rna_total": row.get("RnaTotal")}
        return out

    def genes_of_interest(self, path):
        if not path or not os.path.isfile(path):
            return set()
        with open(path) as fh:
            return {line.strip() for line in fh if line.strip()}

    # -- tiering --------------------------------------------------------------------------
    def classify(self, rec) -> tuple[str, list[str]]:
        reasons = []

        # --- eligibility (never rescued by a good prediction) ---
        if rec["normal_alt_reads"] and int(rec["normal_alt_reads"]) > 0 and not rec["hotspot"]:
            reasons.append("alt reads in matched normal")
            return "Exclude", reasons
        if rec["allele_lost"]:
            reasons.append(f"restricting allele lost in tumour (LILAC CN {rec['lilac_tumour_cn']})")
            return "Exclude", reasons
        if rec["lilac_qc_status"] == "FAIL":
            reasons.append("LILAC typing QC FAIL")
            return "Exclude", reasons

        # --- holds: evidence missing or technically suspect ---
        if rec["ffpe_suspect"]:
            reasons.append(f"possible FFPE deamination ({rec['trinucleotide']}, MUC {rec['unique_fragments']})")
            return "Hold", reasons
        if rec["rna_status"] == "NOT_COVERED":
            reasons.append("RNA not covered at this position (absence not established)")
            return "Hold", reasons
        if rec["allele_disrupted"]:
            reasons.append("restricting allele somatically disrupted")
            return "Hold", reasons

        # --- ranking ---
        strong = rec["ic50_mt"] is not None and rec["ic50_mt"] <= BINDING_STRONG_NM
        strong_pct = rec["percentile_mt"] is not None and rec["percentile_mt"] <= BINDING_PERCENTILE
        clonal = rec["subclonal"] is False
        expressed = rec["rna_status"] == "SUPPORTED"
        fc = rec["fold_change"]
        differential = fc is None or fc >= MIN_FOLD_CHANGE

        if strong:
            reasons.append(f"IC50 {rec['ic50_mt']:.0f} nM")
        if strong_pct:
            reasons.append(f"percentile {rec['percentile_mt']:.2f}")
        if expressed:
            reasons.append("mutant allele expressed in RNA")
        elif rec["rna_status"] == "COVERED_NO_ALT":
            reasons.append("RNA covered, mutant allele not seen")
        if clonal:
            reasons.append("clonal")
        else:
            reasons.append(f"subclonal (SUBCL {rec['subclonal_likelihood']})")
        if rec["neo_likelihood"]:
            reasons.append(f"NEO likelihood {rec['neo_likelihood']}")
        if rec["gene_of_interest"]:
            reasons.append("cancer gene census gene")
        if fc is not None:
            reasons.append(f"wild-type binds {fc:.1f}x less well"
                            if fc >= 1 else f"wild-type binds {1 / fc:.1f}x better")

        if strong and strong_pct and expressed and clonal and differential:
            return "Prioritise", reasons
        if strong or strong_pct:
            return "Review", reasons
        reasons.append("no strong binding prediction")
        return "Hold", reasons

    # -- main -----------------------------------------------------------------------------
    def run(self, genes_file=None):
        agg = self.aggregated_rows()
        ic50_index = self.all_epitope_index()
        neo_pep, neo_var = self.neo_peptides()
        somatic = self.somatic_index()
        rna = self.rna_index()
        lilac = self.lilac_index()
        goi = self.genes_of_interest(genes_file)
        qc = (ni.read_tsv(self.layout.files["lilac_qc_tsv"])[0].get("Status")
              if self.layout.files.get("lilac_qc_tsv") else None)
        escape = ni.immune_escape(self.layout.files["orange_json"]) if self.layout.files.get("orange_json") else {}
        purity = ni.purity_fields(self.layout.files["purple_purity_tsv"])

        rows, concordance = [], []
        tier_counts = {t: 0 for t in TIERS}
        for a in agg:
            key = normalise_variant_key(a["ID"])
            som = somatic.get(key, {}) if key else {}
            rsup = rna.get(key, {}) if key else {}
            allele = ni.normalise_allele(a["Allele"]) if a.get("Allele") else ""
            la = lilac.get(allele, {})
            peptide = a.get("Best Peptide", "")
            npep = neo_pep.get((peptide, allele), {})
            nvar = neo_var.get(key, {}) if key else {}
            subcl = _f(som.get("subclonal_likelihood"), 0.0) or 0.0
            ref, alt = som.get("ref", ""), som.get("alt", "")
            muc = som.get("unique_fragments")
            rec = {
                "case": self.case,
                "variant": a["ID"],
                "chrom": key[0] if key else "",
                "pos": key[1] if key else "",
                "ref": ref,
                "alt": alt,
                "gene": a.get("Gene", ""),
                "aa_change": a.get("AA Change", ""),
                "transcript": a.get("Best Transcript", ""),
                "tsl": a.get("TSL", ""),
                "peptide": peptide,
                "allele": allele,
                "allele_class": ni.allele_class(allele) if allele else "",
                "ic50_mt": _f(a.get("IC50 MT")),
                "ic50_wt": _f(a.get("IC50 WT")),
                "percentile_mt": _f(a.get("%ile MT")),
                "fold_change": None,  # filled below once both affinities are known
                "percentile_wt": _f(a.get("%ile WT")),
                "pvacseq_tier": a.get("Tier", ""),
                "num_passing_peptides": a.get("Num Passing Peptides", ""),
                "dna_vaf": a.get("DNA VAF", ""),
                "rna_vaf": a.get("RNA VAF", ""),
                "rna_depth": a.get("RNA Depth", ""),
                "rna_expression": a.get("RNA Expr", ""),
                "gene_tpm": rsup.get("gene_tpm", ""),
                "transcript_tpm": rsup.get("transcript_tpm", ""),
                "rna_depth_sage": rsup.get("rna_depth_sage", ""),
                "rna_alt_reads": rsup.get("rna_alt_reads", ""),
                "rna_vaf_sage": rsup.get("rna_vaf_sage", ""),
                "rna_status": rsup.get("rna_status", "NOT_COVERED"),
                "sage_tier": som.get("tier_sage", ""),
                "hotspot": som.get("hotspot", False),
                "purple_vcn": som.get("purple_vcn", ""),
                "purple_cn": som.get("purple_cn", ""),
                "subclonal_likelihood": som.get("subclonal_likelihood", ""),
                "subclonal": subcl > CLONAL_SUBCL_MAX,
                "biallelic": som.get("biallelic", False),
                "mappability": som.get("mappability", ""),
                "trinucleotide": som.get("trinucleotide", ""),
                "unique_fragments": muc,
                "normal_alt_reads": som.get("normal_alt_reads", 0),
                "lilac_tumour_cn": la.get("lilac_tumour_cn"),
                "lilac_somatic_hits": la.get("lilac_somatic_hits", 0),
                "lilac_rna_total": la.get("lilac_rna_total", ""),
                "lilac_qc_status": qc,
                "gene_of_interest": a.get("Gene", "") in goi,
                "case_immune_escape": ";".join(k for k, v in escape.items() if v) or "none",
                "purity": purity.get("purity"),
                "ploidy": purity.get("ploidy"),
                **{k: npep.get(k, "") for k in (
                    "neo_likelihood", "neo_likelihood_rank", "neo_exp_likelihood", "neo_effective_tpm",
                    "neo_recog_sim", "neo_other_allele_recog_sim", "neo_allele_cn", "neo_allele_disrupted")},
                **{k: nvar.get(k, "") for k in (
                    "neo_variant_type", "neo_variant_copy_number", "neo_subclonal_likelihood", "neo_nmd_min")},
            }
            if rec["ic50_mt"] and rec["ic50_wt"]:
                rec["fold_change"] = rec["ic50_wt"] / rec["ic50_mt"]
            rec["allele_lost"] = (rec["lilac_tumour_cn"] is not None
                                  and rec["lilac_tumour_cn"] < ALLELE_CN_MIN)
            rec["allele_disrupted"] = (str(rec.get("neo_allele_disrupted", "")).lower() == "true"
                                       or (rec["lilac_somatic_hits"] or 0) > 0)
            rec["ffpe_suspect"] = bool(
                (ref, alt) in FFPE_TRANSITIONS and muc not in (None, "", ".")
                and int(muc) < FFPE_MUC_MIN)
            tier, reasons = self.classify(rec)
            rec["tier"] = tier
            rec["tier_reasons"] = "; ".join(reasons)
            tier_counts[tier] += 1
            rows.append(rec)

            if npep:
                concordance.append({
                    "case": self.case, "variant": a["ID"], "gene": rec["gene"], "peptide": peptide,
                    "allele": allele, "pvacseq_ic50_mt": rec["ic50_mt"],
                    "pvacseq_percentile_mt": rec["percentile_mt"],
                    "neo_likelihood": npep.get("neo_likelihood"),
                    "neo_likelihood_rank": npep.get("neo_likelihood_rank"),
                    "neo_exp_likelihood": npep.get("neo_exp_likelihood"),
                    "neo_recog_sim": npep.get("neo_recog_sim"),
                    "agreement": ("both_strong" if (rec["ic50_mt"] or 1e9) <= BINDING_STRONG_NM
                                  and (_f(npep.get("neo_likelihood"), 0) or 0) > 0 else "differs"),
                })

        # peptides NEO scored that pVACseq also produced anywhere (not only the best per variant)
        shared = sum(1 for k in neo_pep if k in ic50_index)

        cand_path = os.path.join(self.outdir, f"{self.case}.neo_candidates.tsv")
        if rows:
            with open(cand_path, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
                w.writeheader()
                w.writerows(rows)
        conc_path = os.path.join(self.outdir, f"{self.case}.neo_concordance.tsv")
        if concordance:
            with open(conc_path, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(concordance[0].keys()), delimiter="\t")
                w.writeheader()
                w.writerows(concordance)

        summary = {
            "case": self.case,
            "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "rules_version": RULES_VERSION,
            "thresholds": {"binding_strong_nm": BINDING_STRONG_NM, "binding_percentile": BINDING_PERCENTILE,
                            "clonal_subcl_max": CLONAL_SUBCL_MAX, "allele_cn_min": ALLELE_CN_MIN,
                            "ffpe_muc_min": FFPE_MUC_MIN, "min_fold_change": MIN_FOLD_CHANGE},
            "candidates": len(rows),
            "tiers": tier_counts,
            "rna_states": {s: sum(1 for r in rows if r["rna_status"] == s)
                            for s in ("SUPPORTED", "COVERED_NO_ALT", "NOT_COVERED")},
            "pvacseq_tiers": {t: sum(1 for r in rows if r["pvacseq_tier"] == t)
                               for t in sorted({r["pvacseq_tier"] for r in rows})},
            "neo_peptide_rows": len(neo_pep),
            "peptides_scored_by_both": shared,
            "best_peptides_matched_in_neo": len(concordance),
            "case_immune_escape": {k: v for k, v in escape.items() if v} or "none",
            "purity": {k: purity.get(k) for k in ("purity", "ploidy", "msStatus", "tmbPerMb", "tml")},
            "lilac_qc_status": qc,
            "inputs": {"aggregated": getattr(self, "aggregated_path", None),
                        "all_epitopes": getattr(self, "all_epitopes_path", None),
                        **{k: v for k, v in self.layout.files.items() if v}},
            "outputs": {"candidates": cand_path, "concordance": conc_path if concordance else None},
        }
        sum_path = os.path.join(self.outdir, f"{self.case}.neo_candidates.summary.json")
        with open(sum_path, "w") as fh:
            json.dump(summary, fh, indent=1)
        return summary


EXPLANATION = """Tier rules (eligibility first, then ranking):

  Exclude   alt reads in the matched normal (unless a known hotspot); restricting HLA allele lost in the
            tumour (LILAC tumour copy number < %(cn).1f); LILAC typing QC FAIL.
  Hold      possible FFPE deamination (C>T/G>A with fewer than %(muc)d unique fragment coordinates);
            RNA not covered at the position (absence of the mutant allele is NOT established);
            restricting allele somatically disrupted; or no strong binding prediction.
  Review    a strong binding prediction (IC50 <= %(nm).0f nM or percentile <= %(pct).1f) but missing one of
            expression, clonality or mutant-over-wild-type differential.
  Prioritise strong IC50 AND strong percentile AND mutant allele expressed in RNA AND clonal
            (SUBCL <= %(subcl).2f) AND the wild-type peptide binds at least %(fc).1f times less well.

Every row carries tier_reasons, so a tier can always be read back to the evidence that produced it.
""" % {"cn": ALLELE_CN_MIN, "muc": FFPE_MUC_MIN, "nm": BINDING_STRONG_NM,
        "pct": BINDING_PERCENTILE, "subcl": CLONAL_SUBCL_MAX, "fc": MIN_FOLD_CHANGE}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    # none of these are argparse-required, so that --explain can be run on its own to print the
    # tier rules; the combination is validated below instead
    ap.add_argument("--results-case-dir")
    ap.add_argument("--pvacseq-dir", help="pvacseq/MHC_Class_I directory")
    ap.add_argument("--outdir")
    ap.add_argument("--case")
    ap.add_argument("--genes-of-interest-file", default=None)
    ap.add_argument("--selected-vcf", default=None,
                    help="exported <case>.somatic.selected.vcf (default: <outdir>/<case>.somatic.selected.vcf)")
    ap.add_argument("--explain", action="store_true", help="print the tier rules and exit")
    args = ap.parse_args(argv)
    if args.explain:
        print(EXPLANATION)
        return 0
    missing = [n for n in ("results_case_dir", "pvacseq_dir", "outdir", "case")
               if not getattr(args, n)]
    if missing:
        ap.error("the following arguments are required: "
                 + ", ".join("--" + m.replace("_", "-") for m in missing))
    summary = CandidateBuilder(args.results_case_dir, args.pvacseq_dir, args.outdir, args.case,
                               args.selected_vcf).run(args.genes_of_interest_file)
    print(json.dumps({k: summary[k] for k in
                      ("candidates", "tiers", "rna_states", "pvacseq_tiers", "peptides_scored_by_both",
                       "best_peptides_matched_in_neo", "case_immune_escape", "lilac_qc_status")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
