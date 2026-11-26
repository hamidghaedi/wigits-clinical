#!/usr/bin/env python3
"""Clonal haematopoiesis variant curation after Vlasschaert et al.

Applies the selection of Vlasschaert C, Mack T, Heimlich JB, et al., "A practical approach to curate
clonal hematopoiesis of indeterminate potential in human genetic data sets", Blood
2023;141(18):2214-2223, doi:10.1182/blood.2022018825, to variants called by SAGE and annotated by PAVE.

The whitelist is read from the authors' own files, copied unchanged into panels/vlasschaert/ from
https://github.com/briansha/Annovar_Whitelist_Filter_WDL. The rules below follow their
whitelist_filter_rscript_Aug26.R step by step, with the differences that a change of annotator forces
stated where they occur.

  Transcript. Their filter keeps a variant only when ANNOVAR annotates it on the RefSeq accession named
  for the gene. PAVE annotates every Ensembl transcript, so the entry whose RefSeq identifier equals
  that accession is used. Eleven accessions have no Ensembl counterpart in the Hartwig resources; for
  those genes the canonical transcript is used and the row says so. panels/chip_transcripts.tsv records
  the choice and how many whitelist residues agree with the transcript sequence.

  Exon numbers. ANNOVAR numbers ASXL1 and ASXL2 without the three-base third exon, so "exon 11-12" in
  their rule is exons 12 and 13 here. Verified against the exon labels of 30 ASXL1 and 3 ASXL2 variants
  in supplementary tables S2 and S3. PPM1D exons 5 and 6 carry the same numbers in both.

  ASXL1 and ASXL2 are in the loss-of-function gene file, so the R script whitelists their truncating
  variants in every exon and its exon test never removes one. The stated rule restricts them to the
  last two exons. A truncating variant outside those exons is sent to manual review here, with both
  facts in the note.

  Sixteen genes are removed in the paper under "Revisiting the CHIP gene list". Three genes, SRCAP,
  YLPM1 and ZBTB33, were added to the authors' files after publication. panels/chip_gene_status.tsv
  holds the status of every gene.

This module annotates. It removes nothing, consistent with docs/chip_reporting_proposal.md section 1.
"""
from __future__ import annotations

import csv
import os
import re

THREE_TO_ONE = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
                "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
                "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}

# Names used in the authors' loss-of-function file for genes their gene list calls KMT2A and KMT2D.
NORMALISE = {"MLL": "KMT2A", "MLL2": "KMT2D"}

# Missense admitted anywhere within these residue ranges (R script, step 4).
MISSENSE_RANGES = {"TET2": ((1104, 1481), (1843, 2002)), "CBL": ((381, 421),),
                   "CBLB": ((372, 412),), "ZBTB33": ((9, 126), (332, 591))}
# Exons, in Ensembl numbering, to which truncating and splice variants are restricted.
EXON_RULES = {"ASXL1": (12, 13), "ASXL2": (12, 13), "PPM1D": (5, 6)}
EXON_RULE_TEXT = {"ASXL1": "exons 11 and 12 of the published rule, 12 and 13 in Ensembl numbering",
                  "ASXL2": "exons 11 and 12 of the published rule, 12 and 13 in Ensembl numbering",
                  "PPM1D": "exons 5 and 6"}
# Variant classes the R script sends to manual review instead of whitelisting (step 5).
REVIEW_TRUNCATING_OR_SPLICE = {"GATA3": "truncating or splice variants in the zinc-finger domain",
                               "CSF3R": "truncating variants at c.741-791"}
REVIEW_FRAMESHIFT = {"NPM1": "frameshift p.W288fs"}
REVIEW_INFRAME_DELETION = {"CREBBP": "S1680del", "DNMT3A": "F732del and F752del",
                           "EP300": "VF1148_1149del", "KDM6A": "del419", "FLT3": "FY590-591GD and del835",
                           "JAK2": "exon 12 deletions and insertions", "KIT": "ins503, del560, del579, del551-559",
                           "MPL": "del513 and W515-518KT", "ZBTB33": "in-frame indels in the catalytic domains"}
REVIEW_INFRAME_INSERTION = {g: REVIEW_INFRAME_DELETION[g] for g in ("FLT3", "JAK2", "KIT", "MPL", "ZBTB33")}

# Gene-specific rules from the Results of the paper.
BINOMIAL_GENES = ("TET2", "CBL")                    # missense in the catalytic domains
BINOMIAL_EXEMPT = {("TET2", "H1904R"), ("TET2", "I1873T"), ("TET2", "T1884A")}
ASXL1_HOMOPOLYMER_CODONS = range(642, 647)          # the eight-G run encoding G643 to G646
ASXL1_G646W_MIN_ALLELE_FRACTION = 0.10


def _table(path):
    with open(path, newline="") as fh:
        lines = [l for l in fh.read().replace("\r", "").split("\n") if l and not l.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def short_protein(hgvs_protein):
    """p.Thr1932Pro -> T1932P. None when the change is not a single-residue substitution."""
    m = re.match(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$", hgvs_protein or "")
    if not m or m.group(1) not in THREE_TO_ONE or m.group(3) not in THREE_TO_ONE:
        return None
    return f"{THREE_TO_ONE[m.group(1)]}{m.group(2)}{THREE_TO_ONE[m.group(3)]}"


def first_residue(hgvs_protein):
    """p.Gly646fs -> ('G', 646). The first residue named, whatever follows."""
    m = re.match(r"p\.([A-Z][a-z]{2})(\d+)", hgvs_protein or "")
    if not m or m.group(1) not in THREE_TO_ONE:
        return None
    return THREE_TO_ONE[m.group(1)], int(m.group(2))


def parse_pave_ti(value):
    """PAVE_TI: Gene|GeneName|Transcript|Effects|SpliceRegion|HGVS.c|HGVS.p|RefSeqId|Exon|Codon."""
    out = []
    for entry in str(value or "").split(","):
        f = entry.split("|")
        if len(f) >= 10:
            out.append({"gene": f[1], "transcript": f[2], "effects": f[3], "hgvs_coding": f[5],
                        "hgvs_protein": f[6], "refseq": f[7], "exon": f[8], "codon": f[9]})
    return out


def reference_key(gene_accession, annovar):
    """Key for an ANNOVAR variant string from the supplement, e.g.
    DNMT3A:NM_022552:exon23:c.G2645A:p.R882H -> ('sub', 'R882H')."""
    parts = [p for p in re.split(r"[,;]", annovar) if p.strip()]
    part = next((p for p in parts if gene_accession and gene_accession in p), parts[0] if parts else "")
    m = re.search(r"p\.([A-Z])(\d+)(\S*)", part)
    if m:
        ref, pos, rest = m.group(1), m.group(2), m.group(3)
        if "fs" in rest:
            return ("fs", f"{ref}{pos}")
        if rest[:1] in ("X", "*"):
            return ("stop", f"{ref}{pos}")
        if re.fullmatch(r"[A-Z]", rest):
            return ("sub", f"{ref}{pos}{rest}")
        return ("other", f"{ref}{pos}{rest}")
    m = re.search(r"c\.(\d+[-+]\d+)([ACGT])>([ACGT])", part)
    if m:
        return ("splice", f"c.{m.group(1)}{m.group(2)}>{m.group(3)}")
    return ("other", part.strip())


class CuratedDefinition:
    def __init__(self, panel_dir):
        d = os.path.join(panel_dir, "vlasschaert")
        self.genes = {}
        for r in _table(os.path.join(d, "Full_CHIP_gene_list_08262022.txt")):
            self.genes[r["Gene"].strip()] = {"rule": " ".join(r["ReportedMut"].split()),
                                             "accession": r["Accession"].strip()}
        self.missense = {}
        for r in _table(os.path.join(d, "CHIP_missense_vars_cv_08262022.txt")):
            self.missense.setdefault(NORMALISE.get(r["Gene"].strip(), r["Gene"].strip()),
                                     set()).add(r["AAChange"].strip())
        self.aliases = {}
        alias_file = os.path.join(panel_dir, "chip_whitelist_aliases.tsv")
        if os.path.isfile(alias_file):
            for r in _table(alias_file):
                self.missense.setdefault(r["gene"], set()).add(r["aa_change"])
                self.aliases[(r["gene"], r["aa_change"])] = r["whitelist_entry"]
        self.lof_genes = {NORMALISE.get(r["Gene"].strip(), r["Gene"].strip())
                          for r in _table(os.path.join(d, "CHIP_nonsense_FS_vars_cv_08262022.txt"))}
        self.splice_genes = {NORMALISE.get(r["Gene"].strip(), r["Gene"].strip())
                             for r in _table(os.path.join(d, "CHIP_splice_vars_cv_08262022.txt"))}
        self.status, self.evidence = {}, {}
        for r in _table(os.path.join(panel_dir, "chip_gene_status.tsv")):
            self.status[r["gene"]] = r["status"]
            self.evidence[r["gene"]] = " ".join(x for x in (r.get("ukb_age_association", ""),
                                                            r.get("ukb_significance", "")) if x)
        self.transcripts = {r["gene"]: r for r in _table(os.path.join(panel_dir, "chip_transcripts.tsv"))}
        self.artefacts, self.germline, self.validated = {}, {}, {}
        for r in _table(os.path.join(panel_dir, "chip_artefact_variants.tsv")):
            store = self.germline if r["source"].startswith("S4") else self.artefacts
            store.setdefault(r["gene"], set()).add(self._ref_key(r["gene"], r["variant"]))
        for r in _table(os.path.join(panel_dir, "chip_validated_variants.tsv")):
            self.validated.setdefault(r["gene"], set()).add(self._ref_key(r["gene"], r["variant"]))

    def _ref_key(self, gene, annovar):
        return reference_key(self.genes.get(gene, {}).get("accession", ""), annovar)

    def gene_status(self, gene):
        """retained, removed_in_publication, added_after_publication or not_curated."""
        return self.status.get(gene, "not_curated") if gene in self.genes else "not_curated"

    def annotation(self, gene, pave_ti, canonical_transcript):
        """The PAVE annotation on the transcript the curated rule is written against.
        Returns (entry or None, basis)."""
        t = self.transcripts.get(gene)
        entries = [e for e in pave_ti if e["gene"] == gene]
        if t is None:
            hit = [e for e in entries if e["transcript"] == canonical_transcript]
            return (hit[0] if hit else None), "canonical"
        hit = [e for e in entries if e["transcript"] == t["transcript"]]
        return (hit[0] if hit else None), t["basis"]

    def evaluate(self, gene, ann, ref, alt):
        """Apply the whitelist to one annotation. Returns dict(status, rule, note) where status is
        whitelist, manual_review, not_whitelisted, not_on_curated_transcript or not_curated."""
        if gene not in self.genes:
            return {"status": "not_curated", "rule": "", "note": "gene is not in the curated list"}
        if ann is None:
            return {"status": "not_on_curated_transcript", "rule": "",
                    "note": "no annotation on the transcript the curated rule is written against"}
        effects = set(ann["effects"].split("&"))
        short = short_protein(ann["hgvs_protein"])
        try:
            exon = int(ann["exon"])
        except (TypeError, ValueError):
            exon = None
        stop, frameshift = "stop_gained" in effects, "frameshift_variant" in effects
        splice = bool(effects & {"splice_donor_variant", "splice_acceptor_variant"})
        inframe_del, inframe_ins = "inframe_deletion" in effects, "inframe_insertion" in effects
        missense = "missense_variant" in effects and short is not None

        def out(status, rule, note):
            return {"status": status, "rule": rule, "note": note}

        # 1) listed substitutions
        if missense and short in self.missense.get(gene, ()):
            alias = self.aliases.get((gene, short))
            return out("whitelist", "listed_missense",
                       f"{short} is on the whitelist" + (f" as {alias}" if alias else ""))
        # 4) missense anywhere in a named domain
        if missense and gene in MISSENSE_RANGES:
            pos = int(short[1:-1])
            if any(a <= pos <= b for a, b in MISSENSE_RANGES[gene]):
                return out("whitelist", "domain_missense",
                           f"missense at residue {pos}, inside the domain the rule names")
        # 2) and 4) truncating variants
        if stop or frameshift:
            kind = "nonsense" if stop else "frameshift"
            if gene == "ASXL1" and frameshift and abs(len(ref) - len(alt)) == 1:
                res = first_residue(ann["hgvs_protein"])
                if res and res[1] in ASXL1_HOMOPOLYMER_CODONS:
                    if len(alt) > len(ref):
                        return out("whitelist", "asxl1_g646w",
                                   "single-base insertion in the ASXL1 G643-G646 homopolymer, "
                                   "p.G646Wfs*12")
                    return out("not_whitelisted", "asxl1_g645v",
                               "single-base deletion in the ASXL1 G643-G646 homopolymer, "
                               "p.G645Vfs*58, which the paper removes as an artefact")
            if gene in EXON_RULES:
                if exon in EXON_RULES[gene]:
                    return out("whitelist", "truncating_in_exon",
                               f"{kind} in exon {exon}, within {EXON_RULE_TEXT[gene]}")
                if gene in self.lof_genes:
                    return out("manual_review", "truncating_outside_exon",
                               f"{kind} in exon {exon}; the published script admits it and the stated "
                               f"rule restricts {gene} to {EXON_RULE_TEXT[gene]}")
                return out("not_whitelisted", "truncating_outside_exon",
                           f"{kind} in exon {exon}, outside {EXON_RULE_TEXT[gene]}")
            if gene in self.lof_genes:
                return out("whitelist", "truncating", f"{kind} in a gene on the loss-of-function list")
        # 3) and 4) splice-site variants
        if splice:
            if gene in self.splice_genes:
                return out("whitelist", "splice", "splice-site variant in a gene on the splice list")
            if gene in ("ASXL1", "ASXL2") and exon in EXON_RULES[gene]:
                return out("whitelist", "splice_in_exon",
                           f"splice-site variant at exon {exon}, within {EXON_RULE_TEXT[gene]}")
        # 5) classes the script sends to manual review
        if (stop or frameshift or splice) and gene in REVIEW_TRUNCATING_OR_SPLICE:
            return out("manual_review", "review_truncating",
                       f"the rule admits {REVIEW_TRUNCATING_OR_SPLICE[gene]}; needs a reviewer")
        if frameshift and gene in REVIEW_FRAMESHIFT:
            return out("manual_review", "review_frameshift",
                       f"the rule admits {REVIEW_FRAMESHIFT[gene]}; needs a reviewer")
        if inframe_del and gene in REVIEW_INFRAME_DELETION:
            return out("manual_review", "review_inframe",
                       f"the rule admits {REVIEW_INFRAME_DELETION[gene]}; needs a reviewer")
        if inframe_ins and gene in REVIEW_INFRAME_INSERTION:
            return out("manual_review", "review_inframe",
                       f"the rule admits {REVIEW_INFRAME_INSERTION[gene]}; needs a reviewer")
        # not admitted
        if missense:
            same = sorted(x for x in self.missense.get(gene, ()) if x[:-1] == short[:-1])
            if same:
                return out("not_whitelisted", "", f"{short} is not on the whitelist; the same residue "
                                                  f"is listed as {', '.join(same)}")
            n = len(self.missense.get(gene, ()))
            return out("not_whitelisted", "", f"{short} is not among the {n} listed substitutions"
                       if n else "no missense change is admitted in this gene")
        if stop or frameshift:
            return out("not_whitelisted", "", "truncating variants are not admitted in this gene")
        if splice:
            return out("not_whitelisted", "", "splice-site variants are not admitted in this gene")
        label = sorted(effects)[0] if effects else "unannotated"
        return out("not_whitelisted", "", f"{label.replace('_', ' ')} is outside the definition")

    def variant_key(self, ann, ref, alt):
        if ann is None:
            return None
        effects = set(ann["effects"].split("&"))
        short = short_protein(ann["hgvs_protein"])
        if "missense_variant" in effects and short:
            return ("sub", short)
        res = first_residue(ann["hgvs_protein"])
        if "frameshift_variant" in effects and res:
            return ("fs", f"{res[0]}{res[1]}")
        if "stop_gained" in effects and res:
            return ("stop", f"{res[0]}{res[1]}")
        if effects & {"splice_donor_variant", "splice_acceptor_variant"}:
            return ("splice", ann["hgvs_coding"])
        return None

    def reference_cohorts(self, gene, ann, ref, alt):
        """(known artefact, recurrent inherited variant, validated) from tables S2, S4 and S3."""
        key = self.variant_key(ann, ref, alt)
        if key is None:
            return False, False, False
        return (key in self.artefacts.get(gene, ()), key in self.germline.get(gene, ()),
                key in self.validated.get(gene, ()))

    def in_definition(self, gene, evaluation):
        """yes, review or no, combining the whitelist with the status of the gene."""
        status = self.gene_status(gene)
        if status in ("not_curated", "removed_in_publication"):
            return "no"
        return {"whitelist": "yes", "manual_review": "review"}.get(evaluation["status"], "no")
