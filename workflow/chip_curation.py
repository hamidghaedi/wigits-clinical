#!/usr/bin/env python3
"""Variant-level clonal haematopoiesis curation, following Vlasschaert et al.

Implements the curated target definition from the supplementary tables of Vlasschaert C, Mack T,
Heimlich JB, et al., "A practical approach to curate clonal hematopoiesis of indeterminate potential in
human genetic datasets", Blood 2023;141(18):2214-2223, doi:10.1182/blood.2022018825.

The resource files are built from that supplement and live in workflow/panels/:

  chip_genes_vlasschaert.tsv   74 candidate genes, the mutation classes that count in each, and
                               whether the gene is on the retained list
  chip_artefact_variants.tsv   recurrent artefacts and variants failing their binomial test
  chip_validated_variants.tsv  variants associated with age and the TERT risk allele rs7705526

Four points of interpretation, stated because they change results:

  The supplement lists 74 candidate genes in table S1 and gives the evidence for each in table S6,
  whose column "Is on M-CHIP list" is true for 55 and false for 19. The 55 are taken here as the
  retained definition. A variant in one of the other 19 is reported as a candidate gene that was not
  retained, with the variant rule still evaluated so the reader can see what it would have been.

  A gene on the list does not accept every variant in it. Forty-seven of the 74 accept truncating
  variants, and 46 carry an explicit list of missense positions. CEBPA and NF1, for instance, accept
  truncating variants only, so a missense call in either is outside the definition.

  Synonymous variants are outside the definition for every gene.

  Membership of the validated list is positive evidence rather than a requirement. A variant absent
  from it is not thereby an artefact, it is simply unsupported by the cohort evidence.

This module annotates. It does not remove anything, consistent with the reporting rule in
docs/chip_reporting_proposal.md section 1.
"""
from __future__ import annotations

import csv
import os
import re

THREE_TO_ONE = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
                "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
                "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}

TRUNCATING = {"stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant"}


def _rows(path):
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


class CuratedDefinition:
    def __init__(self, genes_tsv, artefacts_tsv=None, validated_tsv=None):
        self.genes = {}
        for r in _rows(genes_tsv):
            self.genes[r["gene"]] = {
                "rule": r["reported_mutation_rule"],
                "accession": r["accession"],
                "admits_truncating": r["admits_truncating"] == "true",
                "missense": set(x for x in r["curated_missense"].split(",") if x),
                "retained": r.get("on_mchip_list", "true") == "true",
                "age_association": " ".join(x for x in (r.get("ukb_age_association_min_ad5", ""),
                                                        r.get("ukb_significance_min_ad5", "")) if x),
            }
        self.artefacts = {}
        if artefacts_tsv and os.path.isfile(artefacts_tsv):
            for r in _rows(artefacts_tsv):
                self.artefacts.setdefault(r["gene"], []).append(r["variant"])
        self.validated = {}
        if validated_tsv and os.path.isfile(validated_tsv):
            for r in _rows(validated_tsv):
                self.validated.setdefault(r["gene"], []).append(r["variant"])

    @staticmethod
    def short_protein(hgvs_protein):
        """p.Thr1932Pro -> T1932P. Returns None when the change is not a simple substitution."""
        m = re.match(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$", hgvs_protein or "")
        if not m:
            return None
        a, pos, b = m.groups()
        if a not in THREE_TO_ONE or b not in THREE_TO_ONE:
            return None
        return f"{THREE_TO_ONE[a]}{pos}{THREE_TO_ONE[b]}"

    def gene_status(self, gene):
        """retained, candidate_not_retained, or not_curated."""
        entry = self.genes.get(gene)
        if entry is None:
            return "not_curated"
        return "retained" if entry["retained"] else "candidate_not_retained"

    def classify(self, gene, consequence, hgvs_protein):
        """Return (in_definition: bool, reason: str)."""
        ok, reason = self._classify_variant(gene, consequence, hgvs_protein)
        entry = self.genes.get(gene)
        if entry is not None and not entry["retained"]:
            verdict = "the variant rule admits this change" if ok else reason
            return False, f"candidate gene that the reference study did not retain; {verdict}"
        return ok, reason

    def _classify_variant(self, gene, consequence, hgvs_protein):
        entry = self.genes.get(gene)
        if entry is None:
            return False, "gene is not in the curated list"
        if consequence == "synonymous_variant":
            return False, "synonymous variants are outside the definition"
        if consequence in TRUNCATING:
            if entry["admits_truncating"]:
                return True, "truncating variant in a gene whose rule admits them"
            return False, "this gene admits only the listed missense positions"
        if consequence == "start_lost":
            return False, "start-lost is not among the admitted classes for any gene"
        if consequence == "missense_variant":
            short = self.short_protein(hgvs_protein)
            if short and short in entry["missense"]:
                return True, f"{short} is on the curated missense list for this gene"
            if not entry["missense"]:
                return False, "this gene admits only truncating variants"
            return False, f"missense not among the {len(entry['missense'])} curated positions"
        return False, f"consequence {consequence} is outside the definition"

    def is_known_artefact(self, gene, hgvs_protein):
        short = self.short_protein(hgvs_protein)
        hits = self.artefacts.get(gene, [])
        if not hits:
            return False
        needle = (short or hgvs_protein or "").lower()
        return bool(needle) and any(needle in h.lower() for h in hits)

    def is_validated(self, gene, hgvs_protein):
        short = self.short_protein(hgvs_protein)
        hits = self.validated.get(gene, [])
        if not hits:
            return False
        needle = (short or hgvs_protein or "").lower()
        return bool(needle) and any(needle in h.lower() for h in hits)


def sample_level_artefact_signature(variants, min_count=5, vaf_window=0.02):
    """Flag the pattern in supplementary table S5: several truncating calls in one sample at a
    similar allele fraction, which indicates a poor-quality library rather than several clones."""
    trunc = [v for v in variants if v.get("consequence") in TRUNCATING]
    if len(trunc) < min_count:
        return False, ""
    fracs = sorted(float(v["blood_allele_fraction"]) for v in trunc
                   if v.get("blood_allele_fraction"))
    if len(fracs) < min_count:
        return False, ""
    for i in range(len(fracs) - min_count + 1):
        window = fracs[i:i + min_count]
        if window[-1] - window[0] <= vaf_window:
            return True, (f"{min_count} or more truncating calls within {vaf_window:.0%} allele "
                          f"fraction of one another, the low-quality sample signature")
    return False, ""
