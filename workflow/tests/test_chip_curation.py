#!/usr/bin/env python3
"""Tests for the clonal haematopoiesis curation rules.

Each case states a variant as PAVE would annotate it on the curated transcript and the verdict the
whitelist of Vlasschaert et al. gives it. Run with pytest, or directly with python.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from chip_curation import CuratedDefinition, reference_key  # noqa: E402

CD = CuratedDefinition(os.path.join(HERE, "..", "panels"))


def ann(effects, protein="", coding="", exon="1", codon="1"):
    return {"gene": "", "transcript": "", "effects": effects, "hgvs_coding": coding,
            "hgvs_protein": protein, "refseq": "", "exon": exon, "codon": codon}


def verdict(gene, a, ref="C", alt="T"):
    e = CD.evaluate(gene, a, ref, alt)
    return e["status"], e["rule"], CD.in_definition(gene, e)


def test_hotspot_substitutions_are_whitelisted():
    for gene, change in (("DNMT3A", "p.Arg882His"), ("DNMT3A", "p.Arg882Cys"), ("JAK2", "p.Val617Phe"),
                         ("SF3B1", "p.Lys700Glu"), ("SRSF2", "p.Pro95His"), ("IDH2", "p.Arg140Gln"),
                         ("U2AF1", "p.Ser34Phe"), ("TP53", "p.Arg175His")):
        assert verdict(gene, ann("missense_variant", change)) == ("whitelist", "listed_missense", "yes"), (gene, change)


def test_a_different_change_at_a_listed_residue_is_outside():
    e = CD.evaluate("DNMT3A", ann("missense_variant", "p.Ile898Asn"), "A", "T")
    assert e["status"] == "not_whitelisted" and "I898T" in e["note"]


def test_synonymous_is_outside_everywhere():
    assert verdict("RUNX1", ann("synonymous_variant", "p.Ser424="))[2] == "no"


def test_domain_missense():
    assert verdict("TET2", ann("missense_variant", "p.His1881Arg"))[:2] == ("whitelist", "domain_missense")
    assert verdict("TET2", ann("missense_variant", "p.Pro1723Ser"))[0] == "not_whitelisted"
    assert verdict("CBL", ann("missense_variant", "p.Cys404Tyr"))[:2] == ("whitelist", "domain_missense")
    assert verdict("CBL", ann("missense_variant", "p.Ala100Thr"))[0] == "not_whitelisted"


def test_truncating_follows_the_gene_lists():
    assert verdict("TET2", ann("stop_gained", "p.Gln1034*"))[:2] == ("whitelist", "truncating")
    assert verdict("DNMT3A", ann("frameshift_variant", "p.Phe732fs"), "CT", "C")[0] == "whitelist"
    assert verdict("JAK2", ann("stop_gained", "p.Gln100*"))[0] == "not_whitelisted"
    assert verdict("CEBPA", ann("missense_variant", "p.Phe33Cys"))[0] == "not_whitelisted"


def test_joined_effects_are_read():
    assert verdict("DNMT3A", ann("splice_donor_variant&intron_variant", "p.?", "c.1429+1G>A", "11"))[:2] == ("whitelist", "splice")


def test_splice_is_admitted_only_in_splice_genes():
    assert verdict("NF1", ann("splice_acceptor_variant&intron_variant", "p.?", "c.100-1G>A"))[0] == "not_whitelisted"
    assert verdict("TET2", ann("splice_donor_variant", "p.?", "c.3409+1G>T"))[0] == "whitelist"


def test_exon_restrictions():
    assert verdict("ASXL1", ann("stop_gained", "p.Arg693*", exon="13"))[:2] == ("whitelist", "truncating_in_exon")
    assert verdict("ASXL1", ann("stop_gained", "p.Arg404*", exon="12"))[0] == "whitelist"
    assert verdict("ASXL1", ann("stop_gained", "p.Gln50*", exon="2")) == ("manual_review", "truncating_outside_exon", "review")
    assert verdict("PPM1D", ann("stop_gained", "p.Arg552*", exon="6"))[0] == "whitelist"
    assert verdict("PPM1D", ann("stop_gained", "p.Arg100*", exon="1"))[0] == "not_whitelisted"


def test_asxl1_homopolymer():
    ins = CD.evaluate("ASXL1", ann("frameshift_variant", "p.Gly646fs", exon="13"), "A", "AG")
    assert (ins["status"], ins["rule"]) == ("whitelist", "asxl1_g646w")
    dele = CD.evaluate("ASXL1", ann("frameshift_variant", "p.Gly645fs", exon="13"), "AG", "A")
    assert (dele["status"], dele["rule"]) == ("not_whitelisted", "asxl1_g645v")
    other = CD.evaluate("ASXL1", ann("frameshift_variant", "p.Glu635fs", exon="13"), "A" * 24, "A")
    assert (other["status"], other["rule"]) == ("whitelist", "truncating_in_exon")


def test_manual_review_classes():
    assert verdict("NPM1", ann("frameshift_variant", "p.Trp288fs"), "C", "CTCTG")[0] == "manual_review"
    assert verdict("FLT3", ann("inframe_insertion", "p.Glu598_Tyr599insAsp"), "C", "CATG")[0] == "manual_review"
    assert verdict("JAK2", ann("inframe_deletion", "p.Asn542_Glu543del"), "CAATGA", "C")[0] == "manual_review"
    assert verdict("CSF3R", ann("stop_gained", "p.Gln776*"))[0] == "manual_review"


def test_gene_status():
    assert CD.gene_status("DNMT3A") == "retained"
    assert CD.gene_status("CBLB") == "retained" and CD.gene_status("SETDB1") == "retained"
    assert CD.gene_status("KMT2D") == "removed_in_publication"
    assert CD.gene_status("ZBTB33") == "added_after_publication"
    assert CD.gene_status("ATM") == "not_curated"
    assert sum(1 for s in CD.status.values() if s == "retained") == 58
    assert sum(1 for s in CD.status.values() if s == "removed_in_publication") == 16
    # admitted by the variant rule, and still outside because the gene was removed
    assert verdict("KMT2D", ann("stop_gained", "p.Leu804*")) == ("whitelist", "truncating", "no")
    assert verdict("ATM", ann("stop_gained", "p.Ser421*")) == ("not_curated", "", "no")


def test_variant_without_the_curated_transcript_is_dropped():
    assert CD.evaluate("RUNX1", None, "C", "T")["status"] == "not_on_curated_transcript"


def test_etnk1_alias_on_current_numbering():
    assert verdict("ETNK1", ann("missense_variant", "p.Asn155Ser"))[0] == "whitelist"


def test_reference_cohort_matching():
    assert reference_key("NM_022552", "DNMT3A:NM_022552:exon23:c.G2645A:p.R882H") == ("sub", "R882H")
    assert reference_key("NM_015338", "ASXL1:NM_015338:exon12:c.1927dupG:p.G646Wfs*12") == ("fs", "G646")
    assert reference_key("NM_015338", "ASXL1:NM_015338:exon12:c.T4368G:p.Y1456X") == ("stop", "Y1456")
    assert reference_key("NM_004380", "NM_001079846:exon26:c.4281-2A>T;NM_004380:exon27:c.4395-2A>T") == ("splice", "c.4395-2A>T")
    art, germ, val = CD.reference_cohorts("DNMT3A", ann("missense_variant", "p.Arg882His"), "C", "T")
    assert (art, germ, val) == (False, False, True)
    assert CD.reference_cohorts("ASXL1", ann("frameshift_variant", "p.Gly646fs"), "A", "AG")[2] is True
    assert CD.reference_cohorts("ASXL1", ann("stop_gained", "p.Tyr1456*"), "A", "C")[0] is True
    assert CD.reference_cohorts("TP53", ann("missense_variant", "p.Arg110Cys"), "G", "A")[1] is True
    # table S2 lists the frameshift p.A1944Gfs*22; the substitution p.A1944G must not match it
    assert CD.reference_cohorts("CREBBP", ann("missense_variant", "p.Ala1944Gly"), "G", "C")[0] is False
    assert CD.reference_cohorts("CREBBP", ann("frameshift_variant", "p.Ala1944fs"), "C", "CG")[0] is True


if __name__ == "__main__":
    failed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn(); print("ok   ", name)
            except AssertionError as exc:
                failed += 1; print("FAIL ", name, exc)
    raise SystemExit(1 if failed else 0)
