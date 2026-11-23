"""Interface checks for docs/neoantigen_interface.md.

Pure-function checks always run. Checks that need real oncoanalyser outputs run against
fixtures under tests/fixtures/neo_interface/ (cut from subject_a) and are skipped when absent.
"""
import json
import os

import pytest

from workflow import neo_interface as ni

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures", "neo_interface")
CASE_DIR = os.path.join(FIXTURES, "subject_a")


# --- allele naming -------------------------------------------------------------------------


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("C1203", "HLA-C*12:03"),
        ("C*12:03", "HLA-C*12:03"),
        ("HLA-C*12:03", "HLA-C*12:03"),
        ("A0201", "HLA-A*02:01"),
        ("B0702", "HLA-B*07:02"),
        ("DPB1*02:01", "HLA-DPB1*02:01"),
        ("DQA1*01:02", "HLA-DQA1*01:02"),
    ],
)
def test_normalise_allele(raw, expected):
    assert ni.normalise_allele(raw) == expected


def test_allele_round_trip():
    for neo in ("A0201", "B3801", "C0702"):
        assert ni.neo_allele(ni.normalise_allele(neo)) == neo


def test_allele_class():
    assert ni.allele_class("HLA-A*02:01") == "I"
    assert ni.allele_class("DPB1*02:01") == "II"


def test_rejects_nonsense_allele():
    with pytest.raises(ValueError):
        ni.normalise_allele("not-an-allele")


# --- RNA coverage three-state rule ------------------------------------------------------------


@pytest.mark.parametrize(
    "depth,alt,expected",
    [
        (50, 12, "SUPPORTED"),
        (10, 2, "SUPPORTED"),
        (50, 0, "COVERED_NO_ALT"),
        (10, 1, "COVERED_NO_ALT"),
        (9, 5, "NOT_COVERED"),
        (0, 0, "NOT_COVERED"),
        (None, None, "NOT_COVERED"),
    ],
)
def test_rna_coverage_status(depth, alt, expected):
    assert ni.rna_coverage_status(depth, alt) == expected


def test_no_alt_at_low_depth_is_not_absence():
    """The distinction the interface spec insists on."""
    assert ni.rna_coverage_status(4, 0) == "NOT_COVERED"
    assert ni.rna_coverage_status(40, 0) == "COVERED_NO_ALT"


# --- PAVE parsing --------------------------------------------------------------------------


def test_parse_impact_from_real_record():
    impact = ("CHD5,ENST00000262450,missense_variant,MISSENSE,false,c.4312C>T,"
              "p.Arg1438Cys,,MISSENSE,1")
    parsed = ni.parse_impact(impact)
    assert parsed["Gene"] == "CHD5"
    assert parsed["Transcript"] == "ENST00000262450"
    assert parsed["CanonicalCodingEffect"] == "MISSENSE"
    assert parsed["HgvsProtein"] == "p.Arg1438Cys"
    assert parsed["GenesAffected"] == "1"


def test_parse_pave_ti_multi_transcript():
    value = ("ENSG00000116254|CHD5|ENST00000262450|missense_variant|false|c.4312C>T|p.Arg1438Cys|NM_015557|29|1438,"
             "ENSG00000116254|CHD5|ENST00000462991|missense_variant|false|c.2458C>T|p.Arg820Cys||18|820")
    rows = ni.parse_pave_ti(value)
    assert [r["Transcript"] for r in rows] == ["ENST00000262450", "ENST00000462991"]
    assert rows[0]["RefSeqId"] == "NM_015557"
    assert rows[1]["RefSeqId"] == ""


def test_strip_transcript_version():
    assert ni.strip_transcript_version("ENST00000262450.11") == "ENST00000262450"
    assert ni.strip_transcript_version("ENST00000262450") == "ENST00000262450"


# --- VCF helpers ---------------------------------------------------------------------------


def test_parse_info_and_sample_values():
    info = ni.parse_info("PURPLE_VCN=1.04;SUBCL=0.031;HOTSPOT;TIER=HIGH_CONFIDENCE")
    assert info["PURPLE_VCN"] == "1.04"
    assert info["HOTSPOT"] is True
    vals = ni.sample_values("GT:AD:DP:AF:RC_CNT", "0/1:178,12:190:0.063:12,0,0,1,178,191")
    assert ni.int_list(vals["AD"]) == [178, 12]
    assert ni.int_list(vals["RC_CNT"])[0] == 12
    assert int(vals["DP"]) == 190


# --- fixture-backed checks -------------------------------------------------------------------

pytestmark_case = pytest.mark.skipif(
    not os.path.isdir(CASE_DIR), reason="no oncoanalyser fixture case; run tools/make_neo_fixtures.sh"
)


@pytestmark_case
def test_case_layout_resolves_samples_and_files():
    layout = ni.CaseLayout.from_results_dir(CASE_DIR)
    assert layout.tumour_dna and layout.normal_dna
    assert layout.tumour_dna != layout.normal_dna
    assert layout.files["purple_somatic_vcf"].endswith(".purple.somatic.vcf.gz")
    assert os.path.isfile(layout.files["purple_purity_tsv"])


@pytestmark_case
def test_sage_append_carries_the_rna_sample_last():
    layout = ni.CaseLayout.from_results_dir(CASE_DIR)
    if not layout.files["sage_append_vcf"]:
        pytest.skip("case has no sage_append output")
    samples = ni.vcf_samples(layout.files["sage_append_vcf"])
    assert samples[:2] == [layout.normal_dna, layout.tumour_dna]
    assert layout.tumour_rna == samples[-1]


@pytestmark_case
def test_lilac_alleles_are_deduplicated_and_prefixed():
    layout = ni.CaseLayout.from_results_dir(CASE_DIR)
    if not layout.files["lilac_tsv"]:
        pytest.skip("case has no LILAC output")
    class_i = ni.lilac_alleles(layout.files["lilac_tsv"], hla_class="I")
    assert class_i, "expected at least one class I allele"
    assert len(class_i) == len(set(class_i))
    assert all(a.startswith("HLA-") for a in class_i)


@pytestmark_case
def test_neo_peptide_alleles_match_lilac_class_i():
    layout = ni.CaseLayout.from_results_dir(CASE_DIR)
    if not (layout.files["neo_peptide_scores_tsv"] and layout.files["lilac_tsv"]):
        pytest.skip("case lacks NEO scorer or LILAC output")
    neo_rows = ni.read_tsv(layout.files["neo_peptide_scores_tsv"])
    neo_set = {ni.normalise_allele(r["Allele"]) for r in neo_rows}
    lilac_i = set(ni.lilac_alleles(layout.files["lilac_tsv"], hla_class="I"))
    assert neo_set <= lilac_i, f"NEO alleles not in LILAC class I: {sorted(neo_set - lilac_i)}"


@pytestmark_case
def test_purity_tsv_has_no_qcstatus_column():
    """Documented correction: QC status lives in *.purple.qc, not the purity TSV."""
    layout = ni.CaseLayout.from_results_dir(CASE_DIR)
    fields = ni.purity_fields(layout.files["purple_purity_tsv"])
    assert "qcStatus" not in fields
    assert {"purity", "ploidy", "msStatus", "tmbPerMb"} <= set(fields)


@pytestmark_case
def test_immune_escape_keys_present():
    layout = ni.CaseLayout.from_results_dir(CASE_DIR)
    if not layout.files["orange_json"]:
        pytest.skip("case has no ORANGE json")
    escape = ni.immune_escape(layout.files["orange_json"])
    assert {"hasHlaEscape", "hasAntigenPresentationPathwayEscape"} <= set(escape)
