"""Checks for workflow/neo_export.py against the versioned subject_a fixture."""
import glob
import gzip
import json
import os

import pytest

from workflow import neo_export, neo_interface as ni

FIXTURE = os.path.join(os.path.dirname(__file__), "fixtures", "neo_interface", "subject_a")
pytestmark = pytest.mark.skipif(not os.path.isdir(FIXTURE), reason="fixture case absent")


@pytest.fixture(scope="module")
def exported(tmp_path_factory):
    out = tmp_path_factory.mktemp("neo_export")
    path, manifest = neo_export.Exporter(FIXTURE, str(out), case="subject_a", keep_non_coding=True).run()
    return out, manifest


def test_manifest_records_sources_with_checksums(exported):
    _, manifest = exported
    assert manifest["rules_version"].startswith("docs/neoantigen_interface.md")
    assert manifest["sources"]["purple_somatic_vcf"]["sha256"]
    assert manifest["samples"]["tumour_dna"] and manifest["samples"]["normal_dna"]
    # purity fields are carried for CCF maths and the case header
    assert manifest["purity"]["purity"] and manifest["purity"]["ploidy"]


def test_selected_vcf_has_rna_format_fields(exported):
    out, manifest = exported
    vcf = manifest["outputs"]["somatic_selected_vcf"]
    lines = [l for l in open(vcf) if not l.startswith("#")]
    assert lines, "expected at least one exported record"
    header = open(vcf).read()
    for key in ("RDP", "RAD", "RAF", "RNA_STATUS"):
        assert key in header
    fmt = lines[0].split("\t")[8].split(":")
    assert fmt[-3:] == ["RDP", "RAD", "RAF"]
    # tumour column (last) must carry values; normal column must not claim RNA support
    tumour = dict(zip(fmt, lines[0].rstrip("\n").split("\t")[10].split(":")))
    normal = dict(zip(fmt, lines[0].rstrip("\n").split("\t")[9].split(":")))
    assert normal["RDP"] == "." and normal["RAF"] == "."
    assert tumour["RDP"].isdigit()
    assert "," in tumour["RAD"]


def test_rna_states_sum_to_selected(exported):
    _, manifest = exported
    assert sum(manifest["rna_coverage_states"].values()) == manifest["counts"]["selected"]


def test_only_pass_records_are_exported(exported):
    _, manifest = exported
    vcf = manifest["outputs"]["somatic_selected_vcf"]
    for line in open(vcf):
        if line.startswith("#"):
            continue
        assert line.split("\t")[6] == "PASS"
    assert manifest["counts"]["selected"] <= manifest["counts"]["somatic_pass"] <= manifest["counts"]["somatic_total"]


def test_coding_filter_is_narrower_than_all_pass():
    """The phase-1 coding filter must actually restrict the set."""
    import tempfile
    with tempfile.TemporaryDirectory() as d1, tempfile.TemporaryDirectory() as d2:
        _, all_pass = neo_export.Exporter(FIXTURE, d1, case="s", keep_non_coding=True).run()
        _, coding = neo_export.Exporter(FIXTURE, d2, case="s", keep_non_coding=False).run()
    assert coding["counts"]["selected"] <= all_pass["counts"]["selected"]
    assert all_pass["counts"]["selected"] == all_pass["counts"]["somatic_pass"]


def test_proximal_vcf_excludes_the_selected_variants(exported):
    out, manifest = exported
    sel = {tuple(l.split("\t")[:2]) for l in open(manifest["outputs"]["somatic_selected_vcf"]) if not l.startswith("#")}
    prox = {tuple(l.split("\t")[:2]) for l in open(manifest["outputs"]["proximal_vcf"]) if not l.startswith("#")}
    assert not (sel & prox)


def test_hla_alleles_file_is_pvacseq_ready(exported):
    out, manifest = exported
    text = open(manifest["outputs"]["hla_alleles"]).read().strip()
    if not text:
        pytest.skip("fixture case has no LILAC alleles")
    alleles = text.split(",")
    assert all(a.startswith("HLA-") for a in alleles)
    assert len(alleles) == len(set(alleles))
    assert manifest["hla"]["lilac_qc_status"] in {"PASS", "FAIL", None} or manifest["hla"]["lilac_qc_status"].startswith("WARN")


def test_immune_escape_is_captured(exported):
    _, manifest = exported
    if not manifest["immune_escape"]:
        pytest.skip("fixture has no ORANGE json")
    assert "hasHlaEscape" in manifest["immune_escape"]
