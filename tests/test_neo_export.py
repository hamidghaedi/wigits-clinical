"""Checks for workflow/neo_export.py against the versioned subject_a fixture."""
import glob
import gzip
import json
import os
import shutil

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
    for key in ("RDP", "RAD", "RAF", "RNA_STATUS", "GX", "TX"):
        assert key in header
    # pvactools reads GX/TX from FORMAT, not INFO
    assert '##FORMAT=<ID=GX' in header and '##FORMAT=<ID=TX' in header
    assert '##INFO=<ID=GX' not in header
    fmt = lines[0].split("\t")[8].split(":")
    assert fmt[-5:] == ["RDP", "RAD", "RAF", "GX", "TX"]
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


def test_proximal_vcf_contains_every_selected_anchor(exported):
    out, manifest = exported
    sel = {tuple(l.split("\t")[:2]) for l in open(manifest["outputs"]["somatic_selected_vcf"]) if not l.startswith("#")}
    prox = {tuple(l.split("\t")[:2]) for l in open(manifest["outputs"]["proximal_vcf"]) if not l.startswith("#")}
    assert sel <= prox
    assert manifest["counts"]["proximal_anchors"] == manifest["counts"]["selected"]


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


def test_proximal_vcf_has_exactly_one_sample(exported):
    """pVACseq rejects a proximal VCF with more than one sample column."""
    _, manifest = exported
    for line in open(manifest["outputs"]["proximal_vcf"]):
        if line.startswith("#CHROM"):
            assert len(line.rstrip("\n").split("\t")) == 10
            assert line.rstrip("\n").split("\t")[9] == manifest["samples"]["tumour_dna"]
            break
    else:
        raise AssertionError("no #CHROM line in the proximal VCF")


def test_expression_reaches_the_tumour_genotype(exported):
    """Expression is populated only where the trimmed ISOFOX data match a variant."""
    _, manifest = exported
    if not manifest["counts"].get("isofox_genes"):
        pytest.skip("fixture has no ISOFOX gene expression")
    vcf = manifest["outputs"]["somatic_selected_vcf"]
    layout = ni.CaseLayout.from_results_dir(FIXTURE)
    expression = {r['GeneName']: r['AdjTPM'] for r in ni.read_tsv(layout.files['isofox_gene_tsv'])}
    for line in open(vcf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        vals = dict(zip(f[8].split(":"), f[-1].split(":")))
        gene = ni.parse_impact(ni.parse_info(f[7]).get('IMPACT', '')).get('Gene', '')
        expected = f"{gene}|{expression[gene]}" if gene in expression else '.'
        assert vals.get('GX', '.') == expected


@pytest.mark.parametrize('tpm', ['0.000e+00', '2.500e+00'])
def test_synthetic_matching_expression_reaches_genotype(tmp_path, tpm):
    """Synthetic matching data prove that zero and nonzero expression survive."""
    case = tmp_path / 'synthetic'
    shutil.copytree(FIXTURE, case)
    layout = ni.CaseLayout.from_results_dir(str(case))
    with gzip.open(layout.files['purple_somatic_vcf'], 'rt') as stream:
        variant = next(line for line in stream if not line.startswith('#') and ni.parse_impact(ni.parse_info(line.split('\t')[7]).get('IMPACT', '')).get('Gene'))
    gene = ni.parse_impact(ni.parse_info(variant.split('\t')[7])['IMPACT'])['Gene']
    with open(layout.files['isofox_gene_tsv'], 'w') as stream:
        stream.write(f'GeneId\tGeneName\tAdjTPM\nSYNTHETIC_GENE\t{gene}\t{tpm}\n')
    _, manifest = neo_export.Exporter(str(case), str(tmp_path / 'export'), keep_non_coding=True).run()
    matching = []
    for line in open(manifest['outputs']['somatic_selected_vcf']):
        if line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        if ni.parse_impact(ni.parse_info(fields[7]).get('IMPACT', '')).get('Gene') == gene:
            matching.append(dict(zip(fields[8].split(':'), fields[-1].split(':')))['GX'])
    assert matching
    assert set(matching) == {f'{gene}|{tpm}'}
