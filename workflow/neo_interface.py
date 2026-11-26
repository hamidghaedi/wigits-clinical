"""Shared helpers for the neoantigen layer's interface with oncoanalyser outputs.

Every rule implemented here is specified in docs/neoantigen_interface.md. Keep the two in step.
Read-only: nothing in this module writes to a run directory.
"""
from __future__ import annotations

import csv
import glob
import gzip
import json
import os
import re
from dataclasses import dataclass, field

# --- allele naming (interface spec section 6) --------------------------------------------------

_NEO_ALLELE = re.compile(r"^(?P<gene>[A-Z]+[0-9]*)(?P<group>\d{2,3})(?P<protein>\d{2})$")
_HLA_ALLELE = re.compile(r"^(?:HLA-)?(?P<gene>[A-Z]+[0-9]*)\*(?P<group>\d{2,3}):(?P<protein>\d{2})$")
CLASS_I_GENES = ("A", "B", "C")
CLASS_II_GENES = ("DPA1", "DPB1", "DQA1", "DQB1", "DRA", "DRB1", "DRB3", "DRB4", "DRB5")


def normalise_allele(allele: str, *, prefix: bool = True) -> str:
    """`C1203`, `C*12:03`, `HLA-C*12:03` -> `HLA-C*12:03` (or without prefix)."""
    allele = allele.strip()
    m = _HLA_ALLELE.match(allele) or _NEO_ALLELE.match(allele)
    if not m:
        raise ValueError(f"unrecognised HLA allele: {allele!r}")
    core = f"{m['gene']}*{m['group']}:{m['protein']}"
    return f"HLA-{core}" if prefix else core


def neo_allele(allele: str) -> str:
    """`HLA-C*12:03` -> `C1203` (the form NEO's peptide_scores.tsv uses)."""
    core = normalise_allele(allele, prefix=False)
    gene, rest = core.split("*")
    return gene + rest.replace(":", "")


def allele_class(allele: str) -> str:
    gene = normalise_allele(allele, prefix=False).split("*")[0]
    if gene in CLASS_I_GENES:
        return "I"
    if gene in CLASS_II_GENES:
        return "II"
    return "other"


# --- RNA coverage state (interface spec section 2) ---------------------------------------------

RNA_MIN_DEPTH = 10
RNA_MIN_ALT = 2


def rna_coverage_status(depth, alt_count, *, min_depth: int = RNA_MIN_DEPTH, min_alt: int = RNA_MIN_ALT) -> str:
    """Three states; an absent mutant allele at low depth is NOT evidence of non-expression."""
    if depth is None:
        return "NOT_COVERED"
    if depth < min_depth:
        return "NOT_COVERED"
    if (alt_count or 0) >= min_alt:
        return "SUPPORTED"
    return "COVERED_NO_ALT"


# --- VCF access --------------------------------------------------------------------------------


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def vcf_samples(path) -> list[str]:
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
    raise ValueError(f"no #CHROM line in {path}")


def parse_info(info: str) -> dict:
    out = {}
    for item in info.split(";"):
        if not item:
            continue
        k, _, v = item.partition("=")
        out[k] = v if _ else True
    return out


def sample_values(fmt: str, column: str) -> dict:
    return dict(zip(fmt.split(":"), column.split(":")))


def int_list(value, default=None):
    if value in (None, "", "."):
        return default
    return [int(x) if x not in (".", "") else 0 for x in str(value).split(",")]


# --- PAVE / IMPACT (interface spec section 1) --------------------------------------------------

IMPACT_FIELDS = (
    "Gene", "Transcript", "CanonicalEffect", "CanonicalCodingEffect", "SpliceRegion",
    "HgvsCoding", "HgvsProtein", "OtherReportableEffects", "WorstCodingEffect", "GenesAffected",
)
PAVE_TI_FIELDS = (
    "GeneId", "GeneName", "Transcript", "Effects", "SpliceRegion", "HgvsCoding", "HgvsProtein",
    "RefSeqId", "Exon", "Codon",
)
# Coding consequence classes carried into pVACseq in phase 1 (fusions stay with NEO/LINX).
PHASE1_CODING_EFFECTS = ("MISSENSE", "NONSENSE_OR_FRAMESHIFT", "INFRAME_INDEL")


def parse_impact(value: str) -> dict:
    parts = value.split(",")
    return dict(zip(IMPACT_FIELDS, parts))


def parse_pave_ti(value: str) -> list[dict]:
    out = []
    for entry in value.split(","):
        if not entry:
            continue
        out.append(dict(zip(PAVE_TI_FIELDS, entry.split("|"))))
    return out


def strip_transcript_version(transcript: str) -> str:
    return transcript.split(".")[0]


# --- run layout --------------------------------------------------------------------------------


@dataclass
class CaseLayout:
    """Resolved paths for one oncoanalyser case directory (results/<case>)."""

    results_case_dir: str
    tumour_dna: str = ""
    normal_dna: str = ""
    tumour_rna: str = ""
    files: dict = field(default_factory=dict)

    @classmethod
    def from_results_dir(cls, results_case_dir: str) -> "CaseLayout":
        def one(pattern, required=True):
            hits = sorted(glob.glob(os.path.join(results_case_dir, pattern)))
            if not hits:
                if required:
                    raise FileNotFoundError(f"{pattern} not found under {results_case_dir}")
                return None
            return hits[0]

        somatic = one("purple/*.purple.somatic.vcf.gz")
        samples = vcf_samples(somatic)
        if len(samples) != 2:
            raise ValueError(f"expected [normal, tumour] in {somatic}, got {samples}")
        normal_dna, tumour_dna = samples
        append = one(f"sage_append/{tumour_dna}/*.sage.append.vcf.gz", required=False)
        tumour_rna = ""
        if append:
            append_samples = vcf_samples(append)
            extra = [s for s in append_samples if s not in (normal_dna, tumour_dna)]
            tumour_rna = extra[-1] if extra else ""
        files = {
            "purple_somatic_vcf": somatic,
            "purple_germline_vcf": one("purple/*.purple.germline.vcf.gz", required=False),
            "purple_purity_tsv": one("purple/*.purple.purity.tsv"),
            "purple_qc": one("purple/*.purple.qc", required=False),
            "sage_append_vcf": append,
            "lilac_tsv": one("lilac/*.lilac.tsv", required=False),
            "lilac_qc_tsv": one("lilac/*.lilac.qc.tsv", required=False),
            "isofox_gene_tsv": one("isofox/*.isf.gene_data.tsv", required=False),
            "isofox_transcript_tsv": one("isofox/*.isf.transcript_data.tsv", required=False),
            "neo_data_tsv": one("neo/finder/*.neo_data.tsv", required=False),
            "neo_peptide_scores_tsv": one("neo/scorer/*.neo.peptide_scores.tsv", required=False),
            "neo_neoepitope_tsv": one("neo/scorer/*.neo.neoepitope.tsv", required=False),
            "linx_neoepitope_tsv": one("linx/somatic_annotations/*.linx.neoepitope.tsv", required=False),
            "orange_json": one("orange/*.orange.json", required=False),
        }
        return cls(results_case_dir, tumour_dna, normal_dna, tumour_rna, files)


def read_tsv(path) -> list[dict]:
    with _open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def lilac_alleles(lilac_tsv, *, hla_class=None, prefix=True) -> list[str]:
    """De-duplicated allele list (homozygous alleles are listed twice by LILAC)."""
    seen, out = set(), []
    for row in read_tsv(lilac_tsv):
        a = normalise_allele(row["Allele"], prefix=prefix)
        if hla_class and allele_class(a) != hla_class:
            continue
        if a not in seen:
            seen.add(a)
            out.append(a)
    return out


def immune_escape(orange_json) -> dict:
    with open(orange_json) as fh:
        return json.load(fh).get("immuneEscape", {})


def purity_fields(purity_tsv) -> dict:
    return read_tsv(purity_tsv)[0]
