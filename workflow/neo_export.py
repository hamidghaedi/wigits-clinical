#!/usr/bin/env python3
"""Export a native oncoanalyser case into pVACseq inputs.

Implements docs/neoantigen_interface.md section 9. Reads a published `results/<case>` directory and
writes, into an output directory:

  <case>.somatic.selected.vcf.gz   PASS somatic variants with a phase-1 coding effect, carrying
                                   tumour/normal DNA and tumour RNA depths in the FORMAT fields
                                   pVACseq reads (AD/DP/AF and RAD/RDP/RAF) plus GX/TX expression
  <case>.proximal.vcf.gz           somatic variants sharing a SAGE phase set with a selected variant,
                                   plus PASS germline variants within the peptide window
  <case>.hla_alleles.txt           de-duplicated HLA alleles (class I by default)
  <case>.interface_manifest.json   source paths + sha256, record counts, tool versions, rule version

The VEP step is deliberately NOT run here: it needs a container and is submitted separately by
workflow/neo_export.sbatch, which runs this script first (`prepare`) and then VEP with the arguments
this script records in the manifest.

Read-only with respect to the run directory.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import sys
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import neo_interface as ni  # noqa: E402

RULES_VERSION = "docs/neoantigen_interface.md@2026-11-23"

# Phase 1 carries point mutations and small indels; fusions stay with NEO/LINX (spec section 7).
CODING_EFFECTS = {"MISSENSE", "NONSENSE_OR_FRAMESHIFT", "INFRAME_INDEL"}
# Window either side of a variant in which a germline variant can change the peptide (spec section 5).
PEPTIDE_WINDOW_BP = 100


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _dp_ad_af(values: dict):
    dp = values.get("DP")
    dp = int(dp) if dp not in (None, "", ".") else None
    ad = ni.int_list(values.get("AD"), default=[0, 0])
    alt = ad[1] if len(ad) > 1 else 0
    af = values.get("AF")
    try:
        af = float(af)
    except (TypeError, ValueError):
        af = (alt / dp) if dp else 0.0
    return dp, alt, af


class Exporter:
    def __init__(self, results_case_dir: str, outdir: str, case: str | None = None,
                 hla_class: str = "I", keep_non_coding: bool = False):
        self.layout = ni.CaseLayout.from_results_dir(results_case_dir)
        self.outdir = outdir
        self.case = case or os.path.basename(os.path.normpath(results_case_dir))
        self.hla_class = hla_class
        self.keep_non_coding = keep_non_coding
        os.makedirs(outdir, exist_ok=True)
        self.counts: dict[str, int] = {}
        self.rna_states: dict[str, int] = {"SUPPORTED": 0, "COVERED_NO_ALT": 0, "NOT_COVERED": 0}

    # -- RNA genotypes ---------------------------------------------------------------------
    def read_rna_support(self) -> dict:
        """key (chrom,pos,ref,alt) -> (rdp, rad, raf, status) from the SAGE_APPEND RNA column."""
        path = self.layout.files.get("sage_append_vcf")
        if not path:
            return {}
        samples = ni.vcf_samples(path)
        rna = self.layout.tumour_rna or samples[-1]
        idx = samples.index(rna)
        out = {}
        with ni._open(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                vals = ni.sample_values(f[8], f[9 + idx])
                dp, alt, af = _dp_ad_af(vals)
                status = ni.rna_coverage_status(dp, alt)
                out[(f[0], f[1], f[3], f[4])] = (dp if dp is not None else 0, alt, af, status)
        self.counts["sage_append_records"] = len(out)
        return out

    # -- expression ------------------------------------------------------------------------
    def read_expression(self):
        genes, transcripts = {}, {}
        gp = self.layout.files.get("isofox_gene_tsv")
        if gp:
            for row in ni.read_tsv(gp):
                genes[row["GeneId"]] = row["AdjTPM"]
        tp = self.layout.files.get("isofox_transcript_tsv")
        if tp:
            for row in ni.read_tsv(tp):
                transcripts[ni.strip_transcript_version(row["TransName"])] = row["AdjTPM"]
        self.counts["isofox_genes"] = len(genes)
        self.counts["isofox_transcripts"] = len(transcripts)
        return genes, transcripts

    # -- main VCF pass ---------------------------------------------------------------------
    def write_vcf(self):
        rna = self.read_rna_support()
        genes, transcripts = self.read_expression()
        src = self.layout.files["purple_somatic_vcf"]
        samples = ni.vcf_samples(src)
        n_idx, t_idx = samples.index(self.layout.normal_dna), samples.index(self.layout.tumour_dna)
        out_vcf = os.path.join(self.outdir, f"{self.case}.somatic.selected.vcf")
        kept_keys, phase_sets = [], set()
        total = passed = selected = 0

        extra_headers = [
            '##FORMAT=<ID=RDP,Number=1,Type=Integer,Description="Tumour RNA depth (SAGE_APPEND RNA sample)">',
            '##FORMAT=<ID=RAD,Number=R,Type=Integer,Description="Tumour RNA allelic depths (ref,alt)">',
            '##FORMAT=<ID=RAF,Number=1,Type=Float,Description="Tumour RNA allele frequency">',
            '##INFO=<ID=RNA_STATUS,Number=1,Type=String,Description="SUPPORTED|COVERED_NO_ALT|NOT_COVERED per docs/neoantigen_interface.md">',
            '##INFO=<ID=GX,Number=.,Type=String,Description="Gene expression (ISOFOX AdjTPM) as GeneId|TPM">',
            '##INFO=<ID=TX,Number=.,Type=String,Description="Transcript expression (ISOFOX AdjTPM) as TransName|TPM">',
            f'##neoExport=<rules="{RULES_VERSION}",case="{self.case}",source="{os.path.basename(src)}">',
        ]

        with ni._open(src) as fh, open(out_vcf, "w") as out:
            for line in fh:
                if line.startswith("##"):
                    out.write(line)
                    continue
                if line.startswith("#CHROM"):
                    for h in extra_headers:
                        out.write(h + "\n")
                    cols = line.rstrip("\n").split("\t")
                    out.write("\t".join(cols[:9] + [samples[n_idx], samples[t_idx]]) + "\n")
                    continue
                total += 1
                f = line.rstrip("\n").split("\t")
                if f[6] != "PASS":
                    continue
                passed += 1
                info = ni.parse_info(f[7])
                impact = ni.parse_impact(info["IMPACT"]) if "IMPACT" in info else {}
                effect = impact.get("CanonicalCodingEffect", "")
                if not self.keep_non_coding and effect not in CODING_EFFECTS:
                    continue
                selected += 1
                key = (f[0], f[1], f[3], f[4])
                kept_keys.append(key)
                for ps in ni.int_list(info.get("LPS"), default=[]) or []:
                    phase_sets.add(ps)

                rdp, rad_alt, raf, status = rna.get(key, (0, 0, 0.0, "NOT_COVERED"))
                self.rna_states[status] += 1
                gene_id = impact.get("Gene", "")
                tx = ni.strip_transcript_version(impact.get("Transcript", ""))
                gx = f"{gene_id}|{genes[gene_id]}" if gene_id in genes else ""
                txv = f"{tx}|{transcripts[tx]}" if tx in transcripts else ""
                add = [f"RNA_STATUS={status}"]
                if gx:
                    add.append(f"GX={gx}")
                if txv:
                    add.append(f"TX={txv}")
                f[7] = ";".join([f[7]] + add)

                fmt = f[8].split(":")
                norm = ni.sample_values(f[8], f[9 + n_idx])
                tum = ni.sample_values(f[8], f[9 + t_idx])
                rad_ref = max(rdp - rad_alt, 0)
                fmt_new = fmt + ["RDP", "RAD", "RAF"]
                norm_new = [norm.get(k, ".") for k in fmt] + [".", ".", "."]
                tum_new = [tum.get(k, ".") for k in fmt] + [str(rdp), f"{rad_ref},{rad_alt}", f"{raf:.4f}"]
                out.write("\t".join(f[:8] + [":".join(fmt_new), ":".join(norm_new), ":".join(tum_new)]) + "\n")

        self.counts.update({"somatic_total": total, "somatic_pass": passed, "selected": selected,
                            "phase_sets": len(phase_sets)})
        return out_vcf, set(kept_keys), phase_sets

    # -- proximal variants -----------------------------------------------------------------
    def write_proximal(self, kept_keys, phase_sets):
        """Single-sample VCF: pVACseq rejects a proximal VCF carrying more than one sample."""
        out_vcf = os.path.join(self.outdir, f"{self.case}.proximal.vcf")
        positions = {(c, int(p)) for c, p, _, _ in kept_keys}
        n_som = n_germ = 0
        src = self.layout.files["purple_somatic_vcf"]
        samples = ni.vcf_samples(src)
        t_idx = samples.index(self.layout.tumour_dna)
        with open(out_vcf, "w") as out:
            with ni._open(src) as fh:
                for line in fh:
                    if line.startswith("##"):
                        out.write(line)
                        continue
                    if line.startswith("#CHROM"):
                        cols = line.rstrip("\n").split("\t")
                        out.write("\t".join(cols[:9] + [self.layout.tumour_dna]) + "\n")
                        continue
                    f = line.rstrip("\n").split("\t")
                    if f[6] != "PASS":
                        continue
                    key = (f[0], f[1], f[3], f[4])
                    if key in kept_keys:
                        continue
                    info = ni.parse_info(f[7])
                    lps = set(ni.int_list(info.get("LPS"), default=[]) or [])
                    if lps & phase_sets:
                        out.write("\t".join(f[:9] + [f[9 + t_idx]]) + "\n")
                        n_som += 1
            germ = self.layout.files.get("purple_germline_vcf")
            if germ:
                g_samples = ni.vcf_samples(germ)
                g_idx = g_samples.index(self.layout.tumour_dna) if self.layout.tumour_dna in g_samples else 0
                with ni._open(germ) as fh:
                    for line in fh:
                        if line.startswith("#"):
                            continue
                        f = line.rstrip("\n").split("\t")
                        if f[6] != "PASS":
                            continue
                        pos = int(f[1])
                        if any(f[0] == c and abs(pos - p) <= PEPTIDE_WINDOW_BP for c, p in positions):
                            out.write("\t".join(f[:9] + [f[9 + g_idx]]) + "\n")
                            n_germ += 1
        self.counts.update({"proximal_somatic": n_som, "proximal_germline": n_germ})
        return out_vcf

    # -- alleles ---------------------------------------------------------------------------
    def write_alleles(self):
        path = os.path.join(self.outdir, f"{self.case}.hla_alleles.txt")
        lilac = self.layout.files.get("lilac_tsv")
        alleles: list[str] = []
        qc_status = None
        if lilac:
            alleles = ni.lilac_alleles(lilac, hla_class=None if self.hla_class == "all" else self.hla_class)
        qc = self.layout.files.get("lilac_qc_tsv")
        if qc:
            qc_status = ni.read_tsv(qc)[0].get("Status")
        with open(path, "w") as fh:
            fh.write(",".join(alleles) + "\n")
        self.counts["hla_alleles"] = len(alleles)
        return path, alleles, qc_status

    # -- manifest --------------------------------------------------------------------------
    def write_manifest(self, outputs: dict, alleles, qc_status):
        sources = {}
        for name, path in self.layout.files.items():
            if path and os.path.isfile(path):
                sources[name] = {"path": path, "sha256": sha256(path), "bytes": os.path.getsize(path)}
        versions = {}
        vfile = os.path.join(os.path.dirname(os.path.dirname(self.layout.results_case_dir)),
                             "pipeline_info", "nf_core_oncoanalyser_software_mqc_versions.yml")
        if os.path.isfile(vfile):
            versions = {"software_versions_yml": {"path": vfile, "sha256": sha256(vfile)}}
        purity = ni.purity_fields(self.layout.files["purple_purity_tsv"])
        escape = {}
        if self.layout.files.get("orange_json"):
            escape = ni.immune_escape(self.layout.files["orange_json"])
        manifest = {
            "case": self.case,
            "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "rules_version": RULES_VERSION,
            "samples": {"tumour_dna": self.layout.tumour_dna, "normal_dna": self.layout.normal_dna,
                         "tumour_rna": self.layout.tumour_rna},
            "counts": self.counts,
            "rna_coverage_states": self.rna_states,
            "hla": {"alleles": alleles, "lilac_qc_status": qc_status, "class": self.hla_class},
            "purity": {k: purity.get(k) for k in
                        ("purity", "ploidy", "status", "wholeGenomeDuplication", "msStatus",
                         "tml", "tmlStatus", "tmbPerMb", "tmbStatus")},
            "immune_escape": escape,
            "outputs": outputs,
            "sources": sources,
            "oncoanalyser": versions,
        }
        path = os.path.join(self.outdir, f"{self.case}.interface_manifest.json")
        with open(path, "w") as fh:
            json.dump(manifest, fh, indent=1)
        return path, manifest

    def run(self):
        vcf, kept, phase_sets = self.write_vcf()
        prox = self.write_proximal(kept, phase_sets)
        allele_file, alleles, qc_status = self.write_alleles()
        outputs = {"somatic_selected_vcf": vcf, "proximal_vcf": prox, "hla_alleles": allele_file}
        manifest_path, manifest = self.write_manifest(outputs, alleles, qc_status)
        return manifest_path, manifest


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-case-dir", required=True, help="published oncoanalyser results/<case> directory")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--case", default=None)
    ap.add_argument("--hla-class", default="I", choices=["I", "II", "all"])
    ap.add_argument("--keep-non-coding", action="store_true",
                    help="keep every PASS variant instead of the phase-1 coding classes")
    args = ap.parse_args(argv)
    path, manifest = Exporter(args.results_case_dir, args.outdir, args.case,
                              args.hla_class, args.keep_non_coding).run()
    print(json.dumps({"manifest": path, "counts": manifest["counts"],
                      "rna_coverage_states": manifest["rna_coverage_states"],
                      "hla_alleles": manifest["hla"]["alleles"]}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
