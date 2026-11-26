#!/usr/bin/env python3
"""Build the clonal haematopoiesis hotspot VCF for SAGE.

SAGE evaluates every position in the file passed to -hotspots whatever the panel says. The
cancer-predisposition list shipped with the Hartwig resources therefore injects BRCA1, BRCA2, ATM, NF1
and similar genes into a CHIP call set, and forces calls on one or two reads at those positions. This
replaces it with the Hartwig somatic hotspots that the curated whitelist recognises.

A hotspot is kept when its gene is curated and its protein change is a listed substitution, a missense
inside a named domain, or a truncating change in a gene on the loss-of-function list. Hartwig annotates
each hotspot on one transcript, which need not be the curated one, so the substitution test accepts
either numbering; the authors' whitelist carries both for genes such as GNAS. The hotspot tier only
raises sensitivity at a site. Whether a call counts is decided later by chip_curation.py on the curated
transcript.

Usage: build_chip_hotspots.py <KnownHotspots.somatic.38.vcf(.gz)> <panel_dir> <out.vcf>
"""
import gzip
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chip_curation import CuratedDefinition, MISSENSE_RANGES, THREE_TO_ONE  # noqa: E402


def one_letter(p):
    p = p[2:] if p.startswith("p.") else p
    m = re.match(r"^([A-Z])(\d+)([A-Z])$", p)
    if m:
        return p
    m = re.match(r"^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$", p)
    if m and m.group(1) in THREE_TO_ONE and m.group(3) in THREE_TO_ONE:
        return f"{THREE_TO_ONE[m.group(1)]}{m.group(2)}{THREE_TO_ONE[m.group(3)]}"
    return None


def main(source, panel_dir, out_vcf):
    cd = CuratedDefinition(panel_dir)
    opener = gzip.open if source.endswith(".gz") else open
    kept, genes = [], {}
    n = {"read": 0, "curated_gene": 0, "listed_missense": 0, "domain_missense": 0, "truncating": 0}
    with opener(source, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            n["read"] += 1
            f = line.rstrip("\n").split("\t")
            m = re.search(r"input=([^|;]+)\|([^|;]+)\|([^;]+)", f[7])
            if not m:
                continue
            gene, transcript, prot = m.group(1), m.group(2), m.group(3)
            if gene not in cd.genes:
                continue
            n["curated_gene"] += 1
            why = None
            short = one_letter(prot)
            if prot.endswith("*") or "fs" in prot:
                if gene in cd.lof_genes:
                    why = "truncating"
            elif short:
                if short in cd.missense.get(gene, ()):
                    why = "listed_missense"
                elif (gene in MISSENSE_RANGES
                      and transcript == cd.transcripts.get(gene, {}).get("transcript")
                      and any(a <= int(short[1:-1]) <= b for a, b in MISSENSE_RANGES[gene])):
                    why = "domain_missense"
            if why:
                n[why] += 1
                kept.append(f)
                genes[gene] = genes.get(gene, 0) + 1
    with open(out_vcf, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write('##INFO=<ID=input,Number=1,Type=String,Description="Source annotation">\n')
        out.write('##INFO=<ID=sources,Number=1,Type=String,Description="Evidence sources">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for f in kept:
            out.write("\t".join(f) + "\n")
    print(f"hotspot records read {n['read']}  in a curated gene {n['curated_gene']}")
    print(f"kept {len(kept)}: listed missense {n['listed_missense']}, domain missense "
          f"{n['domain_missense']}, truncating {n['truncating']}, across {len(genes)} genes")
    print("per gene: " + " ".join(f"{g}={c}" for g, c in sorted(genes.items())))
    if not kept:
        raise SystemExit("no hotspot kept")


if __name__ == "__main__":
    main(*sys.argv[1:4])
