#!/usr/bin/env python3
"""Build a clonal haematopoiesis hotspot VCF for SAGE.

SAGE evaluates every position in the file passed to -hotspots regardless of the panel, so the
cancer-predisposition list shipped with the Hartwig resources injects BRCA1, BRCA2, MSH2, TSC1 and
similar genes into a CHIP call set. This script replaces it with hotspots that belong to the CHIP
panel and sit at positions the curated definition of Vlasschaert et al. recognises.

Input  hotspots_all.tsv            Hartwig KnownHotspots.somatic.38, header plus records
       chip_genes_vlasschaert.tsv  the curated definition
       panel82.txt                 the panel
Output KnownHotspots.chip.38.vcf   positions to be evaluated at hotspot sensitivity
"""
import csv, re, sys

THREE = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G","His":"H",
         "Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S","Thr":"T","Trp":"W",
         "Tyr":"Y","Val":"V"}
NORM = {"MLL/KMT2A":"KMT2A","MLL2/KMT2D":"KMT2D"}

curated, panel = {}, set()
with open(sys.argv[2]) as fh:
    for row in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
        curated[NORM.get(row["gene"], row["gene"])] = {
            "trunc": row["admits_truncating"] == "true",
            "missense": set(x for x in row["curated_missense"].split(",") if x)}
with open(sys.argv[3]) as fh:
    panel = {l.strip() for l in fh if l.strip()}

def to_short(p):
    p = p.lstrip("p.")
    m = re.match(r"^([A-Z])(\d+)([A-Z])$", p)
    if m:
        return m.group(0)
    m = re.match(r"^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$", p)
    if m and m.group(1) in THREE and m.group(3) in THREE:
        return f"{THREE[m.group(1)]}{m.group(2)}{THREE[m.group(3)]}"
    return None

kept, seen_genes = [], set()
stats = {"total": 0, "in_panel": 0, "curated_missense": 0, "truncating": 0}
with open(sys.argv[1]) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        stats["total"] += 1
        f = line.rstrip("\n").split("\t")
        m = re.search(r"input=([^|]+)\|([^|]+)\|(\S+?)(?:;|$)", f[7])
        if not m:
            continue
        gene, prot = NORM.get(m.group(1), m.group(1)), m.group(3)
        if gene not in panel:
            continue
        stats["in_panel"] += 1
        rule = curated.get(gene)
        if rule is None:                       # one of the eight genes outside the definition
            continue
        if prot.endswith("*") or "fs" in prot or "del" in prot or "ins" in prot:
            if rule["trunc"]:
                stats["truncating"] += 1
                kept.append((f, gene, prot)); seen_genes.add(gene)
            continue
        short = to_short(prot)
        if short and short in rule["missense"]:
            stats["curated_missense"] += 1
            kept.append((f, gene, prot)); seen_genes.add(gene)

with open(sys.argv[4], "w") as out:
    out.write("##fileformat=VCFv4.2\n")
    out.write('##INFO=<ID=input,Number=1,Type=String,Description="Source annotation">\n')
    out.write('##INFO=<ID=sources,Number=1,Type=String,Description="Evidence sources">\n')
    out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for f, _, _ in kept:
        out.write("\t".join(f) + "\n")

print(f"hotspot records read          {stats['total']}")
print(f"in the CHIP panel             {stats['in_panel']}")
print(f"kept, curated missense        {stats['curated_missense']}")
print(f"kept, truncating in a gene that admits them  {stats['truncating']}")
print(f"written                       {len(kept)} across {len(seen_genes)} genes")
print("genes: " + " ".join(sorted(seen_genes)))
