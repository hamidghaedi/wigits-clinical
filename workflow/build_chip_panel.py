#!/usr/bin/env python3
"""Write the CHIP gene panel in the driver gene panel format SAGE and PAVE read.

SAGE calls germline-mode variants only in genes whose reportGermlineVariant is set, which in the
Hartwig panel excludes every core CHIP gene except TP53. This takes the genes of panels/chip_genes.tsv,
copies the Hartwig row where one exists, adds a row where none does, and switches germline reporting on.

Usage: build_chip_panel.py <chip_genes.tsv> <DriverGenePanel.38.tsv> <out.tsv>
"""
import csv
import sys


def main(genes_tsv, driver_panel, out_tsv):
    want = []
    for line in open(genes_tsv):
        if line.startswith("#") or line.startswith("gene\t"):
            continue
        g = line.split("\t")[0].strip()
        if g:
            want.append(g)
    with open(driver_panel) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    fields = list(rows[0].keys())
    out, seen = [], set()
    for r in rows:
        if r["gene"] in want:
            r = dict(r, reportGermlineVariant="ANY", reportGermlineHotspot="ANY")
            out.append(r)
            seen.add(r["gene"])
    for g in sorted(set(want) - seen):
        r = {k: "" for k in fields}
        r.update({"gene": g, "reportMissense": "true", "reportNonsense": "true", "reportSplice": "true",
                  "reportDeletion": "false", "reportHetDeletion": "false", "hetDeletionThreshold": "0",
                  "reportLoh": "false", "reportDisruption": "false", "reportAmplification": "false",
                  "amplificationRatio": "0", "reportSomaticHotspot": "true", "likelihoodType": "TSG",
                  "reportGermlineVariant": "ANY", "reportGermlineHotspot": "ANY",
                  "reportGermlineDisruption": "NONE", "reportGermlineDeletion": "NONE",
                  "reportGermlineAmplification": "NONE", "reportPGX": "false",
                  "reportHighExpression": "false", "reportNovelSpliceJunction": "false"})
        out.append(r)
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(out)
    print(f"panel genes {len(out)}  from the Hartwig panel {len(seen)}  added {len(out) - len(seen)}")
    if len(out) != len(set(want)):
        raise SystemExit("panel is missing genes")


if __name__ == "__main__":
    main(*sys.argv[1:4])
