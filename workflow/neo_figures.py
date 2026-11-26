#!/usr/bin/env python3
"""Render the neoantigen figures for the ORANGE supplement.

Three panels are produced, following docs/neoantigen_figure_style.md:

  binding.png    mutant and wild-type binding affinity for each candidate
  anchors.png    position of the mutated residue within each peptide
  criteria.png   the prioritisation criteria met by each candidate

Every colour comes from ORANGE ReportResources. Hue carries the tier, the ORANGE_1..6 ladder is reserved
for ordered quantities, and the fill of a marker carries the state of the RNA evidence. Type matches the
report at 8 pt titles, 7 pt axis text and 6 pt annotation, with axes drawn in mid grey at 0.25 pt so they
sit at the same weight as the rules between table cells.

Requires matplotlib. The supplement renderer itself additionally requires reportlab and pypdf.
"""
from __future__ import annotations

import argparse
import csv
import os

import matplotlib as mpl
import numpy as np

mpl.use("Agg")
import matplotlib.pyplot as plt                                  # noqa: E402
from matplotlib.lines import Line2D                              # noqa: E402

ORANGE = (242 / 255, 139 / 255, 31 / 255)
LADDER = [(255 / 255, 165 / 255, 0), (235 / 255, 155 / 255, 0), (215 / 255, 145 / 255, 0),
          (195 / 255, 135 / 255, 0), (175 / 255, 125 / 255, 0), (155 / 255, 115 / 255, 0)]
BLUE = (38 / 255, 90 / 255, 166 / 255)
DARK = (39 / 255, 47 / 255, 50 / 255)
MID = (101 / 255, 106 / 255, 108 / 255)
LIGHT = (211 / 255, 211 / 255, 211 / 255)
GAINS = (220 / 255, 220 / 255, 220 / 255)
SMOKE = (240 / 255, 240 / 255, 240 / 255)

TIER_COLOUR = {"Prioritise": ORANGE, "Review": BLUE, "Hold": MID, "Exclude": LIGHT}
TIER_ORDER = ["Prioritise", "Review", "Hold", "Exclude"]
WIDTH_IN = 7.43          # 535 pt of content width

RC = {"font.family": "sans-serif",
      "font.sans-serif": ["Nimbus Sans", "Nimbus Sans L", "Helvetica", "Arial", "DejaVu Sans"],
      "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 8, "xtick.labelsize": 7,
      "ytick.labelsize": 7, "axes.edgecolor": MID, "axes.linewidth": 0.25, "xtick.color": MID,
      "ytick.color": MID, "text.color": DARK, "axes.labelcolor": DARK, "legend.fontsize": 6,
      "xtick.major.width": 0.25, "ytick.major.width": 0.25, "figure.dpi": 300,
      "savefig.facecolor": "white"}


def _open_frame(ax, left=True):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if not left:
        ax.spines["left"].set_visible(False)


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def read_candidates(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def read_aggregated(path):
    if not path or not os.path.isfile(path):
        return {}
    return {r["ID"]: r for r in csv.DictReader(open(path), delimiter="\t")}


def panel_binding(cands, outfile, limit=22):
    """Mutant and wild-type binding affinity, one row per candidate."""
    rows = [r for r in cands if r["tier"] in ("Prioritise", "Review") and _f(r["ic50_wt"]) and _f(r["ic50_mt"])]
    rows = sorted(rows, key=lambda r: _f(r["ic50_mt"]))[:limit]
    with mpl.rc_context(RC):
        fig, ax = plt.subplots(figsize=(WIDTH_IN, 3.4))
        y = np.arange(len(rows))[::-1]
        for yi, r in zip(y, rows):
            mt, wt = _f(r["ic50_mt"]), _f(r["ic50_wt"])
            colour = TIER_COLOUR[r["tier"]]
            ax.annotate("", xy=(mt, yi), xytext=(wt, yi),
                        arrowprops=dict(arrowstyle="-|>", color=colour, lw=0.7, shrinkA=0,
                                        shrinkB=0, mutation_scale=6, alpha=0.85))
            ax.scatter(wt, yi, s=11, facecolor="white", edgecolor=MID, lw=0.6, zorder=3)
            ax.scatter(mt, yi, s=18, facecolor=colour, edgecolor=colour, lw=0.6, zorder=4)
        ax.set_xscale("log")
        ax.set_xlim(8, 2e4)
        ax.axvline(500, color=MID, lw=0.5, ls=(0, (4, 3)))
        ax.text(530, len(rows) - 0.5, "500 nM", fontsize=6, color=MID, va="top")
        ax.set_yticks(y)
        ax.set_yticklabels([f"{r['gene']} {r['aa_change']}" for r in rows], fontsize=6)
        ax.tick_params(axis="y", length=0)
        ax.set_xlabel("Predicted binding affinity (nM, logarithmic scale)")
        _open_frame(ax, left=False)
        ax.set_title("Mutant And Wild-Type Binding Affinity", loc="left", pad=6)
        handles = [Line2D([], [], marker="o", ls="", mfc="white", mec=MID, ms=3.4,
                           label="wild-type peptide"),
                   Line2D([], [], marker="o", ls="", mfc=ORANGE, mec=ORANGE, ms=4.2,
                           label="mutant peptide, prioritised"),
                   Line2D([], [], marker="o", ls="", mfc=BLUE, mec=BLUE, ms=4.2,
                           label="mutant peptide, under review")]
        ax.legend(handles=handles, frameon=False, loc="lower right")
        fig.tight_layout()
        fig.savefig(outfile, dpi=300)
        plt.close(fig)
    return outfile


def panel_anchors(cands, aggregated, outfile, limit=14):
    """Position of the mutated residue within each peptide, against the canonical anchor positions."""
    chosen = [r for r in cands if r["tier"] == "Prioritise"]
    chosen += sorted([r for r in cands if r["tier"] == "Review"], key=lambda r: _f(r["ic50_mt"]))
    recs = []
    for r in chosen:
        a = aggregated.get(r["variant"])
        if not a:
            continue
        pos = _f(a.get("Pos"))
        recs.append((r, a["Best Peptide"], int(pos) if pos else None))
        if len(recs) >= limit:
            break
    if not recs:
        return None
    longest = max(len(p) for _, p, _ in recs)
    with mpl.rc_context(RC):
        fig, ax = plt.subplots(figsize=(WIDTH_IN, 3.3))
        for i, (r, peptide, pos) in enumerate(recs):
            yi = len(recs) - 1 - i
            anchors = {2, len(peptide)}
            for j, residue in enumerate(peptide, start=1):
                mutated = pos == j
                face = ORANGE if mutated else (GAINS if j in anchors else "white")
                ax.add_patch(plt.Rectangle((j - 0.5, yi - 0.42), 1, 0.84, facecolor=face,
                                           edgecolor=MID, lw=0.25, zorder=2))
                ax.text(j, yi, residue, ha="center", va="center", fontsize=6, zorder=3,
                        color="white" if mutated else DARK,
                        fontweight="bold" if mutated else "normal")
            ax.text(longest + 1.1, yi, f"{r['gene']} {r['aa_change']}", fontsize=6, va="center",
                    color=DARK)
            ax.text(-0.2, yi, r["allele"], fontsize=6, va="center", ha="right",
                    color=TIER_COLOUR[r["tier"]] if r["tier"] != "Exclude" else MID)
        ax.set_xlim(-3.6, longest + 5.2)
        ax.set_ylim(-0.8, len(recs) - 0.2)
        ax.set_xticks(range(1, longest + 1))
        ax.set_xticklabels([f"P{i}" for i in range(1, longest + 1)], fontsize=6)
        ax.set_yticks([])
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        handles = [plt.Rectangle((0, 0), 1, 1, facecolor=ORANGE, edgecolor=MID, lw=0.25,
                                 label="mutated residue"),
                   plt.Rectangle((0, 0), 1, 1, facecolor=GAINS, edgecolor=MID, lw=0.25,
                                 label="canonical anchor position"),
                   plt.Rectangle((0, 0), 1, 1, facecolor="white", edgecolor=MID, lw=0.25,
                                 label="position exposed to the T-cell receptor")]
        ax.legend(handles=handles, frameon=False, loc="lower center", bbox_to_anchor=(0.45, -0.26),
                  ncols=3, handlelength=0.9, columnspacing=1.2)
        ax.set_title("Position Of The Mutated Residue Within The Peptide", loc="left", pad=6)
        fig.tight_layout()
        fig.savefig(outfile, dpi=300)
        plt.close(fig)
    return outfile


CRITERIA = [
    ("binding affinity within threshold", lambda r: _f(r["ic50_mt"]) is not None and _f(r["ic50_mt"]) <= 500),
    ("percentile within threshold", lambda r: _f(r["percentile_mt"]) is not None and _f(r["percentile_mt"]) <= 2),
    ("mutant allele present in RNA", lambda r: r["rna_status"] == "SUPPORTED"),
    ("clonal", lambda r: r["subclonal"] != "True"),
    ("wild-type binds less well", lambda r: (_f(r.get("fold_change")) >= 2.0) if _f(r.get("fold_change")) else None),
    ("restricting allele retained", lambda r: not (r["allele_lost"] == "True" or r["allele_disrupted"] == "True")),
    ("absent from the matched normal", lambda r: int(r["normal_alt_reads"] or 0) == 0),
    ("transcript support level 1", lambda r: r["tsl"] == "1"),
    ("listed in the Cancer Gene Census", lambda r: r["gene_of_interest"] == "True"),
]


def panel_criteria(cands, outfile, limit=16):
    """The prioritisation criteria met by each candidate."""
    rows = [r for r in cands if r["tier"] == "Prioritise"]
    rows += sorted([r for r in cands if r["tier"] == "Review"], key=lambda r: _f(r["ic50_mt"]))
    rows = rows[:limit]
    with mpl.rc_context(RC):
        fig, ax = plt.subplots(figsize=(WIDTH_IN, 3.6))
        for i, r in enumerate(rows):
            yi = len(rows) - 1 - i
            for j, (_name, test) in enumerate(CRITERIA):
                try:
                    value = test(r)
                except Exception:
                    value = None
                if value is True:
                    face, edge, mark = TIER_COLOUR[r["tier"]], TIER_COLOUR[r["tier"]], ""
                elif value is False:
                    face, edge, mark = SMOKE, LIGHT, ""
                else:
                    face, edge, mark = "white", MID, "\u2013"
                ax.add_patch(plt.Rectangle((j - 0.4, yi - 0.4), 0.8, 0.8, facecolor=face,
                                           edgecolor=edge, lw=0.3, zorder=2))
                if mark:
                    ax.text(j, yi, mark, ha="center", va="center", fontsize=6, color=MID, zorder=3)
            ax.text(-0.75, yi, f"{r['gene']} {r['aa_change']}", ha="right", va="center", fontsize=6,
                    color=TIER_COLOUR[r["tier"]] if r["tier"] != "Exclude" else MID)
        ax.set_xlim(-5.4, len(CRITERIA) - 0.3)
        ax.set_ylim(-0.8, len(rows) - 0.2)
        ax.set_xticks(range(len(CRITERIA)))
        ax.set_xticklabels([c[0] for c in CRITERIA], rotation=35, ha="left", fontsize=6)
        ax.xaxis.set_ticks_position("top")
        ax.tick_params(length=0)
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        handles = [plt.Rectangle((0, 0), 1, 1, facecolor=ORANGE, edgecolor=ORANGE,
                                 label="criterion met, coloured by tier"),
                   plt.Rectangle((0, 0), 1, 1, facecolor=SMOKE, edgecolor=LIGHT,
                                 label="criterion unmet"),
                   plt.Rectangle((0, 0), 1, 1, facecolor="white", edgecolor=MID,
                                 label="criterion not assessed")]
        ax.legend(handles=handles, frameon=False, loc="lower left", bbox_to_anchor=(-0.08, -0.22),
                  ncols=3, handlelength=0.9, columnspacing=1.2)
        fig.tight_layout()
        fig.savefig(outfile, dpi=300)
        plt.close(fig)
    return outfile


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidates-tsv", required=True)
    ap.add_argument("--aggregated-tsv", default=None,
                    help="pVACseq aggregated report, for the mutated-residue position")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--prefix", default="neo")
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    cands = read_candidates(a.candidates_tsv)
    aggregated = read_aggregated(a.aggregated_tsv)
    written = [panel_binding(cands, os.path.join(a.outdir, f"{a.prefix}_binding.png")),
               panel_anchors(cands, aggregated, os.path.join(a.outdir, f"{a.prefix}_anchors.png")),
               panel_criteria(cands, os.path.join(a.outdir, f"{a.prefix}_criteria.png"))]
    for path in written:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
