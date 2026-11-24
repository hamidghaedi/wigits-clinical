#!/usr/bin/env python3
"""Render a neoantigen supplement in the ORANGE report's own visual language.

ORANGE is left completely untouched: this writes separate pages and, optionally, concatenates them
after the original PDF into <sample>.orange.extended.pdf. Nothing in oncoanalyser or in the ORANGE
container is modified, so a pipeline upgrade cannot break it and the original report stays byte-identical.

Style is taken from ORANGE 5.0.1 `ReportResources.java` (read at tag orange-v5.0.1) and verified against a
rendered report:

  fonts      Nimbus Sans L Regular/Bold.  Nimbus Sans L is URW's clone of Helvetica and is metrically
             identical to it, so reportlab's built-in Helvetica reproduces the native line breaks and
             column widths exactly. If the real TTFs are present they are used instead.
  palette    PALETTE_ORANGE (242,139,31), ORANGE_1..6, DARK/MID/LIGHT/GAINSBORO/SMOKE grey, BLUE.
  type       chapter title bold 10 orange · table title bold 8 orange · table header regular 7 mid-grey,
             upper-cased, fixed 23 pt row · content regular 7 dark-grey with a 0.25 pt mid-grey rule
             below · key/value regular 7 mid-grey · subtext 6 · disclaimer 6 mid-grey.
  page       A4; margins top 100 / left 30 / right 30 / bottom 40; header band 60; footer 20.
  furniture  orange-slice logo top-left, "ORANGE Report (Research Use Only)" with ORANGE in orange,
             orange side panel carrying SAMPLE and the sample name in white, page number bold 7 orange
             and the RUO disclaimer in the footer.

Candidate-grade rows are shaded SMOKE_GREY, which is ORANGE's own idiom for "candidate, not confirmed"
(`candidateDriverBackground`); here it marks the Review and Hold tiers.
"""
from __future__ import annotations

import argparse
import csv
import json
import os

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.utils import ImageReader
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (BaseDocTemplate, Frame, Image, KeepTogether, PageTemplate,
                                Paragraph, Spacer, Table, TableStyle)
from reportlab.lib.styles import ParagraphStyle

# -- palette (ReportResources) -------------------------------------------------------------------
def _rgb(r, g, b):
    return colors.Color(r / 255.0, g / 255.0, b / 255.0)

PALETTE_ORANGE = _rgb(242, 139, 31)
ORANGE_LADDER = [_rgb(255, 165, 0), _rgb(235, 155, 0), _rgb(215, 145, 0),
                 _rgb(195, 135, 0), _rgb(175, 125, 0), _rgb(155, 115, 0)]
PALETTE_DARK_GREY = _rgb(39, 47, 50)
PALETTE_MID_GREY = _rgb(101, 106, 108)
PALETTE_LIGHT_GREY = _rgb(211, 211, 211)
PALETTE_SMOKE_GREY = _rgb(240, 240, 240)
PALETTE_BLUE = _rgb(38, 90, 166)

MARGIN_TOP, MARGIN_LEFT, MARGIN_RIGHT, MARGIN_BOTTOM = 100, 30, 30, 40
HEADER_H, FOOTER_H = 60, 20
HEADER_ROW_H = 23          # Cells.createHeader fixed height
CONTENT_W = A4[0] - MARGIN_LEFT - MARGIN_RIGHT

DISCLAIMER = ("This report is for 'Research Use Only (RUO)' and is not suitable for diagnostic or "
              "clinical applications. No rights can be derived from the contents of this report.")


def register_fonts(font_dir=None):
    """Use the real Nimbus Sans L if available, else Helvetica (metrically identical clone)."""
    if font_dir:
        reg = os.path.join(font_dir, "NimbusSansL-Regular.ttf")
        bold = os.path.join(font_dir, "NimbusSansL-Bold.ttf")
        if os.path.isfile(reg) and os.path.isfile(bold):
            pdfmetrics.registerFont(TTFont("NimbusSans", reg))
            pdfmetrics.registerFont(TTFont("NimbusSans-Bold", bold))
            return "NimbusSans", "NimbusSans-Bold", "Nimbus Sans L (report fonts)"
    return "Helvetica", "Helvetica-Bold", "Helvetica (metric clone of Nimbus Sans L)"


class Supplement(BaseDocTemplate):
    def __init__(self, path, sample, logo=None, font_dir=None, first_page=1, total_pages=None):
        super().__init__(path, pagesize=A4, leftMargin=MARGIN_LEFT, rightMargin=MARGIN_RIGHT,
                         topMargin=MARGIN_TOP, bottomMargin=MARGIN_BOTTOM, title=f"{sample} neoantigens")
        self.sample = sample
        self.logo = logo
        self.first_page = first_page
        self.total_pages = total_pages
        self.font, self.bold, self.font_note = register_fonts(font_dir)
        frame = Frame(MARGIN_LEFT, MARGIN_BOTTOM, CONTENT_W,
                      A4[1] - MARGIN_TOP - MARGIN_BOTTOM, id="body",
                      leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0)
        self.addPageTemplates([PageTemplate(id="supp", frames=[frame], onPage=self._furniture)])
        self.styles = self._styles()

    def _styles(self):
        s = {}
        s["chapter"] = ParagraphStyle("chapter", fontName=self.bold, fontSize=10,
                                      textColor=PALETTE_ORANGE, spaceAfter=10, leading=12)
        s["tableTitle"] = ParagraphStyle("tableTitle", fontName=self.bold, fontSize=8,
                                         textColor=PALETTE_ORANGE, spaceAfter=4, leading=10)
        s["content"] = ParagraphStyle("content", fontName=self.font, fontSize=7,
                                      textColor=PALETTE_DARK_GREY, leading=8.5)
        s["key"] = ParagraphStyle("key", fontName=self.font, fontSize=7,
                                  textColor=PALETTE_MID_GREY, leading=8.5)
        s["sub"] = ParagraphStyle("sub", fontName=self.font, fontSize=6,
                                  textColor=colors.black, leading=7.5, spaceBefore=2)
        s["note"] = ParagraphStyle("note", fontName=self.font, fontSize=6,
                                   textColor=PALETTE_MID_GREY, leading=7.5, spaceBefore=2)
        return s

    def _furniture(self, canv, doc):
        canv.saveState()
        top = A4[1]
        # orange-slice logo, top left
        if self.logo and os.path.isfile(self.logo):
            canv.drawImage(ImageReader(self.logo), MARGIN_LEFT, top - HEADER_H + 4,
                           width=42, height=42, mask="auto")
        # "ORANGE Report (Research Use Only)"
        x = MARGIN_LEFT + 52
        y = top - 34
        canv.setFont(self.bold, 11)
        canv.setFillColor(PALETTE_ORANGE)
        canv.drawString(x, y, "ORANGE")
        w = canv.stringWidth("ORANGE ", self.bold, 11)
        canv.setFillColor(PALETTE_DARK_GREY)
        canv.drawString(x + w, y, "Report (Research Use Only)")
        canv.setFont(self.font, 7)
        canv.setFillColor(PALETTE_MID_GREY)
        canv.drawString(x, y - 11, "Neoantigen supplement \u2014 wigits-clinical")
        # orange side panel, top right
        panel_w, panel_h = 210, 46
        px, py = A4[0] - MARGIN_RIGHT - panel_w, top - HEADER_H + 2
        canv.setFillColor(PALETTE_ORANGE)
        canv.rect(px, py, panel_w, panel_h, stroke=0, fill=1)
        canv.setFillColor(colors.white)
        canv.setFont(self.bold, 7)
        canv.drawString(px + 8, py + panel_h - 15, "SAMPLE")
        canv.setFont(self.bold, 10)
        canv.drawString(px + 8, py + panel_h - 30, self.sample)
        # footer: page number then disclaimer
        n = self.first_page + canv.getPageNumber() - 1
        label = f"{n}/{self.total_pages}" if self.total_pages else f"S{canv.getPageNumber()}"
        canv.setFont(self.bold, 7)
        canv.setFillColor(PALETTE_ORANGE)
        canv.drawString(MARGIN_LEFT, MARGIN_BOTTOM - 14, label)
        canv.setFont(self.font, 6)
        canv.setFillColor(PALETTE_MID_GREY)
        canv.drawString(MARGIN_LEFT + 26, MARGIN_BOTTOM - 14, DISCLAIMER)
        canv.restoreState()


def header_cells(doc, labels):
    return [Paragraph(f'<font color="#656A6C">{l.upper()}</font>', doc.styles["key"]) for l in labels]


def orange_table(doc, labels, rows, widths, shade_rows=()):
    """A table with ORANGE's geometry: 23 pt uppercase header, hairline rules, smoke-grey candidates."""
    data = [header_cells(doc, labels)] + rows
    t = Table(data, colWidths=widths, rowHeights=[HEADER_ROW_H] + [None] * len(rows),
              repeatRows=1, hAlign="LEFT")
    style = [
        ("FONTNAME", (0, 0), (-1, -1), doc.font),
        ("FONTSIZE", (0, 0), (-1, -1), 7),
        ("TEXTCOLOR", (0, 1), (-1, -1), PALETTE_DARK_GREY),
        ("VALIGN", (0, 0), (-1, 0), "BOTTOM"),
        ("VALIGN", (0, 1), (-1, -1), "MIDDLE"),
        ("LINEBELOW", (0, 1), (-1, -1), 0.25, PALETTE_MID_GREY),
        ("TOPPADDING", (0, 1), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 1), (-1, -1), 3),
        ("LEFTPADDING", (0, 0), (-1, -1), 0),
        ("RIGHTPADDING", (0, 0), (-1, -1), 6),
    ]
    for i in shade_rows:
        style.append(("BACKGROUND", (0, i + 1), (-1, i + 1), PALETTE_SMOKE_GREY))
    t.setStyle(TableStyle(style))
    return t


def key_value_table(doc, pairs, widths=(170, 365)):
    rows = [[Paragraph(k, doc.styles["key"]), Paragraph(str(v), doc.styles["content"])] for k, v in pairs]
    t = Table(rows, colWidths=widths, hAlign="LEFT")
    t.setStyle(TableStyle([
        ("FONTNAME", (0, 0), (-1, -1), doc.font), ("FONTSIZE", (0, 0), (-1, -1), 7),
        ("TOPPADDING", (0, 0), (-1, -1), 1), ("BOTTOMPADDING", (0, 0), (-1, -1), 1),
        ("LEFTPADDING", (0, 0), (-1, -1), 0), ("VALIGN", (0, 0), (-1, -1), "TOP"),
    ]))
    return t


def _num(v, nd=0, dash="-"):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return dash
    return f"{f:.{nd}f}"


def build(candidates_tsv, summary_json, outfile, sample, logo=None, figures=(), font_dir=None,
          first_page=1, total_pages=None, max_review=18):
    rows = list(csv.DictReader(open(candidates_tsv), delimiter="\t"))
    S = json.load(open(summary_json))
    doc = Supplement(outfile, sample, logo=logo, font_dir=font_dir,
                     first_page=first_page, total_pages=total_pages)
    st = doc.styles
    story = [Paragraph("Neoantigens", st["chapter"])]

    tiers = S["tiers"]
    rna = S["rna_states"]
    thr = S["thresholds"]
    story += [
        Paragraph("Neoantigen Summary", st["tableTitle"]),
        key_value_table(doc, [
            ("Coding candidates assessed", S["candidates"]),
            ("Prioritise / Review / Hold / Exclude",
             f"{tiers['Prioritise']} / {tiers['Review']} / {tiers['Hold']} / {tiers['Exclude']}"),
            ("Mutant allele in tumour RNA",
             f"{rna['SUPPORTED']} supported, {rna['COVERED_NO_ALT']} covered without the mutant allele, "
             f"{rna['NOT_COVERED']} not covered"),
            ("Peptides scored by both predictors", f"{S['peptides_scored_by_both']:,}"),
            ("Class I typing QC", S["lilac_qc_status"]),
            ("Genetic immune escape", S["case_immune_escape"] if isinstance(S["case_immune_escape"], str)
             else ", ".join(S["case_immune_escape"]))  ,
            ("Binding thresholds",
             f"{_num(thr['binding_strong_nm'])} nM and {_num(thr['binding_percentile'],1)} percentile"),
            ("Class II neoepitopes", "Not assessed \u2014 class II prediction is not performed"),
            ("Fusion and splice neoepitopes", "Not assessed from RNA \u2014 DNA fusions only, via NEO/LINX"),
            ("Rules version", S["rules_version"]),
        ]),
        Spacer(1, 20),
    ]

    def cand_rows(subset):
        out = []
        for r in subset:
            wt = _num(r["ic50_wt"])
            out.append([
                Paragraph(r["gene"], st["content"]),
                Paragraph(r["aa_change"], st["content"]),
                Paragraph(r["peptide"], st["content"]),
                Paragraph(r["allele"], st["content"]),
                Paragraph(f'{_num(r["ic50_mt"])} ({wt})', st["content"]),
                Paragraph(_num(r["percentile_mt"], 2), st["content"]),
                Paragraph({"SUPPORTED": _num(r["rna_vaf_sage"], 2),
                           "COVERED_NO_ALT": "no alt",
                           "NOT_COVERED": "no cov"}.get(r["rna_status"], "-"), st["content"]),
                Paragraph(_num(r["gene_tpm"], 1), st["content"]),
                Paragraph("Subclonal" if r["subclonal"] == "True" else "Clonal", st["content"]),
                Paragraph(r["neo_likelihood"] or "-", st["content"]),
            ])
        return out

    labels = ["Gene", "Variant", "Peptide", "Allele", "Binding mut (wt) nM", "%ile",
              "RNA", "TPM", "Clonality", "NEO"]
    # widths are tuned so no header word wraps mid-word at font size 7; CLONALITY is the binding one
    widths = [46, 40, 74, 56, 70, 26, 32, 30, 52, 46]

    pri = sorted([r for r in rows if r["tier"] == "Prioritise"], key=lambda r: float(r["ic50_mt"]))
    rev = sorted([r for r in rows if r["tier"] == "Review"], key=lambda r: float(r["ic50_mt"]))
    story += [
        Paragraph(f"Prioritised Candidates ({len(pri)})", st["tableTitle"]),
        orange_table(doc, labels, cand_rows(pri), widths),
        Paragraph("Strong predicted binding on both measures, mutant allele expressed in tumour RNA, "
                  "clonal, and the wild-type peptide binds less well. Binding is the median mutant IC50 "
                  "with the wild-type in brackets.", st["note"]),
        Spacer(1, 20),
    ]
    shown = rev[:max_review]
    story += [
        Paragraph(f"Candidates Under Review ({len(rev)})", st["tableTitle"]),
        orange_table(doc, labels, cand_rows(shown), widths,
                     shade_rows=range(len(shown))),
        Paragraph("Shaded as candidate-grade: a strong binding prediction with one of expression, "
                  "clonality or the mutant-over-wild-type differential missing. "
                  + (f"Strongest {len(shown)} of {len(rev)} shown; the full set is in the candidate table."
                     if len(shown) < len(rev) else "Human review decides these."), st["note"]),
    ]

    held = [r for r in rows if r["tier"] == "Hold"]
    excl = [r for r in rows if r["tier"] == "Exclude"]
    story += [
        Spacer(1, 16),
        Paragraph("Not Carried Forward", st["tableTitle"]),
        key_value_table(doc, [
            (f"Held ({len(held)})",
             "RNA never covered the position, so absence of the mutant transcript is not established, "
             "or no strong binding prediction. Held, not excluded."),
            (f"Excluded ({len(excl)})",
             "Ineligible on evidence that a binding score cannot override \u2014 variant alt reads in the "
             "matched normal, restricting allele lost, or typing QC failure."),
        ]),
    ]

    for fig in figures:
        if os.path.isfile(fig):
            iw = CONTENT_W
            from PIL import Image as PILImage
            w, h = PILImage.open(fig).size
            story += [Spacer(1, 18), Image(fig, width=iw, height=iw * h / w)]

    doc.build(story)
    return outfile, doc.font_note


def append_to_orange(orange_pdf, supplement_pdf, outfile):
    """Concatenate without touching the original: ORANGE pages first, byte-identical."""
    from pypdf import PdfReader, PdfWriter
    w = PdfWriter()
    for p in PdfReader(orange_pdf).pages:
        w.add_page(p)
    n_orange = len(w.pages)
    for p in PdfReader(supplement_pdf).pages:
        w.add_page(p)
    with open(outfile, "wb") as fh:
        w.write(fh)
    return outfile, n_orange, len(w.pages) - n_orange


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidates-tsv", required=True)
    ap.add_argument("--summary-json", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--logo", default=None, help="orange-slice PNG (extractable from any ORANGE report)")
    ap.add_argument("--font-dir", default=None, help="directory holding NimbusSansL-{Regular,Bold}.ttf")
    ap.add_argument("--figure", action="append", default=[], help="PNG to embed; repeatable")
    ap.add_argument("--orange-pdf", default=None, help="if given, also write <out>.extended.pdf")
    ap.add_argument("--first-page", type=int, default=1)
    ap.add_argument("--total-pages", type=int, default=None)
    a = ap.parse_args(argv)
    out, note = build(a.candidates_tsv, a.summary_json, a.out, a.sample, logo=a.logo,
                      figures=a.figure, font_dir=a.font_dir,
                      first_page=a.first_page, total_pages=a.total_pages)
    print(json.dumps({"supplement": out, "fonts": note}))
    if a.orange_pdf:
        ext = a.out.replace(".pdf", "") + ".extended.pdf"
        _, n_o, n_s = append_to_orange(a.orange_pdf, out, ext)
        print(json.dumps({"extended": ext, "orange_pages": n_o, "supplement_pages": n_s}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
