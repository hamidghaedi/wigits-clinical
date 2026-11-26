# Extending the ORANGE report with the neoantigen layer

**Approach: ORANGE is left untouched and the neoantigen content is rendered as additional pages in
ORANGE's own visual language, appended after the original report.** Gap G28 (neoantigen layer);
gap G10 (immune-escape table) is addressed separately below.

Written from the ORANGE 5.0.1 report sources read at tag `orange-v5.0.1`
(`orange/src/main/java/com/hartwig/hmftools/orange/report/`) **and from a rendered report** —
`KQ0004-2023-AAAD-2-DNA.orange.pdf`, 11 pages, 13 MB, produced by the native oncoanalyser 3.0.0 run.
Reading the source alone was not enough: the page furniture, the orange-slice logo, the side panel and
the footer disclaimer are only apparent in the rendered document.

## 1. Why a supplement rather than a patched chapter

An earlier draft of this document proposed overlaying `ImmunologyChapter.java` — compiling one modified
class against the installed JAR and putting it first on the classpath, the pattern
`workflow/patches/oncoanalyser-3.0.0-qc-report.patch` already uses for the front page. That remains
technically viable (`javac` and OpenJDK 21 are present inside the ORANGE container), but it is the wrong
default here:

| | Patched chapter | **Appended supplement** |
|---|---|---|
| ORANGE output | modified | **byte-identical** (verified by sha256 before and after) |
| Breaks on an ORANGE upgrade | yes — the overlay is compiled against one JAR | no |
| Needs the container rebuilt or classpath-ordered | yes | no |
| Review bundle can be regenerated without re-running ORANGE | no | yes |
| Content limited to what `OrangeRecord` carries | effectively yes | no |

The supplement costs one thing in exchange: the appended pages are not in ORANGE's table of contents and
are not in the ORANGE JSON. The second point matters and is handled in §5.

## 2. The real document, as rendered

| Page | Chapter |
|---|---|
| 1 | Front page — primary tumour, purity, ploidy, fit method, QC; driver summary; genome-wide biomarkers; CUPPA call |
| 2 | Somatic Findings — small variants, amplifications/deletions, fusions, disruptions, viruses |
| 3 | Arm copy-number aberrations; signature allocations; structural driver plots |
| 4–5 | LINX structural driver plots and the whole-genome circos (landscape) |
| 6 | Germline Findings, including pharmacogenetics |
| **7** | **Immunology — HLA class I alleles, then HLA class II alleles** |
| 8 | RNA Findings — RNA QC, high gene expression, RNA-only fusions, novel splice junctions |
| 9 | Tissue of Origin (CUPPA) |
| 10 | Quality Control |
| 11 | Purity and Ploidy |

Page 7 is the chapter the neoantigen content is a continuation of, and it is **half empty** — the two
allele tables occupy the top third. Class II *is* rendered for a whole-genome experiment (DP/DQ, since
LILAC does not type DRB1/3/4/5), with columns ALLELE, QC STATUS, REF FRAGS, TUMOR FRAGS, RNA FRAGS,
TUMOR CN, SOMATIC MUTATIONS. The supplement's allele vocabulary matches that table exactly.

**No immune-escape table is rendered**, confirming G10 from the output rather than only from the source.

## 3. Page furniture to reproduce

From the rendered page, top to bottom:

- **orange-slice logo**, top left, roughly 42 pt square. Extracted from any ORANGE PDF rather than
  re-drawn (`pypdf` → `page.images`), so it is the identical asset, not an imitation.
- **title** — `ORANGE` in bold `PALETTE_ORANGE`, then `Report (Research Use Only)` in bold dark grey.
- **side panel** — a filled `PALETTE_ORANGE` rectangle at the top right carrying `SAMPLE` (bold 7, white)
  over the sample name (bold 10, white).
- **chapter title** — bold 10 `PALETTE_ORANGE`, at the top of the text frame.
- **footer** — page number `n/total` in bold 7 orange, then the RUO disclaimer in regular 6 mid-grey:
  *"This report is for 'Research Use Only (RUO)' and is not suitable for diagnostic or clinical
  applications. No rights can be derived from the contents of this report."*

Page geometry: A4; margins top 100, left 30, right 30, bottom 40; header band 60; footer 20; content
width 535 pt.

## 4. Type, palette and table geometry

Every value below is from `ReportResources.java` and is reproduced verbatim by the renderer.

| Style | Font | Size | Colour |
|---|---|---|---|
| chapter title | bold | 10 | `PALETTE_ORANGE` |
| table title | bold | 8 | `PALETTE_ORANGE` |
| table header | regular | 7 | `PALETTE_MID_GREY`, **upper-cased**, fixed 23 pt row |
| table content | regular | 7 | `PALETTE_DARK_GREY`, 0.25 pt mid-grey rule below |
| key / value | regular | 7 | `PALETTE_MID_GREY` |
| sub-text | regular | 6 | black |
| disclaimer | regular | 6 | `PALETTE_MID_GREY` |

Palette: `ORANGE (242,139,31)`, `ORANGE_1..6 (255,165,0) → (155,115,0)`, `BLUE (38,90,166)`,
`DARK_GREY (39,47,50)`, `MID_GREY (101,106,108)`, `LIGHT_GREY (211,211,211)`,
`GAINSBORO_GREY (220,220,220)`, `SMOKE_GREY (240,240,240)`.

**Fonts.** ORANGE embeds Nimbus Sans L (`fonts/nimbus-sans/NimbusSansL-{Regular,Bold}.ttf`). Nimbus Sans L
is URW's clone of Helvetica and is **metrically identical** to it, so reportlab's built-in Helvetica
reproduces the same advance widths, line breaks and column fits. `--font-dir` takes the real TTFs when
they are available; the renderer reports which was used.

**Shading.** `SMOKE_GREY` behind a row is ORANGE's own idiom for "candidate, not confirmed"
(`candidateDriverBackground`). The supplement uses it for the Review tier and nothing else, so a reader
who knows the native report already knows what the shading means.

## 5. What the supplement contains

**Page 1 — Neoantigens**

1. *Neoantigen Summary* (key/value): candidates assessed; the four tier counts; the three RNA states
   spelled out; peptides scored by both predictors; class I typing QC; immune escape; binding
   thresholds; and explicit **"Not assessed"** rows for class II and for RNA-derived fusion/splice
   neoepitopes, so a reader cannot mistake absence for a negative result.
2. *Prioritised Candidates*: gene, variant, peptide, allele, mutant IC50 with the wild-type in brackets,
   percentile, RNA state or mutant allele fraction, gene TPM, clonality, NEO likelihood.
3. *Candidates Under Review*: the same columns, smoke-grey shaded, strongest first.
4. *Not Carried Forward*: the Hold and Exclude counts with their reason classes — never a silent drop.

**Page 2 — figures**, rendered in matplotlib to the same contract (Nimbus Sans/Helvetica at 7 pt body and
8 pt title, palette colours only, mid-grey 0.25 pt spines matching the cell rules, 300 dpi, content
width):

- *binding strength versus RNA support* — percentile against mutant allele fraction, coloured by tier,
  hollow markers where the mutant allele was not seen, and a separate stripe for positions RNA never
  covered, so the three-state rule is visible rather than collapsed to zero;
- *candidates per restricting allele, stacked by tier* — deliberately allele-centric, because that is how
  the native Immunology chapter is organised.

**The JSON.** The appended pages are not in the ORANGE JSON, so the layer writes its own machine-readable
record — `<case>.neo_candidates.summary.json` plus the full candidate table — and the release manifest
references both alongside `<sample>.orange.json`. The rule that the machine-readable record carries the
same content as the PDF is preserved; it is simply carried by a sibling file.

## 6. Gap G10 — the immune-escape table

In ORANGE 5.0.1 the immune-escape block is **present but commented out** in `ImmunologyChapter`
(the `addImmuneEscapeData` method and its call site), while `ImmuneEscapeTable` is a complete class. The
data is in the JSON — verified on KQ0004, all six flags present and all false. Two routes:

1. render the six flags in the supplement, from the JSON, with no Java involved (consistent with §1);
2. un-comment the native block in a single-class overlay, which renders Hartwig's own table.

Route 1 is recommended for the same reasons as §1; route 2 stays documented because it is the only way to
get the *native* table, should that ever be wanted.

## 7. Implementation

`workflow/orange_supplement.py` — reads the candidate TSV and summary JSON, renders the pages, and
optionally concatenates: ORANGE pages first, then the supplement, into
`<sample>.orange.extended.pdf`. The original PDF is opened read-only and its sha256 is unchanged by the
operation. Page numbering of the appended pages continues the original's (`12/13`, `13/13`) when
`--total-pages` is given, and falls back to `S1`, `S2` for a standalone supplement.

```
python3 workflow/orange_supplement.py \
  --candidates-tsv <case>.neo_candidates.tsv \
  --summary-json   <case>.neo_candidates.summary.json \
  --sample  <tumour_sample_id> \
  --logo    <logo extracted from the ORANGE pdf> \
  --figure  fig_binding_rna.png --figure fig_alleles.png \
  --orange-pdf <sample>.orange.pdf \
  --first-page 12 --total-pages 13 \
  --out <case>.neoantigens.supplement.pdf
```

## 8. What this must not do

- Not modify, re-render or re-interpret anything ORANGE produces.
- Not present a computational shortlist as a clinical recommendation: the appended pages inherit the same
  RUO disclaimer, and the review SOP's recorded decisions are what carry weight.
- Not silently drop candidates — excluded ones are counted with their reason class.
- Not let visual harmony imply equal provenance: everything on these pages comes from the neoantigen
  layer, not from Hartwig's validated pipeline, and the chapter sub-title says so
  ("Neoantigen supplement — wigits-clinical").

## 9. Sources

- `ReportResources.java`, `ReportWriter.java`, `chapters/ImmunologyChapter.java`, `util/Cells.java`,
  `util/Tables.java`, `util/Images.java`, `tables/ImmuneEscapeTable.java` at `orange-v5.0.1`.
- The rendered report `KQ0004-2023-AAAD-2-DNA.orange.pdf` (oncoanalyser 3.0.0, 2026-11-15).
- This repository: `docs/ORANGE_QC_REPORT_RECOVERY.md`,
  `workflow/patches/oncoanalyser-3.0.0-qc-report.patch`.
