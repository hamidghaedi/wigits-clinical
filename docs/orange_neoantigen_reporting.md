# Extending the ORANGE report with the neoantigen layer

**Status: design, not yet implemented.** Gap G28 (neoantigen layer) and G10 (immune-escape table).
Written from the ORANGE 5.0.1 report sources read on 2026-11-23 at tag `orange-v5.0.1`
(`orange/src/main/java/com/hartwig/hmftools/orange/report/`), the installed container
`quay.io-biocontainers-hmftools-orange-5.0.1--hdfd78af_0.img`, and the repository's existing
front-page overlay. Line references are to that tag.

The aim is that an added chapter is **indistinguishable in style from a native one**: same fonts, same
palette, same type sizes, same table geometry, same page furniture. Nothing is re-styled, nothing is
invented; the added content borrows the native resources.

## 1. How ORANGE renders, and where an addition goes

`ReportWriter.writePdf` builds a fixed chapter list and hands each to `writePdfChapters`:

```
FrontPage → SomaticFindings → GermlineFindings (not tumour-only) → Immunology
          → RnaFindings (if isofox) → Cuppa (not tumour-only) → QualityControl (if qSee plot)
          → PurplePlots
```

Each chapter implements `ReportChapter` (`name()`, `pageSize()`, `render(Document)`), gets its own
`AreaBreak`, and the page header/footer are drawn by `PageEventHandler` from the chapter name. So an added
chapter needs no layout work of its own: give it a name and render tables into the document.

**Insertion point.** A neoantigen chapter belongs immediately after `Immunology` — it is the same
biology (HLA presentation) and the native Immunology chapter already loads LILAC. Two options:

| Option | What it means | Cost |
|---|---|---|
| **A. Extend `ImmunologyChapter`** | add `addNeoantigenData(document)` after `addHLAData(document)`; also un-comment `addImmuneEscapeData` (G10) | one overlaid class, exactly the pattern `FrontPageChapter.java` already uses |
| B. New `NeoantigenChapter` | own chapter, own TOC entry, needs `ReportWriter` overlaid too | two overlaid classes, more divergence from upstream |

**Recommend A.** The repository already compiles one overlaid class against the installed JAR and puts it
first on the classpath (`modules/local/orange/main.nf`, `assets/orange_qc_report/FrontPageChapter.java`);
the same mechanism takes `ImmunologyChapter.java` with no further plumbing. `javac` is present in the
ORANGE container (`/usr/local/bin/javac`, OpenJDK 21), which is what makes this work.

**G10 comes free.** In ORANGE 5.0.1 the immune-escape block is *present but commented out* in
`ImmunologyChapter` (the `addImmuneEscapeData` method and its call site), and `ImmuneEscapeTable` is a
complete class. Un-commenting it in the same overlay renders HLA / antigen-presentation / IFN-γ / CD58 /
PD-L1 / SETDB1 escape, which the JSON already carries (verified on KQ0004: all six flags present, all
false). That closes G10 with no new rendering code.

## 2. The native style, as constants (`ReportResources`)

Every value below is read from the 5.0.1 source; an added table must use these and nothing else.

**Fonts** — Nimbus Sans L, loaded from classpath resources by `ReportResources.create()`:
`fonts/nimbus-sans/NimbusSansL-Regular.ttf` and `…-Bold.ttf`, encoded `IDENTITY_H`. An overlay must
obtain fonts through the injected `ReportResources` instance, never load its own.

| Style | Font | Size | Colour |
|---|---|---|---|
| `chapterTitleStyle()` | bold | 10 | `PALETTE_ORANGE` |
| `tableTitleStyle()` | bold | 8 | `PALETTE_ORANGE` |
| `tableHeaderStyle()` | regular | 7 | `PALETTE_MID_GREY` |
| `tableContentStyle()` | regular | 7 | `PALETTE_DARK_GREY` |
| `keyStyle()` / `valueStyle()` | regular | 7 | `PALETTE_MID_GREY` |
| `subTextStyle()` | regular | 6 | `PALETTE_BLACK` |
| `disclaimerStyle()` | regular | 6 | `PALETTE_MID_GREY` |
| `qcWarningStyle()` | bold | 7 | `PALETTE_DARK_GREY` |
| `urlStyle()` | regular | 7 | `PALETTE_BLUE` |
| `pageNumberStyle()` | bold | 7 | `PALETTE_ORANGE` |

**Palette** (`DeviceRgb`): `ORANGE (242,139,31)`, `ORANGE_1..6 (255,165,0) → (155,115,0)`,
`BLUE (38,90,166)`, `DARK_GREY (39,47,50)`, `MID_GREY (101,106,108)`, `LIGHT_GREY (211,211,211)`,
`GAINSBORO_GREY (220,220,220)`, `SMOKE_GREY (240,240,240)`, `WHITE`, `BLACK`.

**Page**: A4; margins top 100 (the chapter title is drawn in the header), right 30, bottom 40, left 30;
header band 60; footer 20; full-page image 750 × 430; front-page circos 400.

**Table geometry** (`util/Tables`, `util/Cells`):
- content tables are `UnitValue.createPercentArray(widths).setWidth(contentWidth())` with `setFixedLayout()`;
- bottom margin 20 between tables;
- header cells: fixed height 23 pt, **text upper-cased**, no border, bottom-aligned, `tableHeaderStyle`;
- content cells: no border except `SolidBorder(PALETTE_MID_GREY, 0.25F)` on the bottom, `tableContentStyle`;
- long tables use `Tables.createWrapping(...)`, which adds "THE TABLE CONTINUES ON THE NEXT PAGE" /
  "CONTINUED FROM THE PREVIOUS PAGE" in `subTextStyle`;
- an empty or unavailable table is `Tables.createEmpty(title, width)` → `NONE`, or `createNotAvailable` → `NA`,
  never a blank space;
- `shadeCandidateCells` paints `SMOKE_GREY` behind rows that are candidates rather than confirmed —
  **the native idiom for "lower confidence", and the right treatment for our Review/Hold tiers**.

## 3. Proposed tables

### 3.1 `Neoantigen Candidates` (one row per prioritised candidate)

Columns, chosen to mirror the native driver tables' economy (gene first, evidence last) and to fit A4
content width at font size 7:

| Column | Source | Note |
|---|---|---|
| Gene | candidate table `gene` | CGC genes are not bolded — the native tables do not bold within a row |
| Variant | `aa_change` | e.g. `S156I` |
| Peptide | `peptide` | mutant epitope |
| Allele | `allele` | `HLA-C*12:03` |
| Binding (nM) | `ic50_mt` | mutant IC50; wild-type shown in brackets when it differs by < 2-fold |
| %ile | `percentile_mt` | |
| RNA | `rna_status` + `rna_vaf` | `SUPPORTED 0.45`, `COVERED_NO_ALT`, `NOT_COVERED` — never a bare 0 |
| Clonality | `subclonal_likelihood` → `Clonal` / `Subclonal` | from PURPLE `SUBCL`, not VAF |
| NEO likelihood | `neo_likelihood` | Hartwig's own presentation score, for comparison |
| Tier | `tier` | Prioritise / Review / Hold |

Rendering rules: `Prioritise` rows plain; `Review` and `Hold` rows shaded with `shadeCandidateCells`
(the native candidate idiom); `Exclude` rows are **not rendered** but counted in the caption. Table title
`tableTitleStyle`, so it reads as a sibling of "HLA Class I Alleles".

### 3.2 `Neoantigen Summary` (key/value, `createKey`/`createValue`)

Coding candidates assessed · peptides predicted · Prioritise/Review/Hold/Exclude counts ·
restricting alleles retained (from LILAC) · class I typing QC · predictors and versions ·
"not assessed" rows for class II and RNA-only fusion/splice candidates, which this version does not cover.

### 3.3 `Genetic Immune Escape` (G10)

The native `ImmuneEscapeTable` as written upstream, un-commented. No new code.

## 4. Figures, if any

ORANGE embeds pre-rendered PNGs (`util/Images`, `PlotPathResolver`); it draws nothing itself. So a
neoantigen figure must be produced upstream as a PNG and passed in. To look native it must match the
report, not a matplotlib default:

- figure text in **Nimbus Sans** (available on the cluster as a Ghostscript/URW font; otherwise the
  metrically-identical Helvetica) at **7 pt body / 8 pt title**, matching `tableContent`/`tableTitle`;
- colours **only** from the palette above — `PALETTE_ORANGE` for the highlighted series,
  `ORANGE_1..6` for ordered categories (that is what the native ladder exists for),
  `MID_GREY` for axes and ticks, `SMOKE_GREY` for fills, `BLUE` only for links;
- no chart junk: no grid unless it carries reading precision, no frame on legends, ticks outward,
  spines in `MID_GREY` at 0.25 pt to match the cell borders;
- sized to the native image boxes (full page 750 × 430 pt, i.e. 10.4 × 6.0 in at 72 dpi) and exported at
  300 dpi so the embedded raster is sharp at A4.

A single candidate figure is proposed: **binding strength versus RNA support**, IC50 percentile on one
axis and mutant-allele RNA VAF on the other, points coloured by tier and sized by clonality, with the
`NOT_COVERED` candidates drawn on a separate "no coverage" stripe rather than at zero — the one plot that
makes the three-state RNA rule visible. Everything else in the chapter is better as a table.

Implementation note: the repository's figure work should load the `figure-style` skill and then override
the font and palette to the values above; the house style is a starting point, the ORANGE palette is the
contract.

## 5. Data path into ORANGE

ORANGE reads `OrangeRecord` from its JSON/config inputs; it has no plug-in slot for an arbitrary table.
Two honest options:

1. **Overlay reads our TSV directly** (chosen for the first version): the ORANGE task is given
   `-neo_candidates_tsv <path>` through an extra argument in `modules/local/orange/main.nf`, and the
   overlaid `ImmunologyChapter` reads it with the same CSV reader hmftools already uses. Minimal, and the
   file is exactly the artefact the review step signs off.
2. Extend the `OrangeRecord` data model (`hmftools-datamodel`) with a neoepitope section. Cleaner, but it
   changes a published schema and every consumer of the JSON — out of scope until the layer is settled.

Whichever path, the **JSON must gain the same content as the PDF** (the repository's own rule that the
JSON is the machine-readable record): a `neoantigens` block with the candidate rows and the summary
counts, written by our own step rather than by ORANGE, and referenced from the release manifest.

## 6. What this must not do

- Not re-fit, re-call or re-interpret anything ORANGE already reports.
- Not present a computational shortlist as a clinical recommendation: the chapter carries the same
  research-use framing as the rest of the report, and the review SOP's decisions are what matter.
- Not silently drop candidates: excluded ones are counted in the caption with the reason class.
- Not diverge from upstream any further than the one overlaid class the mechanism already requires;
  every added line is traceable to this document.

## 7. Sources

- `ReportResources.java`, `ReportWriter.java`, `chapters/ImmunologyChapter.java`, `util/Cells.java`,
  `util/Tables.java`, `util/Images.java`, `tables/ImmuneEscapeTable.java` at `orange-v5.0.1`.
- This repository: `docs/ORANGE_QC_REPORT_RECOVERY.md`, `workflow/patches/oncoanalyser-3.0.0-qc-report.patch`,
  `workflow/overlays/oncoanalyser-3.0.0/assets/orange_qc_report/FrontPageChapter.java`.
- Verified on the KQ0004 ORANGE JSON: `immuneEscape` carries all six flags; `lilac.alleles` carries
  class I and class II.
