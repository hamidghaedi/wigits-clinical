# Figure and table style for the extended (neoantigen) section

Companion to `docs/orange_neoantigen_reporting.md`. That document says *where* the extended section goes
and how the page furniture is reproduced; this one says *how the content is drawn*, so that every plot
and table added by this layer looks like one family and like ORANGE.

Grounded in a short literature review (§1), with the palette contract in §2, the recommended panels in
§3 and what was deliberately rejected in §5.

## 1. What the field actually shows, and what that implies

**pVACview** (Xia et al., *Genome Medicine* 2024;16:132, doi:10.1186/s13073-024-01384-7) is the reference
point: the first interactive tool built specifically for neoantigen prioritisation, and the tool our own
review SOP already uses. Four things from it drive the design here.

1. **Three levels of detail — variant, transcript, peptide.** Candidate information is organised by level,
   not as one flat table, because a variant can sit in several transcripts that yield different peptides.
   *Implication:* the extended section should be layered the same way, and should never present a single
   "best" peptide as if the choice were obvious. The paper is explicit that static reports
   "over-simplify the outputs to make them tractable", leaving the results "fraught with assumptions
   about what 'best' means".
2. **Anchor position is first-class.** pVACview renders an allele- and length-specific anchor probability
   heatmap per peptide, with the mutated residue marked and mutant/wild-type pairs separated, because a
   mutation at an anchor changes *binding* whereas one at an exposed position changes *recognition*.
   *Implication:* a position map belongs in the report (§3, panel C). This is the single largest gap in
   the first version of our supplement.
3. **Failing criteria are shown, not just the verdict.** The tool "highlights the specific failing
   criteria" so a candidate can be rescued when a case has too few candidates.
   *Implication:* an evidence matrix (§3, panel D) is more useful than a tier label alone.
4. **Manufacturability and self-matching are prioritisation criteria too** — problematic residues
   (cysteine is their example, for synthesis and stability) and reference-proteome matches.
   *Implication:* reserve columns for both even though this version does not compute them.

Their Table 1 also states the thresholds we already use — binding affinity ≤ 500 nM and percentile ≤ 2 —
and recommends prioritising transcripts with `TSL = 1` and a `protein_coding` biotype, and genes in the
Cancer Gene Census.

**pVACtools v6** (arXiv 2606.26659) is now published, which is newer than the 5.4.0 we have pinned; it
adds safety evaluation and vaccine-manufacturing support. Worth a look before the next release, but not
a reason to move mid-development.

Note on provenance: the PubMed connector required an approval that was not granted in this session, so the
above comes from web search plus the pVACview full text fetched directly via its DOI. A broader review —
NeoFox, vaxrank, LENS, NeoPredViz, TESLA — is still outstanding.

## 2. Palette contract

**Every colour used by this layer already exists in ORANGE `ReportResources`.** No new hue is introduced.
Three roles, three encodings, and they must not be mixed:

| Role | Encoding | Values |
|---|---|---|
| **Tier** (unordered, one focal) | hue | Prioritise `PALETTE_ORANGE (242,139,31)` · Review `PALETTE_BLUE (38,90,166)` · Hold `PALETTE_MID_GREY (101,106,108)` · Exclude `PALETTE_LIGHT_GREY (211,211,211)` |
| **Ordered quantity** | single-hue ramp | `ORANGE_1..6`, light `(255,165,0)` to dark `(155,115,0)` — expression deciles, anchor probability, cohort percentile |
| **Evidence state** | fill, not hue | solid = mutant allele seen in RNA · hollow = covered, mutant not seen · hatched = not covered · `SMOKE_GREY` row shading = candidate-grade · `GAINSBORO_GREY` = grid or absent |

Rules that follow from this:

- Tier colour and evidence fill are **orthogonal channels**: a reader takes in both at once, and no
  candidate needs two legend lookups.
- The `ORANGE_1..6` ladder is **only** for ordered quantities. Using it for categories wastes the one
  ordinal encoding available and reads as a false ranking.
- **Never red or green** — neither is in the ORANGE palette, and a red/green opposition fails in
  deuteranopia. Where pVACview marks a mutated residue in red, we use `PALETTE_ORANGE`; where it marks a
  driver gene in a green box, we use bold text plus the census column.
- Orange against blue is the one strong categorical contrast in the palette and survives colour-vision
  deficiency; it is reserved for the Prioritise/Review distinction, which is the one that decides
  reviewer effort.
- White with a mid-grey outline means **not assessed** and never zero. A hatched marker means the
  measurement was not made. Absence of evidence gets its own ink.

Typography follows the report: Nimbus Sans L (Helvetica is its metric clone), 8 pt panel titles, 7 pt
axis and body, 6 pt annotation and legend; axes and spines in `PALETTE_MID_GREY` at 0.25 pt, matching the
table cell rules exactly; 300 dpi; content width 535 pt.

## 3. Recommended panels

Each is a mock rendered from the real KQ0004 data, not a sketch.

| Panel | Shows | Why it earns the space |
|---|---|---|
| **A — tier-ordered strip** | all 69 candidates, binding percentile, tier band, fill = RNA evidence | the only view that shows the whole shortlist and its shape; makes the Hold block's size obvious |
| **B — mutant versus wild-type** | arrow from wild-type to mutant IC50 on a log axis, per candidate | makes differential agretopicity visible, which a single IC50 column cannot |
| **C — anchor / position map** | peptide grid, mutated residue in orange, canonical anchors shaded | answers "does this mutation change binding or recognition?" — pVACview's central question |
| **D — evidence matrix** | candidate × criterion, met / not met / not assessed | shows *why* a tier was assigned and which single criterion is blocking a rescue |
| **E — allele coverage** | candidates per restricting allele, stacked by tier | mirrors the native Immunology chapter's allele-centric organisation |

Panel B produced the most consequential finding on this case, and it becomes visible only when the two
affinities are drawn together. Of the six prioritised candidates, **three show a ratio below two**: *GSX1*
at 1.15, *PPM1F* at 1.22 and *ZSWIM9* at 1.60, where the ratio is the wild-type affinity divided by the
mutant affinity. The remaining three are *MAPK13* at 43.3, *RCOR2* at 35.2 and *PAPLN* at 18.2. A ratio
near one means the two peptides bind comparably well, so the substitution gives a T cell little to
distinguish, and a ratio below one means the wild-type peptide is the stronger binder. Eight candidates
under review fall below one, among them *TP53* F270C at 0.51.

The rule in force until now asked only that the wild-type peptide bind less well, which a ratio of 1.15
satisfies. A minimum ratio has therefore been added to the tier rules, set at two and recorded with every
candidate as `fold_change`. Published work does not settle on a single cut-off, so the value is stated in
the summary block of every report and can be revised.

Panel D surfaced a second gap: *PAPLN* T825M is Prioritise with `TSL = 2`. Transcript support level is in
pVACview's prioritisation criteria but is not in our tier rules — the pVACseq TSL filter runs upstream and
then the information is dropped. Either add it to the rules or mark the column informational.

## 4. Tables

Unchanged from the native geometry (23 pt upper-cased mid-grey headers, 0.25 pt rules, smoke-grey for
candidate-grade rows). Three additions:

- a **Cancer Gene Census** marker as bold gene text, not a coloured box;
- **mutant and wild-type binding in one cell** — `205 (3725)` — so the differential is readable in the
  table as well as the figure;
- reserved, explicitly empty columns for **problematic residues** and **reference-proteome match**, shown
  as `not assessed` rather than omitted.

## 5. Rejected

- **Violin plots of feature distributions** (pVACview's NeoFox panel). With 69 candidates on one case the
  distribution is the data; a strip plot shows it honestly and a violin implies a density that is not there.
- **A dynamic scatter with user-chosen axes.** That is what pVACview is for. A static report should make
  one defensible choice and name it.
- **A radar or spider plot of candidate features.** Area is not a meaningful quantity across
  non-commensurable axes, and the shape changes with axis order.
- **A single composite score.** The repository's standing decision is interpretable tiers with stated
  reasons rather than one opaque number, and nothing in the review changes that.
