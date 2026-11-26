# Neoantigen candidate review — standard operating procedure

Applies to the neoantigen discovery layer (gap G28) running downstream of nf-core/oncoanalyser 3.0.0 on
a tumour/normal WGS + tumour RNA case. Version 1, written against KQ0004.

**Scope.** This procedure produces a reviewed, reproducible *research* shortlist of candidate neoantigens.
It is not a vaccine design, not a manufacturing specification and not a clinical recommendation. No
candidate leaves this process without a named reviewer and a recorded decision.

## 1. What the reviewer is given

Each case release contains a review bundle:

```
wigits_releases/<run_id>/neo_review/
├── <case>.neo_candidates.tsv              tiered candidates, one row per (variant, allele, peptide)
├── <case>.neo_candidates.summary.json     counts, thresholds, rules version, input checksums
├── <case>.neo_concordance.tsv             pVACseq versus Hartwig NEO on peptides both scored
├── <case>.interface_manifest.json         sha256 of every oncoanalyser input consumed
├── pvacview/
│   ├── <case>.all_epitopes.aggregated.tsv       pVACview input
│   └── <case>.all_epitopes.aggregated.metrics.json
├── igv/
│   ├── <case>.review_sites.bed            one interval per candidate under review
│   └── igv_session.xml                    tumour DNA, normal DNA, tumour RNA, the selected VCF
├── decisions/
│   └── <case>.review_decisions.tsv        WRITTEN BY THE REVIEWER (template provided)
└── README.md                              case facts: purity, ploidy, TMB, MSI, HLA, immune escape
```

The bundle is self-contained: everything needed to re-reach a decision is inside it, and the manifest
ties it to exact oncoanalyser outputs.

## 2. Before looking at any candidate — case-level gates

Read `README.md` and stop if any of these holds; the shortlist is not interpretable otherwise.

| Gate | Source | Why it stops review |
|---|---|---|
| LILAC QC not `PASS` | `lilac.qc.tsv` | the restricting alleles themselves are uncertain |
| Tumour purity below ~0.20 | PURPLE | subclonality and VAF become uninformative |
| RNA library failed or very low depth | ISOFOX / MultiQC | `NOT_COVERED` would dominate and mean nothing |
| An immune-escape flag is set | ORANGE `immuneEscape` | B2M / HLA / IFN-γ / CD58 / SETDB1 / PD-L1 escape changes how presentation evidence should be read — review continues, but every decision is annotated |
| MSI-high | PURPLE | indel neoantigen load is large; the shortlist is then a sample, not a ranking |

Record which gates were checked, with their values, at the top of the decisions file.

## 3. Tier triage

Tiers come from `workflow/neo_candidates.py` (`--explain` prints the rules in force, and every row carries
`tier_reasons`). The reviewer's effort is allocated accordingly:

- **Prioritise** — review every row in full (sections 4 and 5).
- **Review** — review every row in full; these are the ones where a single piece of evidence is missing
  and human judgement actually decides.
- **Hold** — read the reason only. Promote to full review when the reason is `RNA not covered` *and* the
  variant is in a gene of interest; otherwise leave held.
- **Exclude** — do not review. Spot-check two per case to confirm the exclusion logic behaved.

Never promote a candidate past its eligibility reason (normal alt reads, allele lost, QC fail). Those are
not judgement calls.

## 4. pVACview session

```
pvactools pvacview run <bundle>/pvacview
```

Open the aggregated report with its metrics JSON and, for each candidate under review, check:

1. **Transcript set** — does the peptide arise from a transcript that is actually expressed, or only from a
   minor isoform? The per-transcript panel shows this; a candidate supported only by an unexpressed
   isoform drops to Hold.
2. **Peptide-level evidence** — is the mutant epitope strong across *several* lengths and registers, or a
   single 9-mer outlier? Breadth across registers is the more reliable signal.
3. **Mutant versus wild-type** — the differential agretopicity matters more than the absolute IC50. A
   mutant that binds at 300 nM while the wild-type binds at 350 nM is not a neoantigen worth pursuing.
4. **Predictor agreement** — MHCflurry, MHCflurryEL and MHCnuggetsI are independent enough that unanimity
   is informative. Disagreement is a reason to look harder, not a reason to drop.
5. **Anchor positions** — a mutation at a canonical anchor residue changes binding; a mutation at a
   solvent-exposed position changes recognition. pVACview's anchor view distinguishes them; record which.

Mark each candidate `Accept`, `Reject` or `Pending` in pVACview and export the session; the exported
decisions are then transcribed into the decisions file with the free-text reason.

## 5. IGV confirmation

Every candidate that survives section 4 is looked at in IGV. This is where false positives that no filter
catches are found.

Load the provided session (tumour DNA, normal DNA, tumour RNA, the selected VCF) against the **HMF masked
reference** (`GRCh38_masked_exclusions_alts_hlas.fasta`) — the same reference the calls were made on, never
GRCh38.d1.vd1.

Confirm, for each site:

- supporting reads are **not** all on one strand, and not all starting at the same coordinate (both are
  signatures of PCR or alignment artefact rather than a real variant);
- the alt is absent from the normal at adequate depth — eyeballing it catches low-level contamination that
  the VAF filter passes;
- the region is not a tandem repeat, homopolymer or segmental duplication; indel candidates in
  homopolymers are rejected by default;
- the RNA track agrees with the recorded `rna_status`: for `SUPPORTED`, mutant reads are visible and
  spliced consistently with the claimed transcript; for `COVERED_NO_ALT`, coverage is genuinely present;
- for FFPE material, C>T / G>A calls with few unique fragment start coordinates are treated as deamination
  until shown otherwise.

Record the IGV verdict per candidate, including the rejections — a rejected candidate with its reason is
the most useful record this process produces.

## 6. The decisions file

`decisions/<case>.review_decisions.tsv` — written by the reviewer, one row per reviewed candidate:

| Column | Content |
|---|---|
| `case`, `variant`, `gene`, `peptide`, `allele` | copied from the candidate table (the join key) |
| `assigned_tier` | the tier the pipeline produced |
| `reviewer` | name or identifier |
| `review_date` | ISO date |
| `pvacview_call` | Accept / Reject / Pending |
| `igv_call` | Confirmed / Artefact / Ambiguous / Not reviewed |
| `final_call` | Carry forward / Drop / Defer |
| `reason` | free text; **required** whenever `final_call` disagrees with `assigned_tier` |
| `notes` | anchor position, isoform caveat, escape context, anything a second reader would need |

Two reviewers for anything intended to go further than an internal report; disagreements are resolved by
discussion and the discussion is recorded, not overwritten.

## 7. Reproducibility requirements

- The decisions file references `rules_version` from the summary JSON. If the rules change, prior decisions
  are **not** silently re-tiered: the case is re-run and re-reviewed, and both versions are kept.
- Tool versions come from `VERSIONS.md` in the release; predictions made with different predictor versions
  are not pooled.
- The `interface_manifest.json` checksums identify the exact oncoanalyser outputs consumed. A re-run whose
  checksums differ is a different analysis.

## 8. Known limits of this version

Stated here so a reviewer does not read absence as evidence.

- **Class I only.** Class II prediction is not performed; LILAC types DP/DQ natively but not DRB1/3/4/5, so
  class II coverage would be partial even if prediction were added.
- **SNV and indel only.** Fusion-derived and splice-derived neoepitopes from RNA (pVACfuse, pVACsplice) are
  not covered; DNA fusion neoepitopes come from Hartwig NEO/LINX only.
- **Licence-free predictors only.** The NetMHCpan family requires an IEDB standalone licence that is not
  held here, so the most widely benchmarked predictor is absent from the ensemble.
- **No proteasomal cleavage or TAP transport modelling**, and no TCR-recognition modelling beyond NEO's
  similarity scores.
- **`NOT_COVERED` is not evidence of absence.** It means the RNA did not interrogate the position.
- Clonality uses PURPLE `SUBCL` and variant copy number, not the pVACseq VAF heuristic, which assumes a
  clonal VAF of 0.5 and is wrong at this purity.
