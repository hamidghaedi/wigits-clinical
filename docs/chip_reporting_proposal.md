# Reporting Clonal Haematopoiesis Alongside the ORANGE Report

Gap G29. Companion to `docs/orange_neoantigen_reporting.md`, which established the pattern this layer
follows: ORANGE runs untouched and its report stays byte-identical, while additional pages are rendered
in the same visual language and appended.

Written against measurements taken from the KQ0004 run on 2026-11-24. Every number below was read from
that run rather than assumed.

## 1. Reporting Principle

**Every variant the caller produces within the panel is reported.** Nothing is removed for being small,
weakly supported or close to the limit of detection. A variant that falls below a stated threshold is
reported with a warning attached, and the evidence behind that warning travels with it: depth at the
position, supporting read count, allele fraction, strand balance and the confidence band it falls in.

The reason is that clone size carries the risk gradient rather than acting as a gate on it. Vlasschaert
and colleagues report an annualised incidence of myeloid malignancy of 0.03% per year for individuals
with small clones, rising to 0.5% per year for those with very large clones
(doi:10.1182/blood.2022018825). A report that silently dropped small clones would be hiding the lower end
of a seventeen-fold risk range. A reader who can see a 3% clone, together with the depth and read count
behind it and a clear statement of what that evidence will support, is better served than a reader handed
a shorter list.

Limitations are reported in the same spirit. The report states the depth achieved in the sample, the
allele fraction the depth can resolve, which genes were interrogated, and which genes were not.

## 2. What the Existing Run Already Provides

| Measurement | KQ0004 | Source |
|---|---|---|
| Blood sample | `KQ0004-2023-AAAD-gDNA`, confirmed peripheral blood | user, 2026-11-24 |
| Mean depth, blood | 127.7x | BamTools coverage histogram |
| Mean depth, tumour | 102.3x | BamTools coverage histogram |
| Bases at 30x or more | 96.5% blood, 96.2% tumour | BamTools coverage histogram |
| Depth across CHIP genes | 100x to 200x bins for DNMT3A, TET2, ASXL1, PPM1D, JAK2, SF3B1, SRSF2, TP53, CHEK2 and CBL | BamTools gene coverage |
| Blood alignment | `KQ0004-2023-AAAD-gDNA.redux.bam`, 186 GB, duplicate-marked with base quality recalibration and microsatellite jitter parameters | REDUX output |

The depth is the reason this is worth doing. At 127.7x a variant needs an allele fraction of 4.8% to
yield three or more supporting reads with 95% power, and 7.0% for five or more reads. The same
calculation at the 30x typical of a tumour and normal pair gives 19.5%, which would leave almost all
clonal haematopoiesis out of reach. At the conventional 2% definition the expected supporting count here
is 2.6 reads and the probability of observing three is 0.47.

## 3. Why The Existing Output Cannot Be Reused

Three layers were examined and none of them carries what a CHIP report needs.

**PURPLE germline.** 833 PASS records, none below 10% allele fraction, the lowest at 0.221. The filters
are tuned for constitutional variants, which sit near 0.5 and 1.0.

**PAVE germline.** 843 records, all PASS, one below 10%. Already downstream of the same filtering, and
therefore no better.

**SAGE germline.** 1,192 records, of which 241 fall below 10% allele fraction with a median of 0.018.
The low-fraction calls do survive here, and they are rejected on `minTumorVAF`, `minTumorSupport` and
`minTumorQual`, which are precisely the criteria a CHIP caller inverts. This file has no gene annotation.

The decisive constraint is scope. `DriverGenePanel.38.tsv` lists 514 genes, and 40 of the 50 canonical
CHIP drivers appear in it, but only 63 genes carry a `reportGermlineVariant` value other than `NONE`.
Every core CHIP gene is set to `NONE`:

| Gene | reportGermlineVariant | likelihoodType |
|---|---|---|
| DNMT3A | NONE | TSG |
| TET2 | NONE | TSG |
| ASXL1 | NONE | TSG |
| PPM1D | NONE | ONCO |
| JAK2 | NONE | ONCO |
| SRSF2 | NONE | ONCO |
| TP53 | ANY | TSG |

SAGE germline consequently produced three records inside DNMT3A and **none at all** inside TET2 or PPM1D.
Those genes are present in the panel for somatic purposes and invisible to germline calling. Since
DNMT3A and TET2 together account for the majority of clonal haematopoiesis, a report assembled from the
existing files would miss most of what it set out to find.

## 4. Design

A dedicated variant-calling step over the blood alignment, followed by annotation, classification and
rendering. No new sequencing is required and the tumour is not involved in the calling.

```
KQ0004-2023-AAAD-gDNA.redux.bam  (186 GB, duplicate-marked, BQR and jitter applied)
   |
   | SAGE, blood treated as the sample of interest, restricted to the CHIP panel,
   | thresholds lowered for allele fraction and supporting reads
   v
chip/<sample>.chip.raw.vcf.gz          every call, nothing discarded
   |
   | PAVE annotation, then annotation against gnomAD, ClinVar and the curated hotspot list
   v
chip/<sample>.chip.annotated.vcf.gz
   |
   | workflow/chip_candidates.py
   v
<case>.chip_variants.tsv               one row per variant, every variant, with its confidence band
<case>.chip_summary.json               depth achieved, detection floor, panel scope, counts per band
<case>.chip.supplement.pdf             appended to the ORANGE report, pages numbered in continuation
```

The tumour alignment is read only to annotate whether a variant is also present there, which separates a
blood clone from contamination of the blood by tumour. For a glioma this is expected to be a formality,
since tumour DNA rarely reaches peripheral blood, and the annotation is carried regardless so that the
question is answered rather than assumed.

## 5. Confidence Bands

Bands describe the evidence. They do not control whether a variant appears.

| Band | Criteria at the depth achieved | Presentation |
|---|---|---|
| Confident | allele fraction 5% or above, five or more supporting reads, strand balanced | reported plainly |
| Supported | allele fraction 2% to 5%, three or four supporting reads | reported with a caution noting it sits below the 95% power floor for this depth |
| Provisional | allele fraction below 2%, or fewer than three supporting reads, or strand imbalance | reported with a prominent caution noting the evidence approaches the error floor of the assay |

The thresholds are derived from the depth of the sample in hand rather than fixed in advance. A sample
sequenced more deeply moves the floor down, and the report states the floor it was computed against.

## 6. Known Artefact Sources To Annotate

Vlasschaert and colleagues record that large-scale CHIP ascertainment unmasks recurrent artefactual
variants and that TET2 and ASXL1 in particular call for gene-specific filtering. Rather than removing
such calls, this layer annotates them:

- recurrent artefact positions in ASXL1 and TET2, flagged by position;
- homopolymer and short tandem repeat context, taken from the `RC_REPS` and `RC_REPC` fields SAGE already
  emits, and from the REDUX microsatellite jitter parameters computed for this sample;
- strand bias, from SAGE's own `readStrandBias` assessment;
- mapping quality, from `MQF` and `mapQualRefAltDiff`;
- population frequency, so that a common inherited variant is distinguishable from a clone;
- allele fraction near 0.5 or 1.0, which indicates a constitutional variant rather than a clone.

## 7. Reporting Scope And Its Limits

The report carries a scope statement naming every gene interrogated and the depth achieved in each. Ten
of the fifty canonical CHIP genes are absent from the Hartwig driver panel and will be added to the CHIP
panel directly: GNB1, BCORL1, ETV6, BRCC3, PRPF8, LUC7L2, DDX41, YLPM1, SRCAP and ZBTB33.

Three limits are stated on the face of the report.

The allele fraction floor means a negative finding says that no clone above the stated size was observed.
It does not establish that clonal haematopoiesis is absent.

The sample is whole blood rather than sorted cells, so a clone's lineage is not determined and the
measured allele fraction reflects the mixture of cell types present.

Clonal haematopoiesis of indeterminate potential is a definition that also requires the absence of a
haematologic malignancy and of unexplained cytopenia. Sequencing alone cannot establish either, so this
layer reports somatic variants in clonal haematopoiesis genes together with their clone sizes, and the
designation itself rests with a clinician holding the blood counts and the clinical picture.

## 8. Consent And Incidental Findings

Clonal haematopoiesis carries cardiovascular and haematologic implications for the person sequenced, and
a finding here is incidental to the reason the tumour was sequenced. This layer therefore sits under the
same consent model the repository already defines for germline findings, and the supplement is generated
only when consent covers it. The research framing of the ORANGE report carries across unchanged.

## 9. Sources

- Vlasschaert C, Mack T, Heimlich JB, et al. A practical approach to curate clonal hematopoiesis of
  indeterminate potential in human genetic datasets. *Blood* 2023;141(18):2214-2223.
  doi:10.1182/blood.2022018825. Open abstract consulted; the full text is not open access and the
  curated gene list and filtering parameters in its supplement have not been read.
- `DriverGenePanel.38.tsv` and `KnownBlacklist.germline.38.bed` from
  `hmf_pipeline_resources.38_v3.0.0--8`.
- Measurements from the KQ0004 oncoanalyser 3.0.0 run as cited in section 2.
