# KQ0004 — neoantigen candidate report (development case)

**Research output. Not a vaccine design, not a manufacturing specification, not a clinical
recommendation.** No candidate here has been through the expert review in
`docs/neoantigen_review_sop.md`; the tiers below are what the pipeline proposes *for* that review.

Generated 2026-11-23T20:25:48+00:00 · rules `docs/neoantigen_extension_proposal.md@2026-11-23` · branch `neo/discovery-layer`
Slurm job `07deef02` on Frontenac · inputs checksummed in `KQ0004.interface_manifest.json`

## 1. Case

| Field | Value | Source |
|---|---|---|
| Tumour DNA / normal DNA / tumour RNA | `KQ0004-2023-AAAD-2-DNA` / `-gDNA` / `-2-RNA` | oncoanalyser 3.0.0, native path |
| Primary tumour | Nervous system; CUPPA **CNS: Glioma (100%)** | ORANGE front page |
| Purity / ploidy | 0.9200 / 2.0600 | PURPLE |
| Microsatellite status | MSS | PURPLE |
| TMB / TML | 4.4243 per Mb (LOW) / 69 (LOW) | PURPLE |
| HLA class I | A\*02:01, A\*26:01, B\*07:02, B\*38:01, C\*07:02, C\*12:03 | LILAC 2.0, QC **PASS** |
| Genetic immune escape | none — all six ORANGE flags false | ORANGE JSON |

The case passes every case-level gate in the review SOP. Note the low mutational burden: at 4.42
mutations per Mb this is a **small** candidate space, and the shortlist should be read as such.

The ORANGE report's driver profile is a textbook glioma: *PTEN* p.Gln261\*, *TP53* p.Phe270Cys, *ATRX*
p.Arg2153Ser and two *MED12* missense variants, with *EGFR*, *KRAS*, *CCND2* and *CHD4* amplified and
*CDKN2A*/*CDKN2B* homozygously deleted, on a chr7-gain / chr10-loss background. Two things in that
profile bear on the neoantigen read:

- **Signature 11 (temozolomide), 3% of the mutation burden.** This tumour has seen an alkylating agent.
  TMZ-induced mutations are genuine somatic neoantigen substrate, but they also arrive late, so several
  will be subclonal — which is what the subclonal flags on the *MED12* and *NPM1* candidates reflect.
- **Signature 1 (5-methylcytosine deamination), 35%.** The dominant process produces exactly the
  C>T/G>A changes the FFPE screen looks for, so the FFPE rule must not be tightened on this case without
  separating the two sources; as run, it flagged nothing.

## 2. What was run

```
PURPLE somatic VCF (16,859)  →  PASS (12,962)  →  coding, PAVE canonical effect (72)
   ↓ SAGE_APPEND RNA support, ISOFOX expression (GX/TX on the tumour genotype)
VEP 115.2 (merged cache, --use_given_ref, Wildtype + Frameshift plugins)
   ↓ 1,213 variant·transcript
pVACseq 5.4.0 · MHCflurry + MHCflurryEL + MHCnuggetsI · 8–11mers · 6 class I alleles
   ↓ 275,838 epitopes  →  binding (500 nM / 2nd pct)  →  coverage  →  TSL  →  top score (30)
aggregated report (69 variants)  →  joined with Hartwig NEO, LILAC, PURPLE, CGC v104
   ↓
69 tiered candidates
```

## 3. Tiers

| Tier | n | Meaning |
|---|---|---|
| **Prioritise** | 6 | strong IC50 **and** strong percentile **and** mutant allele expressed **and** clonal **and** wild-type binds worse |
| **Review** | 28 | strong binding but one of those missing — human judgement decides |
| **Hold** | 29 | evidence missing or binding not strong; 18 because RNA never covered the position, 11 for no strong prediction |
| **Exclude** | 6 | ineligible — all six for alt reads in the matched normal |

RNA support across the 69: **34 SUPPORTED**, 16 covered with no mutant allele,
19 not covered. The third state is the important one — those 19 are *not* evidence of
non-expression, and tiering them as Hold rather than Exclude is deliberate.

![Binding strength versus RNA support]({{artifact:81321163-6279-42c8-a572-dac90e83acc2}})

## 4. The prioritised six

| Gene | Change | Peptide | Allele | IC50 mut (nM) | IC50 wt (nM) | %ile | Gene TPM | RNA VAF | VCN | NEO likelihood |
|---|---|---|---|---|---|---|---|---|---|---|
| GSX1 | S156I | `FTITQLLEL` | HLA-C*12:03 | 24 | 27 | 0.05 | 4.4 | 0.45 | 0.96 | 0.010474 |
| MAPK13 | P306L | `TAAQALTHL` | HLA-C*12:03 | 50 | 2185 | 0.18 | 1.3 | 0.39 | 1.12 | 0.000748 |
| RCOR2 | L294F | `SFKRQVQSM` | HLA-C*07:02 | 84 | 2965 | 0.69 | 44.4 | 0.27 | 1.10 | 0.000032 |
| PPM1F | P227T | `LPTDTEGAL` | HLA-B*07:02 | 132 | 162 | 0.29 | 8.4 | 0.48 | 1.17 | 0.001634 |
| PAPLN | T825M | `QPGASGRSTHM` | HLA-B*07:02 | 205 | 3725 | 0.71 | 1.7 | 0.25 | 1.10 | 0.000024 |
| ZSWIM9 | V91M | `RPAMGPPQPG` | HLA-B*07:02 | 282 | 452 | 0.88 | 3.4 | 0.52 | 0.98 | 0.000056 |

Two observations worth the reviewer's attention:

- **None is a cancer gene, and none of ORANGE's drivers reached Prioritise.** Of the six somatic drivers
  ORANGE reports, four produced a candidate at all — *TP53* F270C and *NPM1* L79F land in Review,
  *MED12* C514S in Hold, *ATRX* R2153S in Exclude (alt reads in the matched normal). *PTEN* p.Gln261\*
  produced **no** candidate: it is a stop-gained variant, and a premature stop yields no mutant peptide
  for class I prediction. That is correct behaviour, not a gap, but it is the kind of absence a reader
  will look for. The CGC v104 genes in the table — *ATRX*, *TP53*, *NPM1*, *MED12*, *PTPN6*,
  *PDGFB*, *BLM*, *ERC1* — all land in Review, Hold or Exclude. Driver status and neoantigen quality are
  independent properties, and the ranking does not reward the former.
- **The differential matters more than the absolute.** *RCOR2* L294F binds at 84 nM while its wild-type
  counterpart binds at 2,965 nM — a 35-fold differential. *PPM1F* P227T binds at 132 nM against a
  wild-type 162 nM, which is barely a differential at all despite comparable absolute affinity. Both are
  Prioritise under the current rules; the reviewer should not treat them as equivalent.

## 5. pVACseq versus Hartwig NEO

Both predictors scored the same case independently. **16,080 peptide·allele pairs were scored by
both**, and all 69 aggregated best-peptides were found in NEO's output.

- Spearman correlation between pVACseq mutant IC50 and NEO presentation likelihood: **−0.49**
  (negative as expected — stronger binding, higher likelihood). Moderate, not redundant.
- NEO likelihoods are small in absolute terms (median 5.2 × 10⁻⁵, maximum 1.5 × 10⁻²): NEO is modelling the
  whole presentation chain, not binding alone, so its scale is not comparable to an IC50 and should not be
  thresholded as though it were.
- The two disagree enough to be worth keeping both. *GALM* V122L is NEO's top-ranked peptide (1.5 × 10⁻²)
  but only 121 nM by pVACseq; *TTC13* is NEO's fourth at 717 nM by pVACseq, which the binding filter would
  discard. Conversely *GSX1* S156I is top-three in both.

Keeping two independent predictors was a design decision, and this is the evidence for it: a single
predictor's ranking would be a different shortlist.

## 6. Tier-versus-pVACseq-tier cross-tab

| | Pass | LowExpr | NoExpr | Poor | Subclonal |
|---|---|---|---|---|---|
| Prioritise | 2 | 3 | 0 | 1 | 0 |
| Review | 8 | 2 | 8 | 8 | 2 |
| Hold | 0 | 2 | 19 | 8 | 0 |
| Exclude | 3 | 0 | 2 | 1 | 0 |

The disagreement in the corners is the point. Three candidates pVACseq calls `Pass` are **excluded**
here for alt reads in the matched normal — a germline/contamination check pVACseq's VAF filters let
through. Three it calls `LowExpr` are **prioritised**, because expression is read from ISOFOX TPM plus
SAGE_APPEND allele-level support rather than pVACseq's internal expression heuristic. And pVACseq's
`Subclonal` calls rest on an estimated clonal VAF of 0.5, which is wrong at purity 0.92; clonality here
comes from PURPLE `SUBCL` and variant copy number instead.

## 7. Caveats

- Class I only. No class II prediction; LILAC types DP/DQ but not DRB1/3/4/5.
- SNVs and indels only. No RNA-derived fusion or splice neoepitopes (pVACfuse, pVACsplice are phase 2).
- The NetMHC family is absent — no IEDB standalone licence — so the most-benchmarked predictor is not in
  the ensemble. See `VERSIONS.md`.
- Zero proximal variants for phasing, verified not a bug: all 12 variants carrying a phase set are the
  only PASS member of their set, and the PURPLE germline VCF holds no PASS record within 100 bp.
- FFPE deamination screening is implemented (C>T/G>A with few unique fragment coordinates) but flagged
  nothing in this case.
- One unit test (`test_expression_reaches_the_tumour_genotype`) still fails on the trimmed `subject_a`
  fixture, which lacks the ISOFOX `GeneName` column the real case has. The production path is verified
  (72/72 records carry gene and transcript expression); the fixture needs regenerating.

## 8. Files

| File | Content |
|---|---|
| `KQ0004.neo_candidates.tsv` | 69 candidates × 65 columns, with `tier` and `tier_reasons` |
| `KQ0004.neo_concordance.tsv` | pVACseq versus NEO on the 69 shared best-peptides |
| `KQ0004.neo_candidates.summary.json` | counts, thresholds, rules version, input paths |
| `KQ0004.interface_manifest.json` | sha256 of all 14 oncoanalyser inputs consumed |
| `pvacseq/MHC_Class_I/*.aggregated.{tsv,metrics.json}` | pVACview session input |

All on the cluster at `/global/scratch/hpc6049/wigits_runs/neo_ext/KQ0004/`.
