# Case inventory — neoantigen-relevant oncoanalyser outputs

Generated 2026-11-23 by `workflow/neo_case_inventory.py` (read-only; Slurm job, 1 CPU). Source runs are native
oncoanalyser 3.0.0 WGTS runs under `/global/scratch/hpc6049/wigits_runs/` (KQ cases) and the nf-core test
run (`subject_a`). Full JSON: `case_inventory.json` beside this file.

| case | purity | ploidy | PURPLE qc | somatic recs | NEO cands | NEO peptide rows | NEO types | LILAC alleles | LILAC qc | RNA recs (DP>=10 / alt>=2) | ISOFOX genes | pass fusions | LINX neoep | ORANGE |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| KQ0004 | 0.9200 | 2.0600 | None | 16859 | 143 | 72168 | {'INFRAME_FUSION': 3, 'OUT_OF_FRAME_FUSION': 62, 'MISSENSE': 78} | DPB1*02:01,DPB1*16:01,DPA1*01:03,DPA1*01:03,A*02:01,A*26:01 | PASS | 16859 (472 / 900) | 39500 | 2922 | 36 | True |
| KQ0005 | 0.4600 | 2.1000 | None | 10038 | 167 | 173988 | {'OUT_OF_FRAME_FUSION': 109, 'INFRAME_FUSION': 17, 'MISSENSE': 39, 'FRAMESHIFT': 2} | DPB1*04:01,DPB1*11:01,DPA1*01:03,DPA1*02:01,A*02:01,A*03:01 | WARN_UNMATCHED_INDEL | 10038 (155 / 246) | 39500 | 1871 | 52 | True |
| KQ0012 | 0.8200 | 2.0200 | None | 12062 | 59 | 19680 | {'OUT_OF_FRAME_FUSION': 2, 'MISSENSE': 56, 'FRAMESHIFT': 1} | DPB1*02:01,DPB1*04:01,DPA1*01:03,DPA1*01:03,A*02:01,A*02:05 | PASS | 12062 (189 / 362) | 39500 | 3501 | 2 | True |
| KQ0013 | 0.8900 | 1.9800 | None | 11696 | 65 | 19890 | {'OUT_OF_FRAME_FUSION': 9, 'INFRAME_FUSION': 3, 'MISSENSE': 51, 'FRAMESHIFT': 2} | DPB1*01:01,DPB1*04:01,DPA1*01:03,DPA1*02:01,A*01:01,A*03:01 | PASS | 11696 (112 / 216) | 39500 | 1286 | 12 | True |
| KQ0017 | 0.7100 | 5.0000 | None | 9576 | 81 | None | {'OUT_OF_FRAME_FUSION': 28, 'INFRAME_FUSION': 7, 'MISSENSE': 42, 'FRAMESHIFT': 3, 'INFRAME_DELETION': 1} |  | None | None (None / None) | None | None | 17 | False |
| KQ0205 | 0.2400 | 2.9000 | None | 7347 | 23 | 9426 | {'OUT_OF_FRAME_FUSION': 7, 'MISSENSE': 14, 'INFRAME_DELETION': 1, 'FRAMESHIFT': 1} | DPB1*01:01,DPB1*04:01,DPA1*01:03,DPA1*02:01,A*01:01,A*02:01 | PASS | 7347 (1227 / 413) | 39500 | 3144 | 2 | True |
| subject_a | 1.0000 | 2.0000 | None | 103 | 6 | 5280 | {'OUT_OF_FRAME_FUSION': 2, 'INFRAME_FUSION': 1, 'MISSENSE': 3} | A*03:01,A*03:01,B*07:02,B*07:02,C*07:02,C*07:02 | FAIL | 103 (1 / 0) | 39500 | 2 | 2 | True |

"NEO cands" = rows of `neo/finder/*.neo_data.tsv` (one per neoepitope-generating event); "NEO peptide rows" =
rows of `neo/scorer/*.neo.peptide_scores.tsv` (peptide × allele). "RNA recs" = records in the tumour
`sage_append` VCF, with the RNA sample column at DP ≥ 10 / alt reads ≥ 2. "PURPLE qc" reads `qcStatus` from the
purity TSV, which does not carry that column in 4.4 (it is in `*.purple.qc`); fix in the next revision.

## Findings

1. **Four complete WGTS cases usable now:** KQ0004, KQ0005, KQ0012, KQ0013 — PURPLE, NEO (finder + scorer),
   LILAC (PASS, except KQ0005 `WARN_UNMATCHED_INDEL`), SAGE_APPEND with the RNA sample appended, ISOFOX, LINX
   neoepitopes and ORANGE JSON all present. Purity 0.46–0.92. Plus KQ0205 (native, purity 0.24, ploidy 2.9,
   best RNA coverage of the set) and subject_a (LILAC `FAIL` — tiny test data, class I only).
2. **KQ0017 was still running** when the inventory was taken (Slurm job `KQ0017-native-WGTS`), which is why only `neo/finder` output existed. Not a failure.
   Excluded until the run is finished or re-run.
3. **LILAC 2.0 types class II DP and DQ natively** (`DPA1/DPB1/DQA1/DQB1`, with tumour copy number and RNA
   fragment counts per allele) in addition to class I `A/B/C`. It does not report `DRB1/3/4/5`. The class II
   gap in the proposal therefore narrows to DR typing plus class II binding prediction; NEO itself scores class I
   alleles only (`A0201`-style names, no `HLA-` prefix, no separator).
4. **NEO candidate classes observed:** MISSENSE, FRAMESHIFT, INFRAME_DELETION, INFRAME_FUSION, OUT_OF_FRAME_FUSION.
   Out-of-frame fusions dominate the counts in KQ0004/KQ0005 (62 and 109); these are DNA-derived (LINX) and
   many will lack RNA support — a ranking question for the candidate table. `AlleleDisrupted` is `false` for every
   peptide row in every case (no HLA LOH/somatic hit called in this set).
5. **RNA support is sparse at the whole-genome level by construction** (most SAGE somatic records are non-coding):
   112–472 records with RNA DP ≥ 10 per KQ case, 1,227 in KQ0205. For coding candidates the fraction will be far
   higher; the adapter must distinguish "not covered" from "covered, mutant allele absent".
6. ISOFOX gene table has 39,500 genes in every case (same transcript annotation); `pass_fusions` 1,286–3,501 per
   case — a large RNA-only fusion space for phase 2 (pVACfuse).
7. SAGE_APPEND sample order is `gDNA, tumour DNA, tumour RNA`; the RNA column is last. Sample IDs follow
   `KQ0004-2023-AAAD-{gDNA,2-DNA,2-RNA}`.

## Consequences for the interface specification (Phase 1)

- Class I alleles for pVACseq come from LILAC (`A*02:01` → `HLA-A*02:01`); DP/DQ alleles are available for class II
  prediction once a predictor is provisioned; DRB1 needs an external typer.
- NEO's allele naming (`A0201`) must be normalised to join with LILAC and pVACseq.
- Variant classes to carry through VEP/pVACseq in Phase 1: missense, frameshift, in-frame indel (DNA fusions stay
  with NEO/LINX; RNA fusions are phase 2).
- Use the SAGE_APPEND RNA column for `RAD/RDP/RAF`; emit an explicit `rna_coverage_status` field.
