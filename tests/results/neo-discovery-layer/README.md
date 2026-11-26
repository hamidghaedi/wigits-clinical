# Neoantigen discovery layer — development record

Branch `neo/discovery-layer` (from `hartwig_adoptation` at `7a8f207`). Gap G28. Development case:
**KQ0004** (`KQ0004_full_native_20261114_onco3.0.0`), a complete native oncoanalyser 3.0.0 WGTS run
with tumour DNA, matched normal DNA and tumour RNA. `subject_a` (nf-core test case) is the unit-test
fixture only.

## What runs today

| Command | Where | Outcome |
|---|---|---|
| `python3 -m pytest tests/test_neo_interface.py tests/test_neo_export.py -q` | Slurm, 1 CPU | **36 passed** (jobs 79c61ccc, 3ccbf043) |
| `python3 workflow/neo_case_inventory.py out.json` | Slurm, 1 CPU | inventory of the WGTS runs (`case_inventory.md`) |
| `python3 workflow/neo_export.py --results-case-dir <results/case> --outdir <dir> --case <id>` | Slurm, 1 CPU, 8 GB, < 1 min on KQ0004 | four interface files + manifest |
| `vep … --offline --cache --merged --cache_version 115 --use_given_ref --fasta GRCh38_masked_exclusions_alts_hlas.fasta` | Slurm, 4 CPU | annotates the HMF-genome VCF (job 0e3c2ebe) |

## KQ0004 export, 2026-11-23 (job 254c81e6)

Inputs read: 14 source files, each with its sha256 in `KQ0004.interface_manifest.json`.
Samples resolved from the VCF header: tumour DNA `KQ0004-2023-AAAD-2-DNA`, normal `…-gDNA`,
tumour RNA `…-2-RNA`.

| Quantity | Value |
|---|---|
| somatic records in the PURPLE VCF | 16,859 |
| PASS | 12,962 |
| selected (phase-1 coding classes: missense, nonsense/frameshift, in-frame indel) | **72** |
| RNA `SUPPORTED` (DP ≥ 10, alt ≥ 2) | 36 |
| RNA `COVERED_NO_ALT` (DP ≥ 10, alt < 2) | 17 |
| RNA `NOT_COVERED` (DP < 10) | 19 |
| SAGE phase sets among the selected variants | 11 (12 variants) |
| proximal somatic partners | 0 — **verified, not a defect**: every PASS record sharing those 11 phase sets *is* one of the selected variants (12 in-selected, 0 outside) |
| proximal germline variants | 0 — the PURPLE germline VCF holds 833 PASS records (reportable genes only); none within 100 bp of the 72 positions |
| class I alleles (LILAC, de-duplicated, pVACseq spelling) | HLA-A\*02:01, HLA-A\*26:01, HLA-B\*07:02, HLA-B\*38:01, HLA-C\*07:02, HLA-C\*12:03 — LILAC QC `PASS` |
| ISOFOX expression joined | 39,500 genes, 214,977 transcripts |
| purity / ploidy | 0.92 / 2.06, `NORMAL`, MSS, TMB 4.42/Mb (`LOW`), TML 69 (`LOW`) |
| ORANGE immune escape | all six flags `false` (no HLA, antigen-presentation, IFN-γ, PD-L1, CD58 or SETDB1 escape) |

Outputs in `/global/scratch/hpc6049/wigits_runs/neo_ext/KQ0004/`: `KQ0004.somatic.selected.vcf` (115 KB),
`KQ0004.proximal.vcf`, `KQ0004.hla_alleles.txt`, `KQ0004.interface_manifest.json`.

Note on scale: 72 coding candidates from a 4.4/Mb TMB case is the expected order; the native NEO run on
the same case produced 143 candidate events (78 missense + 65 fusion), the difference being that NEO adds
DNA fusion neoepitopes, which phase 1 deliberately leaves with NEO/LINX.

## Issues found and resolved

1. `/cvmfs/.../profile/bash.sh` is not `set -u` safe — source it before `set -u` (also in the launcher).
2. The merged VEP cache is **BAM-edited**: VEP silently enables `--use_transcript_ref`, which replaces the
   VCF's REF allele. `--use_given_ref` is mandatory and is set in `workflow/neo_export.sbatch`.
3. pVACseq requires the `Wildtype` and `Frameshift` VEP plugins. They are **not** in
   `ensembl-vep_release_115.2.sif` (they ship with pVACtools) — the sbatch wrapper points
   `--dir_plugins` at `${KB_ROOT}/vep_plugins/pvactools` and fails the job if no `WildtypeProtein`
   annotation appears.
4. `.gitignore` ignored `*.vcf.gz`, silently excluding the test fixtures; a scoped `!tests/fixtures/**`
   negation now versions them (65 KB, public nf-core test data).
5. `grep … | head` under `set -o pipefail` exits 141 (SIGPIPE) and fails an otherwise successful job.

## First prediction on KQ0004, 2026-11-23

pVACtools **5.4.0** (12 GB image in the shared SIF cache) with MHCflurry 2.0.6 (presentation models bundled
at `/opt/mhcflurry_data/2.0.0/`), MHCnuggets 2.4.1, TensorFlow 2.15.1. Predictors available without a
licence: MHCflurry, MHCflurryEL, MHCnuggetsI/II, BigMHC_EL, BigMHC_IM, DeepImmuno. The NetMHC family is
listed but needs the licensed IEDB standalone, which is **not** installed — noted for `VERSIONS.md`.

Command (job 0d1e2629 for prediction, e37bcd56 for filters, 8e65690f for the aggregate):

```
pvacseq run KQ0004.vep.vcf.gz KQ0004-2023-AAAD-2-DNA \
  HLA-A*02:01,HLA-A*26:01,HLA-B*07:02,HLA-B*38:01,HLA-C*07:02,HLA-C*12:03 \
  MHCflurry MHCflurryEL MHCnuggetsI <outdir> \
  --normal-sample-name KQ0004-2023-AAAD-gDNA -e1 8,9,10,11 \
  -p KQ0004.proximal.vep.vcf.gz --n-threads 8 --pass-only \
  --normal-cov 5 --tdna-cov 10 --trna-cov 5 --normal-vaf 0.02 --tdna-vaf 0.05 --trna-vaf 0.0 --expn-val 0.0
```

| Stage | Rows |
|---|---|
| variant–transcript combinations entering prediction | 1,213 |
| all epitopes (6 alleles x 8–11mers x 3 algorithms) | **275,838** |
| binding filter (500 nM, 2nd percentile) | 1,852 |
| coverage filter (DNA ≥ 10x/VAF ≥ 0.05, normal ≥ 5x/VAF ≤ 0.02, RNA ≥ 5x) | 1,573 |
| transcript support level filter | 130 |
| top score filter (best epitope per variant) | **30** |
| aggregated report (pVACview input) | **69 variants**, + 3.4 MB metrics JSON |

The aggregated report carries per-allele peptide counts, `Best Peptide`, `Best Transcript`, `TSL`, anchor
position, IC50 and percentile for mutant and wild type, `RNA Expr`, `RNA VAF`, `Allele Expr`, `RNA Depth`,
`DNA VAF`, `Tier` and an `Evaluation` column — exactly the pVACview review surface. pVACtools estimated the
tumour clonal VAF as 0.5 from the DNA VAF distribution and treats VAF < 0.25 as subclonal; PURPLE's purity
of 0.92 and `SUBCL` per variant are the better source and are carried in our own table.

Top of the shortlist (gene, AA change): GSX1 S156I, ATRX R2153S, FAAP100 G540R, ARMCX2 V368I,
COL18A1 S294L, MAPK13 P306L, MCCC2 V339M, CDH7 R633W.

## Requirements pVACseq imposes on the inputs (discovered by running it)

1. The **proximal variants VCF must itself be VEP-annotated** — pVACtools checks *its* CSQ header, and the
   error message ("Input VCF does not contain a CSQ header") names the main input, which is misleading.
2. The proximal VCF must contain **exactly one sample**; PURPLE's two-sample header is rejected.
   `neo_export.py` now writes the tumour column only.
3. Both VCFs must be **bgzipped and tabixed**; `bgzip` is on neither the pVACtools nor the VEP image, so the
   wrapper borrows htslib from `samtools-1.22.1--h96c455f_0.sif`.
4. VEP must be run with **`--tsl` and `--biotype`** in addition to `--symbol --terms SO --transcript_version
   --hgvs` and the two plugins.
5. pVACtools 5.4.0's image is **missing its own** `supporting_files/cancer_census_hotspot_gene_list.tsv`, so
   `pvacseq run` dies in post-processing after the predictions are complete. Workaround: run the filter and
   aggregate steps separately (`binding_filter`, `coverage_filter`, `transcript_support_level_filter`,
   `top_score_filter`, `generate_aggregated_report`). A genes-of-interest list built from the COSMIC CGC
   v104 already on the cluster (768 symbols, `data/knowledgebases/genes_of_interest/`) is kept for the
   candidate table; `generate_aggregated_report` in this version takes no such option.
6. The harness runs Slurm scripts under `bash -eo pipefail`: unguarded version probes abort a job, and
   `grep | head` exits 141.

## Outstanding

- `workflow/neo_candidates.py`: join the 69-variant aggregated report with Hartwig NEO's peptide scores
  (`Likelihood`, `RecogSim`, `AlleleCN`, `AlleleDisrupted`), LILAC allele status, SAGE_APPEND RNA support,
  PURPLE `SUBCL`/VCN-derived CCF and the CGC gene list, then assign the Prioritise/Review/Hold/Exclude tiers.
- The NEO-versus-pVACseq concordance table on the shared peptides.
- `docs/neoantigen_review_sop.md` (pVACview + IGV) and the `wigits_releases/<run_id>/neo_review/` layout.
- `VERSIONS.md` pinning oncoanalyser 3.0.0, VEP 115.2 (cache 115 merged), pVACtools 5.4.0, MHCflurry 2.0.6,
  MHCnuggets 2.4.1, and recording that the NetMHC family is unavailable (no IEDB standalone licence).
