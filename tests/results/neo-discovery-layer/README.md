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

## Outstanding

- The pVACtools container pull (`docker://griffithlab/pvactools:5.4.0`, Slurm job 5f058ee5) was still
  running after 39 minutes. It gates: the VEP plugin install, the predictor inventory
  (NetMHCpan/NetMHCpanEL/MHCflurry availability and licences), and the first pVACseq run.
- `workflow/neo_candidates.py` (unified tiered table joining pVACseq with NEO, LILAC, SAGE_APPEND and
  PURPLE) and the pVACview/IGV review SOP follow once pVACseq has run on KQ0004.
- `VERSIONS.md` to be created with the pinned oncoanalyser, VEP, pVACtools and predictor versions.
