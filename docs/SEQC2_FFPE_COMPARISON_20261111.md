# SEQC2 FFPE full-genome comparison

Owner authorized submission on 11 November 2026. Execution package version
0.1.0, development branch `hartwig_adoptation`. Existing KQ0205 jobs are untouched.

## Scope

1. Restore the full SEQC2 six-hour FFPE HCC1395 tumour and HCC1395BL normal
   paired FASTQs. These are the previously selected `FFG_IL_T_6h` and
   `FFG_IL_N_6h` libraries, not six separate patients. **Both libraries are FFPE.**
2. Run native oncoanalyser 3.0.0 on the HMF masked GRCh38 reference, with
   native alignment, REDUX, SAGE and downstream DNA analyses. No RNA input is
   supplied; no synthetic RNA or clinical measurements are added.
3. Following successful control completion, run Google DeepSomatic 1.10.0
   `FFPE_WGS` on the exact same full tumour/normal REDUX BAMs. No regions,
   downsampling, new alignment, or changed caller FILTERs. Batch size is 256.

The DeepSomatic caller arm is queued through a dependent dispatcher. It is
**not yet an end-to-end downstream replacement**: genuine read context and
consumer compatibility remain open in the adapter experiment. No unsupported
SAGE fields are invented and no downstream feature is disabled to claim parity.
Existing historical DNA-pipeline results remain a secondary comparison; this
submission does not rerun every historical pipeline version.

## Input evidence

Resolved afresh with the installed SRA metadata bridge using both run accessions.
Two runs, two experiments, two biological samples, four paired FASTQ objects;
WGS / GENOMIC / PAIRED / Illumina HiSeq 4000. All four expected FASTQs are
available from ENA. There are no submitted-format-only or missing FASTQ objects
among these two explicitly scoped runs. This is not a whole-study inventory.

| Role | Run / experiment / sample | Native isolate | Compressed bytes |
|---|---|---|---:|
| Normal | SRR7890958 / SRX4728410 / SAMN10102582 | ATCC-1395-BL-R-6 | 136,011,567,143 |
| Tumour | SRR7890962 / SRX4728406 / SAMN10102586 | ATCC-1395-BrCA-R-6 | 266,814,507,144 |

Total 402,826,074,287 bytes (402.83 GB; 375.16 GiB). The four live ENA URLs,
byte sizes and MD5s match the prior downloader's explicit FFPE manifest. The
archive metadata confirms six-hour FFPE treatment. The prior explicit library
mapping plus archive isolate identifiers supplies the role assignment; accession
order is not used to infer pairing.

Raw files were removed during the previously authorized scratch cleanup, so
restoration is required. The frozen existing downloader uses two concurrent
HTTPS transfers, resumable `.part` files, published byte sizes and MD5s, then
renames completed objects. It does not delete data. Input evidence and original
archive attributes are in `package_v0.1.0/ena_files.tsv` and `input_mapping.json`.

## Submitted jobs and paths

| Stage | Job | Dependency | At submission |
|---|---:|---|---|
| FASTQ restoration | 12452929 | None | Accepted; pending scheduling |
| Native full workflow | 12452930 | afterok:12452929 | Waiting for restoration |
| DeepSomatic dispatch | 12452931 | afterok:12452930 | Waiting for native completion |

Native root:
`/global/scratch/hpc6049/wigits_runs/SEQC2_FFPE6h_native_20261111_onco3.0.0/`

DeepSomatic root:
`/global/scratch/hpc6049/wigits_runs/SEQC2_FFPE6h_deepsomatic_20261111/`

The native root's `submission.json` contains exact sbatch commands and job IDs;
`package_v0.1.0/release.json` fingerprints the frozen source and input metadata.
The package records its base commit plus file hashes, including the new launch
files; it does not assert that the base commit already contained those files.
`package_v0.1.0/samplesheet.csv` is the two-library upstream input. The tumour's
descriptive cancer type is Breast; these are reference cell-line materials,
not a patient report or a new clinical diagnosis.

`workflow/run_reference_native.sbatch` uses the existing input-v2 resource
generator and native launcher. Heavy CPU tasks receive 64 CPUs and memory/time
scaled from input sizes. The generated config overrides the old site 32-CPU
ceiling. Reference resources and containers are reused. The only runtime code
correction is the already tested GRIDSS 2.13.2 long-count Kraken parser fix,
copied into this package with its hash. No KQ0205 Kraken outputs are reused.

`workflow/submit_reference_deepsomatic.py` runs after the native head succeeds,
checks the two published REDUX BAM/index paths, sizes allocations from actual
BAM bytes, and submits example generation → GPU inference → postprocessing.
Those three job IDs will be recorded in the DeepSomatic root's `submission.json`
when dispatch executes. CPU stages use 64 CPUs; inference uses 16 CPUs plus
one GPU. Whole-genome baseline limits of 60/72/12 hours grow for larger BAMs;
they are conservative limits, not completion predictions. The baseline example
stage measured 19h48; baseline inference has not yet completed. The unfinished
distributed-example experiment is not adopted as the default for this run.

## Evaluation after completion

Use the retained SEQC2 high-confidence SNV/indel truth and territory. Confirm
REF/coordinate compatibility between its d1.vd1 reference and HMF masked genome
before scoring; do not silently omit discordant regions after seeing results.
Score SNVs/indels separately, stratify VAF, inspect DeepSomatic-only truth calls,
and measure downstream allele retention and compute time. Somatic truth does
not validate CNV, SV, MSI or reporting.

HCC1395 was used in DeepSomatic development. The official
[1.10 FFPE case study](https://github.com/google/deepsomatic/blob/r1.10/docs/deepsomatic-case-study-ffpe-wgs.md)
identifies chromosome 1 as held out during training. Treat chr1 as the primary
held-out accuracy comparison; whole-genome results assess integration and
concordance and must not be described as independent validation.

## Verification

`python3 -m unittest discover -s tests`: 52 tests, OK, 3 existing skips.
`bash -n workflow/run_reference_native.sbatch`: exit 0.
`git diff --check`: exit 0.
The additional tests cover workload-scaled allocations, resource envelope and
separate CPU/GPU command construction. Actual analysis completion remains pending.
