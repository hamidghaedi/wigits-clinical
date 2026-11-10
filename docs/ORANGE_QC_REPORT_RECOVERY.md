# ORANGE reports after a failed PURPLE fit

Release: Oncoanalyser `3.0.0-qlmp.1`, QC-report patch `0.1.0`, based on native
Oncoanalyser `7c74c87a43749952b38c9a18915947570f0595a0` and ORANGE `5.0.1`.

Status correction after owner review, November9: this is an isolated renderer
experiment, not an established upstream defect or production fix. The restricted
chr17/chr19 input and absent fit satisfy the native skip condition. Successful
full-input rendering must be evaluated with suitable native inputs; this test
does not justify changing native eligibility. The launcher defaults to unmodified
Oncoanalyser3.0.0. Previously submitted experimental job12441408 is preserved;
its results cannot validate the unmodified reporting path.

## Problem and resulting behavior

The KQ0205 chr17/chr19 baseline completed 39 tasks successfully but emitted no
ORANGE report. PURPLE explicitly logged `failed to find copy number fit` and
`insufficient input Amber and Cobalt data`, wrote `FAIL_NO_TUMOR` / `NO_TUMOR`,
and did not create a circos plot. A two-chromosome execution fixture cannot be
treated as a successful whole-genome copy-number analysis.

The native ORANGE subworkflow required both PURPLE and LINX plots before
scheduling. Its front page also dereferenced a missing circos image, although
the other report chapters already support native QC-failure notices.

The local patch makes report availability reflect that distinction:

- Valid inputs and a circos plot: the existing rendering behavior is retained.
- Native PURPLE `FAIL_NO_TUMOR`, `FAIL_CONTAMINATION` or `FAIL_TINC`, with required
  analytical inputs: run ORANGE and retain the original failure status. When
  the circos image is absent, the front page uses the existing native QC-fail
  notice. No plot, fit, copy number or driver result is fabricated.
- Missing required analytical inputs, or a missing circos image without a
  supported native QC-failure status: report an explicit error. Do not finish
  successfully while silently dropping the requested report.
- LINX plots are optional; the installed ORANGE loader already supports their
  absence. Do not pass the filesystem root as an empty LINX plot argument.

No PURPLE, variant-selection or biological interpretation rule changes. This
is native ORANGE execution-test work, not a change to the finalized SomaReport
template or authorization for clinical issuance.

## Maintained source

- `workflow/patches/oncoanalyser-3.0.0-qc-report.patch`: narrow Nextflow module,
  subworkflow and release-version changes.
- `workflow/overlays/oncoanalyser-3.0.0/lib/OrangeReportEligibility.groovy`:
  tested report eligibility.
- `workflow/overlays/oncoanalyser-3.0.0/assets/orange_qc_report/FrontPageChapter.java`:
  native front page with missing-image/QC-failure handling. GPLv3 source and
  `COPYING` are retained separately from the MIT Nextflow wrapper.
- The ORANGE task compiles this one class against the installed 5.0.1 JAR and
  places it first on the classpath. All other installed classes are unchanged;
  the original container remains read-only. Tool versions record patch0.1.0.

The original native baseline and shared Nextflow asset checkout stay unchanged.
The default launcher selects the unmodified upstream release. This experiment
requires an explicit local source override.

## Build and verify the source release

From the repository root:
```bash
source workflow/conf/frontenac_paths.env
bash workflow/prepare_oncoanalyser_qc_report.sh \
  /global/home/hpc6049/.nextflow/assets/.repos/nf-core/oncoanalyser/clones/7c74c87a43749952b38c9a18915947570f0595a0 \
  "$ONCO_SOURCE_DIR"
```

The destination must be new. The helper applies the maintained patch to the
pinned upstream, copies the source overlays, commits the derived tree and tags
`v3.0.0-qlmp.1`. Record its exact resulting SHA with the execution. A released
destination is not overwritten; a further fix requires a new numbered version.

Eligibility regression (13 cases, no genomic work):
```bash
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
module load StdEnv/2023 nextflow/26.04.4
nextflow run tests/test_orange_report_eligibility.nf \
  -lib workflow/overlays/oncoanalyser-3.0.0/lib -ansi-log false
```

## Exact real-data recovery

`samplesheets/KQ0205_chr17_19_orange_qc_report_r1.csv` retains the three original
FASTQ/sample bindings and adds native `filetype` entries for the completed
SAGE, PURPLE, LINX, QSEE, RNA and ancillary directories. They all come from the
same native baseline. No prior full-genome fit or other case's data are mixed in.

Run the owning pipeline in its supported report-only recovery mode. Alignment,
calling and other completed analytical steps are not repeated:

```bash
sbatch --account=def-hpcg6049_cpu --time=02:00:00 --cpus-per-task=2 --mem=8G \
  --export=ALL,ONCO_SOURCE=/global/project/hpcg6049/somatic_pipeline/software/oncoanalyser/3.0.0-qlmp.1 \
  --output=/global/scratch/hpc6049/wigits_runs/KQ0205_chr17_19_orange_qc_report_r1.head.%j.log \
  workflow/run_oncoanalyser_frontenac.sh \
  KQ0205_chr17_19_orange_qc_report_r1 \
  samplesheets/KQ0205_chr17_19_orange_qc_report_r1.csv \
  --processes_manual orange -c "$PWD/workflow/conf/orange_recovery.config"
```

The head uses two CPUs/8GiB, and ORANGE uses two CPUs/8GiB/one hour. This is a
bounded report render over small existing tables/VCFs, not a new genomic
analysis or a claim of a universal whole-genome resource model.

Submit from the repository root. Slurm spools the head script, so the launcher
resolves configuration through `SLURM_SUBMIT_DIR/workflow` when its own directory
does not contain the configuration. `WIGITS_WORKFLOW_DIR` is an explicit override
for submission from another directory. The initial head job12441406 failed
before launching Nextflow because of the former script-relative lookup; the
corrected launcher is tested for direct, spooled, explicit and missing paths
with `python3 tests/test_frontenac_launcher.py` (four paths pass).

Acceptance: successful ORANGE task, nonempty/parseable ORANGE JSON and PDF,
preserved `FAIL_NO_TUMOR`, null/unavailable circos metadata rather than a fake
image, and visible native QC-fail text. The trace must not contain new alignment
or calling jobs. Record job IDs, exit codes, artifact paths and checks below.

Execution is pending until those facts are recorded. The separate month-only
clinical sampling-date and native DOID limitations remain as documented in the
metadata review; this recovery does not invent missing clinical data.

## Sources

- Oncoanalyser `subworkflows/local/orange_reporting/main.nf` and
  `modules/local/orange/main.nf`, pinned revision above.
- [ORANGE front page](https://github.com/hartwigmedical/hmftools/blob/orange-v5.0.1/orange/src/main/java/com/hartwig/hmftools/orange/report/chapters/FrontPageChapter.java).
- [Native QC-fail handling](https://github.com/hartwigmedical/hmftools/blob/orange-v5.0.1/orange/src/main/java/com/hartwig/hmftools/orange/report/chapters/PurplePlotsChapter.java).
- [Optional LINX plot loading](https://github.com/hartwigmedical/hmftools/blob/orange-v5.0.1/orange/src/main/java/com/hartwig/hmftools/orange/algo/linx/LinxDataLoader.java).
