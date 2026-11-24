# Laboratory QC component validation

Executed November24,2026 on Frontenac. Runtime:
`lab-qc/v0.1.2`, `ba968961f82cbff8a4126eb96e934fe92270dd09`.
The separate LAB_RUN_QC entry does not run alignment, REDUX or callers.

## Completed

- Offline suite106 tests OK,3 existing integration skips; adapter suite3 tests
  OK;6 focused laboratory tests. Shell syntax and diff checks pass.
- Nextflow26.04.4 metadata preview exits0 with the v1 syntax parser required
  by this named DSL2 entry. Strict-parser rejection was diagnosed before data
  execution; no tool or environment upgrade was made.
- Installed compute-node tool probes completed0 (Slurm12509836,10s).
  fastp reports1.0.1, Picard reports3.1.0-SNAPSHOT, MultiQC reports1.32.
  The installed image named3.1.0 is not used as the executable version string.
- A two-library, two-sequencing-run synthetic canary completed all10 tasks and
  exited0 (Slurm head12509905,13m06s including queue waits). Four raw FASTQ QC,
  two WGS, two insert-size and two run-report tasks completed successfully.
  Both reports have pipeline-written `_SUCCESS`, exactly one main biological
  sample row per run, all required summary fields, matching embedded branding
  image bytes and50px CSS. Actual emitted source tag/SHA match this release.
- The canary reference/BAM dictionaries include an extra nonprimary contig.
  Generated native interval lists retain both dictionary records while measuring
  only the primary contig. Native BedToIntervalList, rather than a manually
  truncated dictionary, establishes this boundary.

## Direct adversarial review

- [BLOCKER, resolved] Flat MultiQC data ZIPs failed the initial finalizer's
  nested-path assumption. v0.1.1 handles flat/nested archives, with regression
  tests. The v0.1.0 canary produced metrics/HTML but was not accepted/published.
- [BLOCKER, resolved] A subset interval dictionary can disagree with the full
  BAM/reference dictionary. v0.1.2 constructs primary BED from the exact FAI
  and uses native BedToIntervalList with the complete dictionary. The new
  nonprimary-contig canary passes.
- [WARNING] Original instrument PF/undetermined/demultiplexing statistics are
  absent from the two inspected real run sources. Assigned FASTQ yield cannot
  replace them. Both reports disclose their absence; no counts are fabricated.
- [NOTE] Fastp is QC-only; its disabled filtering is not library acceptance.
  BAMTOOLS duplicates and raw fastp duplication estimates keep separate meanings.
  Native Picard metrics are not labelled as Parabricks/Hartwig measurements.

## In progress and remaining

Real processing is submitted for two laboratory sequencing runs (4 and10
libraries,112 lane pairs). Slurm head12509906 is running; initial FASTQ tasks
have completed and real Picard coverage is progressing. Conditional delivery
12509908 uses existing coordinator staging/synchronization after successful
completion. No real report completion, scientific clearance or universal
resource/performance equivalence is claimed yet.

The explicit BAM-ready entry is implemented. A future native-queue trigger
remains unwired; this component does not change the default native source.
The pull request could not be created through the GitHub integration (HTTP403).
The feature branch and numbered tags are pushed; main has not been merged.
Detailed real inputs, generated metrics, logs and timing remain on private HPC
storage; the coordinator's run-level delivery receipt records their exact paths.
