# workflow

Read-only allocation/use reporting: `report_resources.py --run-dir <run>`.
It writes CSV and standalone HTML summaries; see
[resource reporting](../docs/RESOURCE_REPORTING.md). Future site launches run
it once after Nextflow returns, without changing the native exit status.

Site configuration for Frontenac is in `conf/` (`frontenac.config`, `frontenac_limits.config`, `frontenac_paths.env`); the smoke-test launcher is `run_oncoanalyser_test_frontenac.sh`. All HPC paths come from `conf/frontenac_paths.env`, the machine-readable form of [`docs/directory_layout.md`](../docs/directory_layout.md).

Full-size native runs use [the input-adjusted execution entry](../docs/FULL_NATIVE_EXECUTION.md),
`run_full_native.sbatch`, with a frozen package and the already staged FASTQ
lane samplesheet. It generates initial resources before launching stock Oncoanalyser,
including64-core heavy CPU tasks, and overrides the historical32-core ceiling.
The model is a conservative first-full-case extrapolation pending runtime evidence.

Existing d1.vd1 BAMs use the separate [native BAM entry](../docs/NATIVE_BAM_ENTRY.md):
`run_bam_oncoanalyser_frontenac.sbatch`. It supplies a custom reference through
the supported upstream interface and retains full REDUX and native callers.
The v0.3.0 entry takes [a case YAML](../samplesheets/bam_case.example.yaml),
with `prepare_mates` set explicitly per BAM. It performs optional standard
Picard preparation and sorting/indexing, with no custom initial/final audit
or receipt-based launch validation. See the [current contract](../docs/BAM_INTAKE_CONTRACT.md).

Intended layout once work starts.

```
workflow/
  main.nf              thin wrapper that runs oncoanalyser then the clinical layers
  nextflow.config      profiles for the HPC side
  modules/
    serve/             build a SERVE-compatible actionability directory
    protect/           evidence matching against the ORANGE JSON
    rose/              conclusion sentences from the rules table
    reporter/          PDF and JSON report, including QC-fail variants
  assets/
    rules/             ROSE rules table, germline reporting gene table
    templates/         report templates
```

The default pipeline is unmodified `nf-core/oncoanalyser` 3.0.0 in `wgts` mode.
The chromosome-subset run does not establish a production pipeline defect.
An isolated experimental `3.0.0-qlmp.1` source release can be selected explicitly
to render QC-failed ORANGE reports when a failed PURPLE fit has no circos plot. It retains
the native QC failure and native QC-fail notices; it does not invent a fit or an
image. Missing required analytical inputs, or a missing circos plot without a
supported QC-failure status, raise an explicit error instead of silently
discarding the report.

The maintained change is `patches/oncoanalyser-3.0.0-qc-report.patch` plus
`overlays/oncoanalyser-3.0.0/`. The Java overlay is derived from ORANGE 5.0.1 and
retains GPL-3.0 licensing in its `COPYING` file. All other native classes and
scientific outputs are unchanged. `prepare_oncoanalyser_qc_report.sh` builds a
new, committed/tagged local source release from the pinned upstream checkout;
it never modifies the shared Nextflow asset cache or an installed image.

For the experimental renderer only, set `ONCO_SOURCE` to the local versioned
checkout in `ONCO_SOURCE_DIR`. It is not the default or a validated production fix.

See [QC-failed report recovery](../docs/ORANGE_QC_REPORT_RECOVERY.md) for the
exact build/test/launch instructions and real-data execution evidence.
