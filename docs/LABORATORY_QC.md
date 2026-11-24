# Laboratory sequencing-run QC

Component `lab-qc/v0.1.1`, owner direction November24,2026. Separate
`LAB_RUN_QC` Nextflow entry: existing BAMs and original FASTQs to one branded
MultiQC per **sequencing run**, without running callers or changing native
results. Use after all-lane BAM preparation, not after clinical reporting.
The later native AMBER/PURPLE analytical report remains separate.

## Sources and calculations

- Each manifest sample binds its real sequencing run, exact BAM/index, original
  FASTQ lane pairs, native BAMTOOLS flag counts and emitted software versions.
  A tumour and normal from different sequencing runs appear in different reports.
- fastp1.0.1, QC-only, no read output: disable adapter, poly-G trimming, quality
  and length filtering. Read/base/Q30 counts use **before_filtering**. Sum counts
  across lanes; Q30 is summed Q30 bases divided by summed bases. Individual mate
  reads are not read pairs. This does not add fastp preprocessing to analysis.
- Picard image3.1.0 (the installed executable reports3.1.0-SNAPSHOT) CollectWgsMetrics: exact alignment reference, primary chr1-22/X/Y
  intervals, MAPQ20/baseQ20, cap250, COUNT_UNPAIRED=false, USE_FAST_ALGORITHM=true,
  READ_LENGTH151; no STOP_AFTER. Duplicates and overlapping mate observations
  are excluded by native Picard. Fast accumulation is used without substituting
  Hartwig/Parabricks values; cross-tool identity is not claimed.
- Picard CollectInsertSizeMetrics retains native defaults; General Statistics
  includes the measured FR mean. No deduplication is rerun.
- Native BAMTOOLS flag text is copied byte-for-byte for MultiQC's samtools-format
  parser and labelled with its actual producer. Combined duplicate percentage
  is primary-duplicate reads / primary reads. It is not renamed as Picard's
  PERCENT_DUPLICATION or a separate optical count. Both source and formula are saved.
- Demultiplexing outputs can be supplied per run. Missing PF/undetermined/index
  statistics are explicitly unavailable, never reconstructed from BAM or FASTQ
  assigned yield. Their absence does not block analysis. No instrument loading
  result is inferred. No fastp filter removal is interpreted as a failed library.

Original native files and lineage are retained. One top-level General Statistics
row per laboratory sample; lane-level raw QC remains in expandable children and
an explicit yield table. Software versions come from task output, never image
names. Queen's Pathology logo/CSS use the accepted shared asset supplied by path.
Header displays the immutable component tag; SHA and small-source hashes live
in provenance and `_SUCCESS`. New versioned execution roots preserve old reports.
This implements the applicable presentation/source/release requirements of the
coordinator's MultiQC reporting contract; legacy fastp/Falco preprocessing and
Parabricks source requirements do not describe native REDUX inputs.

## Execute

Use a clean exact tagged checkout and the existing Frontenac site configuration.
No new environment, image or reference download is required. Real manifests and
results remain outside Git under the standard `mohcc_runs/runs/<analysis_id>/`.
The launcher selects Nextflow's supported v1 syntax parser for this DSL2 named
entry; Nextflow26's strict parser does not accept `-entry`.

```bash
export WIGITS_WORKFLOW_DIR=/path/to/exact-release/workflow
sbatch --account=def-hpcg6049_cpu \
  "$WIGITS_WORKFLOW_DIR/run_laboratory_qc.sbatch" \
  /path/to/manifest.json /path/to/execution
```

Manifest keys: `reference`, `fai`, `dict`, `intervals` (Picard interval-list,
not BED), `logo`, `css`, and `samples`. Each sample has `run_id`, `sample_id`,
`bam`, `bai`, `bam_gib`, `bamtools_flags`, `source_versions`, and `fastqs` with
`label`, `r1`, `r2`, `input_gib` per lane pair. Optional `demultiplexing` maps
run IDs to original source-file lists; `run_notes` records factual source gaps.
Original lab identities/lanes govern membership; filenames do not invent it.

CPU/RAM are task-specific: fastp4CPU/6GiB (up to12 concurrent tasks), native
WGS4CPU/16GiB (up to4), insert2CPU/8GiB (up to4). Initial walltime uses actual
compressed FASTQ-pair or BAM size with a conservative2x allowance: fastp
12GiB/hour, WGS20GiB/hour (minimum4h), insert40GiB/hour. These are planning
allowances pending the first full-run measurements, not validated throughput.
No running analytical job or source file is modified.

The workflow publishes per-run `qc/`, HTML/data ZIP, source summaries and
`_SUCCESS`. Deliver the successful run HTML through the coordinator's maintained
`stage_multiqc_reports.py` to `<phase>/qc_reports/<sequencing_run>/` and existing
Laptop_sync. Preserve the compact QC/source package persistently too.

## Validation and limitations

`python -m unittest discover -s tests -p test_laboratory_qc.py` checks weighted
Q30, count invariance, correct duplicate denominator, native table parsing and
run membership. Real tool/Nextflow execution and rendered/data membership checks
remain necessary; source inspection/unit tests alone are not a completed report.
Early automatic triggering from future native queues is separate from this
explicit BAM-ready entry and is not yet deployed. Record full-run evidence in
`tests/results/laboratory-run-qc/` before claiming that completion.

Primary tool references:
[fastp](https://github.com/OpenGene/fastp),
[Picard](https://broadinstitute.github.io/picard/command-line-overview.html#CollectWgsMetrics),
[MultiQC custom content](https://docs.seqera.io/multiqc/custom_content).
