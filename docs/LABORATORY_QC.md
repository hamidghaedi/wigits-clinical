# Laboratory sequencing-run QC

Component `lab-qc/v0.4.1`, owner direction November24,2026. Separate
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
- Owner-selected Parabricks4.6.0-1 bammetrics for coverage,32 CPU threads,
  MAPQ20/baseQ20/cap250, unpaired=false, chr1-22/X/Y whole-chromosome BED from
  the exact alignment reference FAI. No GPU, sampling, STOP_AFTER or BAM rewrite.
  Native Picard coverage is historical/comparison evidence, never relabelled as
  the same execution. The matched normal comparison took8m39 versus70m10 and
  differs slightly numerically; preserve the actual source values and labels.
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
names. Phase1 MOHCCN network logo at its original natural size and default MultiQC presentation are
locked in `workflow/assets/lab_multiqc_style.json` and `lab_multiqc.css`.
The source comparison is the coordinator's February20 run19 YAML/HTML.
Headers retain Project/Pipeline/Batch/Run/Tissue/Flowcell/Run Date/Report Date;
values come from the laboratory manifest and report preparation, never old literals.
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

Manifest keys: `reference`, `fai`, `dict`, `logo`, `css`, and `samples`.
The native task generates its measured interval list. Each sample has `run_id`, `sample_id`,
`bam`, `bai`, `bam_gib`, optional `bamtools_summary`, `bamtools_flags`, `source_versions`, and `fastqs` with
`label`, `r1`, `r2`, `input_gib` per lane pair. Optional `demultiplexing` maps
run IDs to original source-file lists; `run_notes` records factual source gaps.
Original lab identities/lanes govern membership; filenames do not invent it.

CPU/RAM: fastp4CPU/6GiB (up to12), insert2CPU/8GiB (up to4),
coverage32CPU with `lab-coverage-depth-v1` initial RAM:64GiB at measured native
BAMTOOLS mean coverage<=50x, plus32GiB per additional50x band. Missing/unreadable
summary uses128GiB, with no new scan or launch gate. Add32GiB per BAM-size band
above512GiB. Thus the141x tumour requests128GiB; a23x normal64GiB. The 64GiB
BAM-size-only floor failed in the full batch at141x (peak64GiB, Slurm OOM).
This depth allowance is a conservative recovery heuristic, not universal validation
or a demonstrated causal model of the tool's memory use. See
`workflow/lib/LabQcResources.groovy` and the metadata-only boundary canary.
Coverage walltime stays1min/GiB BAM with1h floor. There is no artificial library
concurrency cap; the site's existing queue and scheduler control concurrency.
A failed coverage task uses Nextflow `finish`: already-running tasks finish,
the workflow still fails, and no automatic resource escalation/retry is added.
Completed tasks are reusable on explicit resume. Existing issued jobs keep their
original resources. No running analytical job or source file is modified.

The workflow publishes per-run `qc/`, HTML/data ZIP, source summaries and
`_SUCCESS`. Deliver the successful run HTML through the coordinator's maintained
`stage_multiqc_reports.py` to `<phase>/qc_reports/<sequencing_run>/` and existing
Laptop_sync. Preserve the compact QC/source package persistently too.

## Validation and limitations

`python -m unittest discover -s tests -p test_laboratory_qc.py` checks weighted
Q30, count invariance, correct duplicate denominator, native table parsing and
run membership. Real tool/Nextflow execution and rendered/data membership checks
remain necessary; source inspection/unit tests alone are not a completed report.
The current October15/19 full-size QC jobs retain their frozen v0.1.2 sources.
Their completion remains separate from validation of automatic dispatch.
See `tests/results/lab-qc-dispatch/` for dispatch tests and metadata evidence.

Primary tool references:
[fastp](https://github.com/OpenGene/fastp),
[Picard](https://broadinstitute.github.io/picard/command-line-overview.html#CollectWgsMetrics),
[MultiQC custom content](https://docs.seqera.io/multiqc/custom_content).


## Automatic BAM-ready dispatch

Owner correction November24: **no recurring readiness check**. Component
lab-qc/v0.4.1 retains the event route and supersedes v0.2.0's five-minute dispatch. The archive sync still
copies completed reports on its existing schedule; it no longer starts QC.

A small operational subscription to native `bamtools_metrics_dir` receives each
successfully completed BAMTOOLS output. BAMTOOLS already consumes the completed
sample REDUX BAM. `LabQcEvents.bamReady` calls `dispatch_laboratory_qc.py`
v2.2.1 once for that event. The first events record readiness; the final expected
DNA library starts one LAB_RUN_QC job for its sequencing run. No SAGE, PURPLE,
ORANGE or analytical-head completion dependency. Cached outputs re-emit the
same topic on resume; persistent state prevents duplicate submissions.

Use the actual completed task BAM/index/flags and REDUX/BAMTOOLS versions.
Do not wait for asynchronous publish copies or final combined versions.
[Nextflow's publishDir documentation](https://docs.seqera.io/nextflow/reference/process/directives/publish-dir)
explicitly directs downstream consumers to task outputs. Metadata only: no
BAM/FASTQ content reads, counting, new raw-data hashing or scientific thresholds.
The laboratory CSV retains the whole DNA run membership and every recorded lane.
Missing libraries wait visibly; no partial run is labelled complete.

`prepare_lab_qc_event_source.py` creates/reuses an isolated clean source from
an exact base commit, adding only a topic subscription and `LabQcEvents.groovy`.
Native analytical modules/configuration are unchanged. Base and derived commits
and hook digest are recorded. The coordinator's new-package launcher v1.2.0
binds this source and a frozen QC config in the case package, then exports
`WIGITS_LAB_QC_CONFIG`. Original installed sources and issued packages remain
untouched. No plugin, new environment or daemon is required.

Each completed run emits `_SUCCESS` to the report-ready callback immediately;
the QC head also calls `--deliver-run` after its own Nextflow execution succeeds. This stages HTML and compact evidence to the canonical phase/run QC
folder; existing Laptop_sync copying continues. Failed QC does not trigger
successful delivery; ambiguous submissions need operator recovery. File locking
serializes concurrent events without dropping the final event. Per-run saved
intent/state prevents duplicate submissions. Errors in the operational callback
are visible warnings and do not alter native analysis success or scientific holds.

Private configuration retains native tracker queue pointers, laboratory CSVs,
explicit reference profiles, frozen workflow/dispatcher/Python, account,
archive/run/state roots, logo/CSS and no-backfill baseline. It additionally gives
`native_source_root`, the shared installation directory for isolated event sources.
Operational state is `<state_directory>/bam_events/<analysis_id>/<sample>.json`
and `<sequencing_run>.run.json`. `last_tick.json` is the last **event** evaluation,
not evidence of a timer. Existing externally managed QC remains unchanged.

The event CLI accepts only explicit events:

```bash
python /path/to/frozen-release/workflow/dispatch_laboratory_qc.py \
  --config /path/to/frozen-case-config.json --bam-ready /path/to/completed/bamtools_SAMPLE
python /path/to/frozen-release/workflow/dispatch_laboratory_qc.py \
  --config /path/to/frozen-case-config.json --deliver-run SEQUENCING_RUN
```

Scope remains registered native DNA queues. RNA-only queues, ambiguous multi-run
BAMs and repeated chunk labels require an explicit released mapping. A changed
laboratory scope after delivery requires an explicit new report revision, not a
silent overwrite. Full-size v0.1.2 metric jobs continue; event validation is not
universal QC resource/performance validation.


## Reuse and independent report completion

`expected_samples` maps every selected sequencing run to its complete sample
list derived from the laboratory manifest. Both the workflow and report helper
retain exact membership. Each run groups its own expected metric events using
Nextflow groupKey; another run cannot delay its completed report.

Optional per-sample `reuse` has `fastp`, `insert`, `coverage` entries, each a
`files` map from normalized metric filenames to immutable source paths. Raw QC
and insert metrics can be copied by small local tasks inside the Slurm head;
large BAM/FASTQ sources are never copied. Coverage reuse must explicitly name
producer `parabricks.bammetrics`. An optional numeric `dependency` is a Slurm
afterok binding for an already-running producer; a small copy task starts only
after it succeeds. There is no file-readiness polling and no duplicate coverage
calculation. Full source paths and small-file hashes remain in report provenance.

For the October15/19 batch all14 BAMs already exist. Reuse112 completed raw-QC
lane products and14 insert metrics. Reuse the matched Parabricks normal output
and its running tumour producer, then calculate the remaining12 coverage sets.
Original Picard metrics and jobs remain evidence. Superseded queued publication
must not overwrite the owner-selected Parabricks reports. The final two HTMLs
use the existing phase/run delivery paths and exact source/sample receipts.
