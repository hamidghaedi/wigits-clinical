# Full-run runtime recovery

Recovery component **0.1.0**, prepared 10 November 2026. Owner request: repair
both failed workflows using completed work, without repeating alignment or
interrupting the separate full-genome DeepSomatic jobs.

The original nf-core/oncoanalyser3.0.0 source remains
`7c74c87a43749952b38c9a18915947570f0595a0`. No shared container, reference,
samplesheet, caller settings, resource policy or executing task is modified.

## Repairs

| Failure | Evidence | Correction and limits |
|---|---|---|
| Native VirusBreakend | `NumberFormatException` for Kraken count3,576,733,308 in `KrakenReportLine` | GRIDSS2.13.2 declares both read-count fields as long but parses them using Integer.parseInt. Replace only those two calls with Long.parseLong in an isolated copy of its JAR. Taxonomy IDs remain int; counts and scientific logic are unchanged. |
| Custom-reference CIDER | `FileNotFoundException: GRCh38.d1.vd1.fa.dict` while the correct `GRCh38.d1.vd1.dict` is staged | Add a relative filename alias in this task directory before staging. No reference regeneration or content change. |
| Custom-reference ESVEE | JVM22 SIGSEGV in C2-compiled `htsjdk.samtools.AbstractBAMFileIndex.query` | Exclude this method from JIT compilation using JAVA_TOOL_OPTIONS, passed through an Apptainer env-file. Keep the same Java, tool, BAMs and settings. This is a targeted crash workaround, not proof of the underlying JVM defect; successful full execution still needs observation. |
| Custom-reference AMBER | Slurm CANCELLED during workflow shutdown | Resume normally; there is no demonstrated independent analytical failure. |

Parser source:
[GRIDSS v2.13.2 KrakenReportLine.java](https://github.com/PapenfussLab/gridss/blob/v2.13.2/src/main/java/au/edu/wehi/idsv/kraken/KrakenReportLine.java).
Native source files inspected: `modules/local/{virusbreakend,cider,esvee}/main.nf`.
The patch is `workflow/patches/gridss-2.13.2-kraken-long-counts.patch`;
build/test entry is `workflow/build_gridss_count_fix.sbatch`.

## Resume mechanism

`workflow/prepare_run_recovery.py` freezes the original recorded Nextflow command,
adds its explicit original session ID and a final process-specific recovery
configuration. `workflow/resume_completed_run.sbatch` runs that command from the
original run directory. It bypasses BAM preparation and resource regeneration.
The original work directory and Nextflow cache remain the resume authority.

The native resume waits for the old controller to finish, allowing its active
ESVEE job to finish and enter the cache. It cannot race the active controller.
The completed Kraken read-name list (~95GB), report and viral-contig list are
linked read-only by use into the new task's working directory. VirusBreakend's
existing completed-file checks skip that scan; outputs and logs are new. No
counts are capped, rows removed, outputs fabricated or analyses disabled.

External package:
`/global/scratch/hpc6049/wigits_runs/recovery_20261110_v0.1.0/`.
Use `bam_r2/` and `native_r2/`; earlier preparations are retained as evidence.
Each package records its original command, session, sample-sheet hash and files.

## Verification and execution

- Unit suite:34 tests, exit0,3 existing integration skips. Bash syntax and
  `git diff --check` passed.
- Parser regression checks0,50,2^31−1,2^31,3,576,733,308 and5,000,000,000.
  The original JAR reproduces NumberFormatException; the patched JAR preserves
  both counts and original text and passes. Build job12449088 completed0:0.
- Earlier build12449086 passed the count test but treated native Picard help's
  exit1 as failure. Corrected the harness to accept help's documented behaviour
  only when expected argument help is present; no analytical change.
- JVM smoke12449092 completed0:0 and printed the exact method exclusion.
  Earlier smoke12449089 exposed Apptainer's comma parsing in `--env`; env-file
  fixes this. Its dependent head12449090 was cancelled before workflow execution.
- Custom BAM resume head: **12449093**, dependent on successful JVM smoke.
- Native FASTQ resume head: **12449091**, dependent on terminal old head12444653
  and successful parser build12449088. Its currently running ESVEE12446323 is
  preserved.
- DeepSomatic12447879/12447880/12447881 are independent and untouched.

Submission is not completion. Confirm CACHED upstream tasks and actual recovery
outcomes from the new trace/log before claiming either workflow is complete.

Observed startup: BAM resume12449093 reports12 CACHED tasks, including both
REDUX outputs, all three reference preparation tasks, both SAGE callers, BAM
metrics, COBALT, TEAL and VirusBreakend. Only CIDER12449094, AMBER12449095 and
ESVEE12449096 have been submitted at this point (pending scheduler resources).
No alignment or BAM preparation was submitted. Native resume12449091 is still
waiting for the old native head; DeepSomatic preparation12447879 remains running.
