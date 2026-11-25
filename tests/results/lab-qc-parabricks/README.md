# Accelerated laboratory run QC — November25

Owner directs Parabricks coverage, laboratory-manifest membership, completed-BAM
event triggering, Phase1 presentation and two complete October-run reports.
The v0.4.0 candidate keeps native analysis/callers unchanged. It replaces only
the separate lab coverage collector, propagates exact run membership and source
metadata, permits explicit small-metric reuse and groups each run independently.

Reference: coordinator `completed_runs/260220_A01939_0077_AHVG7NDSXF/multiqc_out/260220_A01939_0077_AHVG7NDSXF/`
HTML/YAML, shared network logo, and base YAML. Carry visual settings and real
metadata structure, not old hardcoded software versions. Actual Parabricks,
Picard insert-size, BAMTOOLS and fastp sources remain labelled separately.

Offline suite129 tests pass,3 existing skips; adapter3 pass. Shell/diff checks
pass. Existing tests cover event timing/duplicate protection and missing-library
scope. A new test refuses a shrunk manifest membership. Full rendering and
real-data execution remain pending until the runtime canary/final reports pass.

Direct review: retain whole laboratory scope and all lane pairs. No caller
change, GPU, raw integrity job, inferred demultiplexing count or source BAM
rewrite. Reuse only immutable completed sources or explicit successful upstream
Slurm dependencies; no timer/file watcher. Current Picard work remains intact.
Per-run grouping removes unnecessary cross-run waiting. Production caps no
longer serialize coverage in groups of4; site scheduler limits remain.


## Depth-aware memory recovery

The first v0.4.0 full batch failed when a141x tumour exceeded64GiB after2m51.
Nextflow terminated the other coverage tasks. v0.4.1 reads the already-completed
BAMTOOLS summary for an initial depth-aware allowance (128GiB for this tumour),
retains32CPU/filters and changes failure handling to `finish`, with no retries.
The depth/size heuristic is provisional; measured failure establishes only that
64GiB was insufficient, not the exact cause or a guaranteed upper bound.
The metadata canary tests band boundaries, missing/malformed summary, large BAM
allowance and walltime floor; full-size recovery remains separately recorded.

Metadata canary12511665 completed exit0 in7s under Slurm/Nextflow26.04.4.
The129-test offline suite passes (3 pre-existing skips); the adapter3 tests,
shell syntax and diff checks pass. This proves directive calculation, not
full-size coverage completion. Runtime details remain in the private coordinator
`phase2-laboratory-parabricks-20261125-r1` execution record.

## Full-size delivery verification

The v0.4.1 full batch12511667 completed exit0 in1h32m44. Twelve new coverage
tasks, the dependent reuse of the successful matched tumour and both run-level
render tasks succeeded. All29 completed small-source reuse tasks were cached.
The two reports contain exactly4/10 laboratory samples, all eight required
summary metrics populated, and no duplicate REDUX display rows. The completed
task `_SUCCESS` hashes match the persistent canonical HTMLs; the corresponding
Laptop_sync Drive copies match in size and MD5. Source releases, metrics,
original BAMs and scientific/identity holds are preserved. The coordinator's
`docs/qc/phase2_laboratory_parabricks_20261125/README.md` records the two direct
links and retained exact execution receipts. This is full-size validation for
these14 libraries, not a universal resource/performance bound.
