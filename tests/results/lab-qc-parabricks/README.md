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

Offline suite126 tests pass,3 existing skips; adapter3 pass. Shell/diff checks
pass. Existing tests cover event timing/duplicate protection and missing-library
scope. A new test refuses a shrunk manifest membership. Full rendering and
real-data execution remain pending until the runtime canary/final reports pass.

Direct review: retain whole laboratory scope and all lane pairs. No caller
change, GPU, raw integrity job, inferred demultiplexing count or source BAM
rewrite. Reuse only immutable completed sources or explicit successful upstream
Slurm dependencies; no timer/file watcher. Current Picard work remains intact.
Per-run grouping removes unnecessary cross-run waiting. Production caps no
longer serialize coverage in groups of4; site scheduler limits remain.
