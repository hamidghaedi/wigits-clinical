# Resource reporter v0.1.0

Development branch: `hartwig_adoptation` in `wigits-clinical`.
The complete native BAM, full-input resource and ORANGE experiment histories
are ancestors of this branch. Existing numbered release tags and active frozen
runtime packages remain unchanged.

## Implementation

`workflow/report_resources.py` produces CSV and standalone HTML from native
trace/log IDs, bounded execution metadata and Slurm accounting. It reads no
genomic input content and performs no scheduler mutation. CSV rows preserve
the difference between final parent accounting and provisional live batch
observations. Per-stage summaries are weighted; duplicate references to a job
are counted once. All represented attempts contribute allocation cost,
including unsuccessful historical work. See `docs/RESOURCE_REPORTING.md`.

The shared launcher invokes reporting once after ordinary Nextflow success or
failure. Development tests confirm that reporter success/failure never changes
the native exit status. This is for future source packages; the two existing
KQ0205 packages were not changed or interrupted.

## Checks

- Slurm12444959:COMPLETED,exit0:0,8s; all22 tests passed with no skips.
- Login suite:22 tests passed,3 real-tool tests skipped as expected.
- Shell syntax and whitespace checks passed.
- Reporter live calls for both run roots exited0 with no collection warnings.
- Coverage tests include days/fractional CPU times and memory units, weighted
  percentages, parent/step separation, live-zero handling, signal-killed jobs,
  GPU allocation deduplication, duplicate trace/log IDs, HTML escaping, offline
  reporting and native-exit preservation by the future completion hook.

## Generated snapshots

Both run roots use `results/pipeline_info/resource_usage/`, containing
`resource_usage.csv`, `resource_usage.html`, `sacct.psv` and `sstat.psv`:

- `/global/scratch/hpc6049/wigits_runs/KQ0205_full_native_20261109_onco3.0.0/`
  initially represented16 jobs.
- `/global/scratch/hpc6049/wigits_runs/KQ0205_full_bam_custom_20261109_r1/`
  initially represented9 jobs.

These are refreshable snapshots, not final whole-run performance claims.
The BAM root includes earlier failed REDUX and cancelled audit attempts; the
FASTQ root includes transfer work and zero-runtime cancelled recovery jobs.
Do not compare aggregate totals as matched BAM-versus-FASTQ performance.
Raw run payload/accounting exports remain on scratch, not in Git.
