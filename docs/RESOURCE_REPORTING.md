# Resource-use reporting

Reporter release: **resource-report/v0.1.0**. Source:
`workflow/report_resources.py`. Python standard library only; no installation,
analysis rerun or additional compute allocation is required.

This implements the coordinator's `docs/SLURM_JOB_EFFICIENCY_GUIDE.md` as a
small read-only report. It observes allocations and accounting; it never
changes jobs, resources, scientific QC or analysis exit status.

## Run and refresh

```bash
python3 workflow/report_resources.py --run-dir /absolute/path/to/run
```

Default outputs under `<run>/results/pipeline_info/resource_usage/`:

- `resource_usage.csv`: one row per distinct native job ID.
- `resource_usage.html`: standalone readable summary, per-stage totals and job table.
- `sacct.psv`: parent and step accounting used in that snapshot.
- `sstat.psv`: provisional live batch measurements where available.

Re-run the same command to refresh these report files. It does not read BAM or
FASTQ contents. Job discovery uses native trace/log submissions, bounded
execution metadata, preparation dispatch records and controller log names.
It includes failed/cancelled attempts and historical work represented by that
run directory. It does not scan every job owned by the user or recurse through
the work/payload tree. Additional standalone jobs can be specified with
`--job-id 123,456`. `--out-dir` selects a different report location.

Saved accounting can be reviewed offline:

```bash
python3 workflow/report_resources.py --run-dir /path/to/run \
  --sacct-file /path/to/sacct.psv --sstat-file /path/to/sstat.psv
```

## What the measurements mean

- Final CPU use: `100 × parent TotalCPU / (AllocCPUS × ElapsedRaw)`.
  Average busy CPUs are `TotalCPU / ElapsedRaw`. Parent totals are not added
  to their step totals. Each native job ID is counted once.
- Running jobs use available `sstat` batch CPU times as explicitly provisional
  estimates. Separate `srun` steps may be absent. The usual running-job
  `sacct TotalCPU=0` is not displayed as zero efficiency.
- Memory is the recorded batch MaxRSS observation where available, or a
  labelled maximum step observation. It is not a guaranteed simultaneous
  whole-job/multi-node peak. Allocated RAM comes from AllocTRES, not ReqMem.
  Nextflow CPU/RSS observations are retained separately in the CSV.
- Allocated core-hours include unsuccessful attempts and controller/transfer
  costs; they are occupancy, not site billing units. Completed-job efficiency
  is weighted by allocated core-seconds, with coverage shown. Support jobs
  and native workflow steps have separate stage rows.
- GPU card-hours describe allocation only; no GPU activity is invented.
  Generic and typed TRES entries for the same GPUs are not added together.
- Blank metrics mean unavailable. CPU figures above100% are flagged for review,
  not capped. Signal-killed job accounting may omit descendant CPU time.

Compare equivalent tasks/input sizes. These reports alone do not establish
that one pipeline is faster or that a memory request can safely be reduced.
Recorded allocation workloads are shown when `state/input_resources.json` exists.

Definitions: [Slurm sacct](https://slurm.schedmd.com/sacct.html) and
[sstat](https://slurm.schedmd.com/sstat.html).

## Minimal future wiring

The shared site launcher calls this script once after Nextflow returns, on
ordinary success or failure. A reporting error emits a message and preserves
the native exit code. It adds no launch validation, monitoring daemon or Slurm
analysis job. Its own controller can still be RUNNING while this final command
executes; a later manual refresh finalizes that row and any delayed accounting.
Preparation-only failures before reaching the shared native launcher can be
reported manually with the same command.

This hook applies only to **future packages built from this source**. The two
current KQ0205 packages remain unchanged; their CSV/HTML snapshots were generated
by invoking the reporter separately. Refresh them manually as needed.

## Development branch

`wigits-clinical:hartwig_adoptation` is the canonical development branch for this
native Oncoanalyser work. It contains the complete earlier BAM, resource and
ORANGE experiment histories. The same branch name in `dna-seq-nf` belongs to
that separate workflow; changes there are not required to run this reporter.
Existing release tags and frozen packages retain their original source identity.
