# Task-input sizing for native RNA

Policy **input-v5**, November24,2026. Development is from main v0.1.0 on
`codex/rna-input-resources`. New execution packages pin the tested commit;
the repository-wide main release is unchanged. DNA input-v4 resource rules and
analytical arguments are unchanged. Do not overwrite an executing package.

## Evidence and failure

Five completed full-size native RNA cases and one failed case supply26 task
attempts:22 successful and four failed duplicate-marking attempts.
The measurements below are per completed task, excluding queue time. Inputs
are literal stage-file sizes from existing metadata/cleanup receipts, not new
read scans. MaxRSS is the Slurm batch-step observation, not a guaranteed
simultaneous whole-job maximum. All full cases used64 CPUs for STAR, sort and
ISOFOX, and8 CPUs for Picard.

| Completed step | Jobs | Actual input range | Elapsed range | Peak RSS range |
|---|---:|---|---|---|
| STAR | 6 |55.55–78.57 GiB, all compressed RNA FASTQs |1.14–3.26h |64.85–96.25 GiB |
| samtools sort | 6 |1102.44–2373.29 GiB, unsorted BAM (five known sizes) |0.89–3.19h |53.05–54.01 GiB |
| MarkDuplicates | 5 |111.02–239.16 GiB, sorted BAM (four known sizes) |3.17–8.17h |80.29–93.73 GiB |
| ISOFOX | 5 |138.67–298.47 GiB, marked BAM (four known sizes) |0.32–0.67h |59.27–123.28 GiB |

KQ0017 MarkDuplicates twice ended at11h59m with exit140, using about93.5GiB
of112GiB allocated. The retained command lacks `--READ_NAME_REGEX null`; its
second log remains in optical comparisons. The package predates the validated
[RNA optical-comparison correction](RNA_DUPLICATE_MARKING.md). An earlier
attempt reached output writing before its limit, so do not describe all elapsed
time as optical work. These observations support insufficient walltime and the
missing RNA option; they do not indicate memory exhaustion.

Compressed FASTQ size alone is a poor proxy for downstream RNA: STAR emits
uncompressed, multimapping-containing BAMs. Use the actual input of each task.

## Initial request rules

Let G denote that task's actual input GiB (1024 cubed bytes). Time is rounded
up to whole hours, except duplicate marking uses two-hour increments. Memory
is rounded up to16GiB. The2x multiplier and startup allowance are explicit safety
margins, not a fitted prediction or universal validation.

| Step | Input | CPUs | Memory GiB before rounding | Walltime h before rounding |
|---|---|---:|---|---|
| STAR | Sum of all forward/reverse FASTQs passed to this task |64 |max(96,96+0.5G) |max(2,2×0.042G+0.5) |
| samtools sort | Unsorted BAM |64 |80 |max(2,2×0.0016G+0.5) |
| MarkDuplicates | Sum of input sorted BAMs |8 |max(128,96+0.125G) |max(12,2×0.060G+1) |
| ISOFOX | Duplicate-marked BAM |64 |max(160,96+0.25G) |max(2,2×0.003G+0.5) |

Time rates round above the largest observed elapsed/input ratios with known
sizes: STAR0.0415, sort0.00154, MarkDuplicates0.0591, ISOFOX0.00296h/GiB.
The corrected duplicate-marking group has only one completed full case; older
successful optical-enabled tasks therefore remain a conservative reference,
not proof of corrected-mode scaling. Unknown-size observations inform memory
headroom and runtime floors but are not used to estimate a size rate.

Sort RAM is intentionally tied to its fixed64 threads/default buffer rather
than increasing with total file size: it spills to disk.80GiB exceeds every
observed batch RSS by about48%. No sort/compression settings are changed.
Picard retains8 CPUs; allocating64 does not turn it into a64-thread traversal.
The MarkDuplicates128GiB floor exceeds observed RSS by about37%. ISOFOX keeps
a160GiB floor because one completed full case used123GiB despite its small
runtime. Requests above768GiB/336h stop rather than silently clamping.

Outside the measured input ranges these rules are conservative extrapolations;
use the existing resource reports to review each completed execution. Duplicate
complexity, storage contention and hardware remain unmodelled variables.

## Wiring and recovery

`build_input_resources.py` emits the RNA task closures with the real native
process input names. `resources.json` records these as dynamic rules, without
pretending intermediate BAM sizes are known at launch. WGS-only packages omit
the RNA selectors. The full-native entry already loads the shared optical
comparison policy; use the current source rather than the older frozen queue
checkout for future packages.

For an existing failed run, preserve the original source, FASTQs, command,
session, DNA requests and work/cache. Add the RNA-only recovery overlay:

```bash
python3 workflow/prepare_run_recovery.py \
  --run /path/to/run --out /path/to/run/package_rna_resources \
  --route rna --rna-resources --session ORIGINAL_SESSION_UUID \
  --native-source /path/to/original/native/source
```

For the failed case's166.44GiB sorted BAM, MarkDuplicates resolves to
**8 CPUs,128GiB,22h**. Completed alignment/sort tasks must remain CACHED on
resume. Walltime directives alone do not change an analytical command; memory
affects the native Picard heap as before.

## Validation status

Python tests check workload scope, conservative retrospective envelope coverage,
fixed threaded CPU choices, Picard8 CPUs, limits and original-command recovery.
`tests/run_rna_resource_canary.sbatch` resolves all four actual variable shapes
through real Nextflow under Slurm, using sparse synthetic files and echo-only
tasks. It reads no sequence payload and runs no analytical tool. A passed canary
establishes directive resolution, not full-case speed or scientific equivalence.
Full-case recovery and actual allocations are recorded by the coordinator.

Private measured source table and Slurm snapshots:
`/global/scratch/hpc6049/rna_resource_sizing_20261123/measurements.csv`.
This directory retains the audit label from the preceding session; source
and recovery observations include November24 UTC.
