# KQ0205 full original-FASTQ native execution

Owner instruction November10,2026: submit the full FASTQ route now using the
native Hartwig reference,64 CPUs for heavy steps and input-adjusted resources;
keep the custom-reference BAM run active.

Scientific source: stock nf-core/oncoanalyser3.0.0,
`7c74c87a43749952b38c9a18915947570f0595a0`. Site entry release:
`input-resources/v0.3.0`, policy `input-v2`.

## Inputs and reference

- Tumour DNA:8 original FASTQs,4 lanes with both mates,
  286,799,439,103bytes; flowcell from FASTQ header `HMNFNDSXF`.
- Matched blood-normal DNA:8 original FASTQs,4 lanes with both mates,
  150,532,895,224bytes; flowcell `HNMH2DSXF`.
- Tumour RNA:existing full BAM-derived pair,2 FASTQs,
  59,648,646,551bytes. No new RNA recovery is run and these are not represented
  as original sequencer FASTQs.
- Native CSV:8 separate DNA lane rows plus1 merged RNA row. The native workflow
  performs alignment/merging; no manual DNA lane concatenation or downsampling.
- Native genome:`GRCh38_hmf`, staged masked HMF FASTA and its matching BWA-MEM2,
  STAR/GENCODE38 and other indexes from the staged3.0.0 reference configuration.
  This run does not use d1.vd1 or the experimental ORANGE patch.

Runtime root:
`/global/scratch/hpc6049/wigits_runs/KQ0205_full_native_20261109_onco3.0.0/`.
The existing original-file transfer manifest/state provide the16-file DNA
source/landing scope. `state/original_fastq_scope.json` records the current
scope; `package_fastq_v0.3.0/samplesheet.csv` is the issued native CSV.
The historical `package/` with BAM-recovery paths is preserved and unused.

## Initial allocations

The existing input-v2 policy sizes by measured compressed bytes where no
independently verified read count is supplied; counts are recorded as unknown,
not invented. These are conservative allocation limits, not runtime predictions.

| Native process | CPUs | Memory(GiB) | Walltime(h) |
|---|---:|---:|---:|
| BWA-MEM2 DNA alignment |64|224|56|
| REDUX |64|224|44|
| SAGE somatic |64|336|64|
| SAGE germline |64|224|24|
| ESVEE |64|448|84|
| AMBER/COBALT |64|224|24|
| VIRUSBreakend |64|272|56|
| STAR/ISOFOX |64|160|12|
| RNA coordinate sort |64|128|8|
| Picard RNA MarkDuplicates |8|96|12|

The generated resource configuration is loaded last, overriding the historical
32CPU/180GiB site ceiling. Picard's serial-heavy operation is not represented
as64-threaded. Head allocation:4CPUs/8GiB/13days. Scheduler test-only admission
accepted64CPUs/448GiB/84h on`frnt137`; it did not submit an analysis job.

## Development checks

Slurm12444645:COMPLETED,exit0:0,7s; all14 tests passed with no skips,
including the direct-original-FASTQ entry test and existing resource-scaling,
mate-preparation and transfer tests. Login suite:14 tests,3 real-tool skips.
Shell syntax/diff checks passed. These are development checks, not added
whole-file patient audit stages.

## Submission

The frozen package was submitted as **Slurm head12444653**, observed RUNNING
on`frnt075`. Source release commit:
`26fd5d764dcda6335a8a6c941e6e6646d4301e14`.

Issued command shape:

```bash
sbatch --parsable --account=def-hpcg6049_cpu \
  --output=<run-root>/logs/native_full_fastq_head.%j.log \
  <run-root>/package_fastq_v0.3.0/workflow/run_full_native.sbatch \
  <run-root>/package_fastq_v0.3.0 KQ0205_full_native_20261109_onco3.0.0
```

`state/execution.json` records the exact absolute command, job and submission
time; `package_fastq_v0.3.0/release.json` pins source/reference/scope. The head
log confirms Nextflow26.04.4, Apptainer1.3.5, native3.0.0, mode`wgts`,
genome`GRCh38_hmf` and the generated resource config loaded last. The resource
file explicitly requests64CPUs for BWA-MEM2 and other heavy threaded tasks.
Task dispatch/allocation and completion are not established by head startup.

The
BAM-route controller12444300 and tumour preparation12444301 remain active
and are not dependencies of this separate full FASTQ run. Work/cache and output
directories are separate. Completion and report generation remain pending.

### Alignment dispatch confirmed

Native BWA-MEM2 jobs12444656/12444657/12444658 were observed RUNNING with
64CPUs/224GiB; STAR job12444659 was submitted requesting64CPUs/160GiB.
The first DNA task's `.command.run` contains`#SBATCH -c 64`,56h and229376MB,
and its `.command.sh` passes`-t 64` to`bwa-mem2 mem`. Thus the64-core request
is established in both Slurm allocation and the actual aligner command, not
only the configuration. Additional lane dispatch is ongoing. BAM jobs
12444300/12444301 remain RUNNING independently.
