# Full-input native execution and initial resources

Site execution release: `input-resources/v0.1.0`, policy `input-v1`.
Scientific source: unmodified nf-core/oncoanalyser3.0.0,
`7c74c87a43749952b38c9a18915947570f0595a0`. Its `conf/base.config` supplies the
native labels; local modules receive `task.cpus` through their existing thread
arguments. This release changes preparation and site allocation only.

## Owner-selected route

November9: run complete KQ0205 tumour DNA, blood-normal DNA and tumour RNA
through FASTQ realignment using the native Hartwig reference and resource
bundle. The custom-reference BAM route was not selected. Do not use the
experimental ORANGE patch or override native report eligibility.

The two DNA FASTQ pairs are not present in the bounded local search. Prepare
all three libraries consistently from the available BAMs. These are complete
BAM-derived read sets, not a claim of original raw sequencer FASTQ identity.
No chromosome selection, downsampling, MAPQ, duplicate or QC-fail filtering.
Exclude secondary/supplementary records to emit each sequencing read once.
Prefer OQ where present, otherwise use BAM QUAL.

`workflow/prepare_full_fastqs.py` independently counts source primary records
with samtools flagstat, collates the BAM, exports all primary reads and verifies
every paired FASTQ name/count/record through gzip EOF, including CRC checks.
Require full primary-record accounting and no omitted singleton/other records
before writing `completion.json`. Preserve these records/files if the gate
fails; do not silently shrink the input scope. Source BAMs are unchanged.
Run bulk preparation with `workflow/prepare_full_fastqs.sbatch` under Slurm.

## Resource policy

`workflow/build_input_resources.py` reads the completed recovery records. One
unit is 250 million primary read records, rounded up. Use stage-relevant DNA
pair, single-library or RNA workload; per-library shared selectors conservatively
use the largest library they may process. Heavy threaded stages receive64 CPUs
above250 million reads, otherwise16. Picard MarkDuplicates receives8 because
reserving64 cores does not make its algorithm64-threaded.

For each stage the explicit versioned policy defines memory/time floors and
per-unit increments; the script is the exact coefficient table. Memory includes
a1GiB-per-requested-CPU allowance and rounds up to16GiB. Time rounds up to4h.
Requests above768GiB or14d fail planning instead of silently clamping. No
first-attempt retry ladder substitutes for this plan.

Evidence available before this full run:

- The chromosome-subset run completed39 tasks. At12 CPUs, tumour alignment took
  75m55s; REDUX37m48s; somatic SAGE37m7s; ESVEE6m4s with67.9GB peak RSS.
- Indexed full BAM records are4,594,276,610 tumour,1,894,142,243 normal and
  1,913,624,600 RNA. These include secondary/supplementary records and are
  **planning upper bounds**, not recovered primary-read counts.
- Those upper bounds give BWA64 CPUs/208GiB/40h, REDUX64/192/32h,
  somatic SAGE64/272/40h and ESVEE64/352/52h. The full pair is outside the
  subset's measured range. These allocations are conservative extrapolations,
  not validated runtime predictions or fitted coefficients. Exact allocations
  will be regenerated from actual recovered primary counts before Nextflow.
- The proposed time floors are deliberately much larger than naive linear
  extrapolation for steps such as ESVEE; a partial-chromosome fit failure and
  sparse downstream records cannot calibrate a successful whole-genome run.
  Memory floors exceed the subset peaks and include per-core and workload
  allowances. Collect full-run traces before claiming this model is validated.
- Current Slurm inventory includes64-CPU nodes with498000MB and larger CPU
  nodes. `sbatch --test-only` accepted64CPU/384GiB/72h and64CPU/192GiB/48h
  for `def-hpcg6049_cpu`, routed to `cpubase_14days`; no partition is forced.
- Scratch user accounting has no configured block limit; the filesystem showed
  about1.7PB free. Do not interpret that as a reservation or deletion authority.

FASTQ preparation uses64 CPUs/64GiB, with initial walltime rounded up to4h:
`max(12h, source_BAM_GiB / 8 * 1.5h)`. The1.5 factor is a declared safety margin
on a conservative8GiB/hour allowance for counting, collation, compression and
full pair validation. This is also an unvalidated first-full-case allocation.

## Execution package and entry

Freeze the site scripts/configuration, samplesheet and input bindings under:
`RUN_ROOT/<run_id>/{package,inputs,state,work,results,logs}`. Do not edit an
executing package. Keep reference/container assets at their existing locations.

Submit the three preparation jobs, then submit `workflow/run_full_native.sbatch`
with an `afterok` dependency on all three. The dependent head checks that the
samplesheet matches each completed FASTQ pair, generates
`state/input_resources.{config,json}`, and invokes the native launcher with
that configuration **last**, overriding the historical32-core site ceiling.
The native source remains3.0.0; it receives no process exclusions or manual
stage selection. DNA alignment and RNA alignment use the staged HMF genome,
BWA-MEM2 index and GENCODE38 STAR index.

Actual commands, SHA, job IDs and state belong in the external package's
`execution.json`. The prelaunch upper-bound table is replaced by actual measured
inputs for scheduling. Heavy steps remain CPU execution; no GPU is requested.
Native QC still governs whether ORANGE/PEACH can run. Pipeline completion alone
does not establish scientific validity or clinical report release.

## Verification

`tests/test_full_input_preparation.py` checks synthetic mapped, unmapped and
duplicate primary pairs, excluded secondary/supplementary records, original
quality recovery, workload growth, scope rejection and excessive-request refusal.
`tests/test_frontenac_launcher.py` covers Slurm spool configuration discovery.
Use `bash -n` on the two sbatch scripts and the native launcher. At execution,
retain generated config/evidence and inspect actual task `.command.run` CPU,
memory and walltime requests, then traces/accounting.
