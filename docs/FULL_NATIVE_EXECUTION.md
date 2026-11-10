# Full-input native execution and initial resources

Current staged-FASTQ entry: `input-resources/v0.3.0`, policy `input-v2`.
The earlier `input-resources/v0.1.0`/`input-v1` recovery-based execution below
is historical source, superseded for the current original-FASTQ run.
Scientific source: unmodified nf-core/oncoanalyser3.0.0,
`7c74c87a43749952b38c9a18915947570f0595a0`. Its `conf/base.config` supplies the
native labels; local modules receive `task.cpus` through their existing thread
arguments. This release changes preparation and site allocation only.

## Current original-FASTQ execution, November10

All16 original DNA FASTQs have landed through the recorded two-hop transfer:
four lanes with both mates for tumour and matched blood normal. The full
existing RNA FASTQ pair is also local. The native CSV contains8 DNA lane rows
plus1 RNA row; RNA is the existing full BAM-derived pair, not newly recovered
or claimed to be original sequencer FASTQs. DNA flowcell identifiers are taken
from the original FASTQ headers (`HMNFNDSXF`, `HNMH2DSXF`), not the run-folder
side prefix. Native `read_alignment_dna/main.nf` constructs per-lane read groups
and groups alignment products for subsequent sample processing.

The current entry `workflow/run_full_native.sbatch` accepts the frozen package
and run ID, sizes existing inputs with `build_input_resources.py --samplesheet`,
then calls stock Oncoanalyser3.0.0 with **GRCh38_hmf**. It does not invoke BAM
recovery or require recovery completion records for original sequencer files.
All lane rows remain explicit in the native samplesheet. The existing original
transfer manifest and completion state record the16-file scope.

Resources use policy input-v2, a declared compressed-byte workload proxy, not
invented read counts. Initial heavy threaded steps receive64 CPUs. Examples for
this complete case: BWA-MEM2=224GiB/56h; REDUX=224GiB/44h;
SAGE somatic=336GiB/64h; ESVEE=448GiB/84h; STAR=160GiB/12h.
Memory/walltime are conservative limits, not predicted durations. Picard RNA
MarkDuplicates retains its tool-appropriate8CPUs. The generated configuration
is loaded last, overriding the old32CPU/180GiB site ceiling. The head limit is
13days because the exact14day head request was previously rejected.

The full FASTQ and custom-reference BAM routes use separate run roots and
work/cache directories. Do not interrupt the running BAM route. Runtime
launch/accounting details are in
`tests/results/full-native-fastq-input/EXECUTION.md` once submitted.

## Historical setup and owner correction, November9

**Historical execution status, owner correction November9:** use existing full-size
FASTQs; do not recover them from BAMs. Preparation jobs12441461/12441462/12441463
were cancelled before execution (elapsed0); no full workflow head was submitted.
The prepared recovery-based entry below is preserved source, not the issued
execution route. The supplied inventory identifies the full DNA BAMs and an
existing full RNA FASTQ pair; the three FASTQ pairs prepared for the earlier
test contain only chr17/chr19-selected reads and their mates. Exact full DNA
FASTQ paths remain to be bound to the new samplesheet. Do not infer that the
owner's full FASTQs are absent from all storage, reuse subset files as full
inputs, or restart recovery without an explicit new instruction.

## Owner-selected route

### Existing FASTQ transfer, subsequent November9 instruction

Live RC2 listings established all16 original DNA FASTQs: four lanes and both
mates per tumour/normal library,437,332,334,327 bytes total. They remain under
the respective native run `fastq/` directories (runs0058 and0063). The earlier
claim of missing full DNA FASTQs did not establish NAS availability.

`workflow/transfer_fastqs_two_hop.py` copies one exact FASTQ at a time through
RC2→Fir→Frontenac. It adopts saved task IDs, requests checksum verification on
both hops, checks the landed size, and retains NAS and relay copies. It performs
no analysis, payload cleanup or BAM conversion. Runtime manifest/state:
`/global/scratch/hpc6049/wigits_runs/KQ0205_full_native_20261109_onco3.0.0/package_rc2_fastqs_v1/`
and `state/rc2_fastqs/`. Controller Slurm job12441547 adopts first-hop Globus
task548a63f5-dcb5-11f1-a6d6-0affd629e535. That task was observed ACTIVE/OK with
4,853,858,304 bytes transferred. Completion is pending, not implied by submission.

The existing full RNA FASTQ pair remains local under the recorded QUEEN001
`fastq_from_bam/KQ0205-RNA/` directory. A native full-run launch remains pending
the complete original DNA landing and an existing-FASTQ resource/input package.
Do not dispatch the recovery-based entry documented below for this run.

November9: run complete KQ0205 tumour DNA, blood-normal DNA and tumour RNA
through FASTQ realignment using the native Hartwig reference and resource
bundle. The custom-reference BAM route was not selected. Do not use the
experimental ORANGE patch or override native report eligibility.

The initial proposal was to prepare all three libraries from the available BAMs.
The owner rejected this proposal; the helper description below is historical.
Its outputs would be complete
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
