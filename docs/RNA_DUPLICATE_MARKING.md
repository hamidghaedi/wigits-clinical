# RNA duplicate marking

Release: `rna-duplicates/v0.1.0`. Native Oncoanalyser source, GATK 4.6.1.0,
alignments, caller settings and the input-v4 resource planner are unchanged.

## Problem and correction

An RNA MarkDuplicates task exhausted two 12-hour allocations while repeatedly
comparing optical coordinates inside very large duplicate sets. Logs showed
sets of 67,764 paired-read ends and continuous OpticalDuplicateFinder work.
The second attempt reproduced the problem. More alignment CPUs cannot make
this duplicate-set comparison a 64-thread algorithm.

Picard explicitly recommends `READ_NAME_REGEX=null` for RNA-seq or other data
with extremely large duplicate sets when optical-adjusted library complexity
is not the goal: [the official MarkDuplicates reference](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates).

`workflow/conf/rna_duplicate_marking.config` passes that setting only to
`READ_ALIGNMENT_RNA:GATK4_MARKDUPLICATES`. Ordinary duplicate marking remains
enabled, reads are retained, and the output BAM/index/metrics interfaces are
unchanged. DNA REDUX and all DNA/RNA downstream analyses remain enabled.
Optical-duplicate counts and the optical-adjusted library-size estimate are
not assessed by this RNA mode. A zero in that auxiliary metrics column must
not be described as a measured absence of optical duplicates.

New packages using the released `run_full_native.sbatch` load the shared RNA
configuration before explicit package overrides. Already executing frozen
packages retain their source. No other resource policy is changed here.

## Resume an existing failed RNA stage

Use the original command, source, samplesheet, resource configuration, work
directory and session. The maintained recovery preparer adds one final config;
it does not regenerate resources, reconstruct FASTQs or re-enter alignment:

```bash
python3 workflow/prepare_run_recovery.py \
  --route rna --run /path/to/existing/run \
  --session ORIGINAL_SESSION_UUID \
  --native-source /path/to/recorded/native/source \
  --out /path/to/existing/run/package_rna_recovery

sbatch --account=YOUR_CPU_ACCOUNT \
  /path/to/existing/run/package_rna_recovery/resume_completed_run.sbatch \
  /path/to/existing/run/package_rna_recovery
```

An explicit local source must equal the source in the saved launch. The
recovery manifest records the original command and input-sheet hash. Check
the new trace: completed STAR and sorting tasks must be CACHED; failed or
previously unstarted tasks may execute. Do not count submission as completion.

## Validation

`workflow/compare_rna_duplicate_marking.py`, dispatched by
`workflow/test_rna_duplicate_marking.sbatch`, runs the installed GATK container
twice on 30 synthetic coordinate-sorted records. The test compares all first
11 SAM fields, including ordinary duplicate flags, and the examined-pair,
duplicate-pair and percent-duplication metrics. Only optical attribution is
expected to differ. No real input rows enter the fixture.

See [the runtime validation record](../tests/results/rna-duplicate-marking/README.md).
This is an interface/flag regression test, not clinical validation of a real
case or a whole-cohort performance guarantee.
