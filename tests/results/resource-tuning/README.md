# CPU performance investigation,10 November2026

Branch: `hartwig_adoptation`. Owner clarified that reducing reservations is not
a substitute for making execution faster. Proposed per-tool CPU caps and the
lower-resource BAM dispatcher were withdrawn before commit/release. Current
`build_input_resources.py`, `bam_intake.py` and BAM entry match the prior release;
input-v2 and BAM entry0.3.0 remain in force. There is no tuning release yet.

## Completed

- Inspected both full-case resource CSVs, native task commands and phase logs.
- Retrieved exact hmftools source for REDUX2.0.5 and ISOFOX2.0.1. Commit IDs,
  source functions, measured CPU use and limitations are in
  [the investigation](../../../docs/RESOURCE_TUNING.md).
- Added a separate REDUX performance pilot. No native checkout or live package
  edits; no scientific parameter changes, input transfers or production hooks.
- `python3 -m unittest discover -s tests`:24 tests, exit0,3 real-tool skips on
  login. Final source snapshot on Slurm12446615:24 tests pass with **no skips**,
  COMPLETED0:0,5seconds. Includes the existing synthetic mate/quality tests and
  new checks that trial commands retain scientific arguments and64 CPUs.
- `bash -n workflow/benchmark_redux_performance.sbatch`, exit0;
  `git diff --check`, exit0.
- Verified32 frozen active package files still matched the saved pre-change
  fingerprints. This is development evidence, not a new patient input gate.

## Running / pending

Slurm **12446613**,64 CPUs/192GiB/2h, running onfrnt119 when recorded.
It executes four trials sequentially on one node: baseline, partition ratio8,
ratio8 plus parallel concatenation, repeated baseline. All trials use the same
native pre-REDUX chromosome-test normal BAM, reference, container,120GiB heap
and scientific options. A repeated baseline helps expose cache/order effects.

Frozen pilot source, exact submission script/arguments, logs and measurements:

```text
/global/scratch/hpc6049/wigits_runs/resource_tuning_dev_20261110/performance_pilot/
  benchmark_redux_performance.py
  benchmark_redux_performance.sbatch
  redux.12446613.log
  redux_trials/{baseline,partition8,partition8_parallel,baseline_repeat}/
  redux_trials/measurements.json
```

Pinned sources, private accounting and earlier unissued allocation candidates
are under the parent workspace's `source/` and `evidence/`. Final source tests
are in `final_source_tests/`. Earlier candidate tests are superseded evidence,
not proof that a CPU cap was adopted. Neither full-case controller was cancelled,
requeued or resubmitted; both remained running at the final inspection.

Pending: pilot timing/CPU comparison, output scientific equivalence, a matched
full-input confirmation if promising, and separate sorting/streaming experiments.
No speedup or general scalability improvement has yet been established.
