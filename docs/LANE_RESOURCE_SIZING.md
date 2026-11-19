# Input-v3 lane and sample resource sizing

Release `input-resources/v0.5.0`, owner direction November19,2026.

The preceding input-v2 planner summed every FASTQ lane into a library, selected
the larger tumour/normal library, then assigned that request to every BWA-MEM2
lane. It assigned KQ0254's274GiB whole-tumour request to each34GiB tumour lane
and each5GiB normal lane:64CPUs/224GiB/56h. This was the wrong workload scope.

## Implemented scope

| Step | Workload basis |
|---|---|
| BWA-MEM2 | The current task's `reads_fwd.size() + reads_rev.size()`; individual lane pair, or actual chunk pair if native splitting is enabled |
| REDUX | All input lanes of that specific sample, keyed by native `meta.sample_id` |
| Paired SAGE/ESVEE and other existing steps | Existing input-v2 per-stage scope and coefficients retained |

`build_input_resources.py` emits the task-local closures from
`conf/lane_alignment_resources.config` and individual REDUX sample allocations.
Its JSON retains every lane's size and preview allocation beside complete
library and paired-case plans. File size is a workload proxy, never a read count.
No full payload scan, hashing, lane filtering or scientific tool change.

## Alignment policy and evidence

Let G be compressed R1+R2 GiB for the actual task:

- CPUs:16 for G≤10;64 for G>10, retaining the previous workload-band rule at
  the correct task scope.
- Memory:64GiB fixed index/sort allowance +2G GiB +0.5GiB per CPU, rounded up
  to32GiB. Minimum is96GiB for nonempty input under this rule.
- Time:max(4h, G/8 ×64/CPUs hours), rounded up to2h. This is an8GiB/hour
  allowance at64CPUs, with a conservative linear CPU allowance.
- Refuse planning/evaluation above768GiB or336h; never silently clamp an
  excessive request and rely on retries.

Thirty-two completed64CPU lane alignments from KQ0205, KQ0004, KQ0005 and
KQ0012 took0.69–2.51h. Slurm batch RSS observations span67.68–144.52GiB;
these are step observations, not a guaranteed simultaneous whole-job maximum.
Average busy CPUs span41.57–44.90. For context, KQ0205's four tumour lanes
represent roughly67GiB compressed per lane and took~2.5h; the new8GiB/hour
allowance is substantially slower than that measured throughput. The RAM
formula retains an index/sort floor and workload/thread allowances rather than
dividing full-sample RAM by lane count. These are conservative allocations with
explicit margins, not fitted runtime predictions or universally validated limits.

KQ0254 previews: normal4.78–5.02GiB →16CPU/96GiB/4h;
tumour34.14–34.64GiB →64CPU/192GiB/6h. Its REDUX remains whole-sample:
normal39.23GiB →64CPU/144GiB/16h; tumour274.25GiB →64CPU/224GiB/44h.
Native analytical arguments, reference, all lanes/mates and downstream scope
are unchanged. The existing coordinator's staged-native queue is reused.

## Validation and execution

Python regressions cover lane scope, adding an unrelated/larger lane without
inflating another, workload bands, maximum refusal and per-sample REDUX.
`tests/test_lane_resource_config.nf` is a real Nextflow metadata-only canary:
its synthetic sparse files have known sizes; tasks only echo resolved CPU,
memory and time. A test-only local executor capacity allows these directives
to be inspected without allocating the requested memory or running tools.
Run the canary head under Slurm, with4CPUs/8GiB and no real patient inputs.
Real-case run/restart IDs and canary outcomes are recorded in
`tests/results/lane-resource-sizing/README.md` and the coordinator's private
`phase2_native_20261118` queue. Old packages, work/cache and all input files
remain preserved; new packages record the new resource release and use resume.
