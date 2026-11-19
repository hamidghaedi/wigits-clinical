# Input-v4 measured alignment allowances

Release `input-resources/v0.6.0`, owner direction November19,2026: keep64CPUs
fixed; reduce memory with input size/headroom and tighten walltime using the
observed complete-task throughput. This supersedes v0.5.1's fixed224GiB RAM and
four-hour walltime floor. Preserve completed/running tasks; replace pending
alignment requests through the original native workflow/cache.

The preceding input-v2 planner summed every FASTQ lane into a library, selected
the larger tumour/normal library, then assigned that request to every BWA-MEM2
lane. It assigned KQ0254's274GiB whole-tumour request to each34GiB tumour lane
and each5GiB normal lane:64CPUs/224GiB/56h. This was the wrong workload scope.

## Implemented scope

| Step | Workload basis |
|---|---|
| BWA-MEM2 RAM/time | The current task's `reads_fwd.size() + reads_rev.size()`; individual lane/chunk pair; CPU64 fixed |
| REDUX | All input lanes of that specific sample, keyed by native `meta.sample_id` |
| Paired SAGE/ESVEE and other existing steps | Existing input-v2 per-stage scope and coefficients retained |

`build_input_resources.py` emits the task-local closures from
`conf/lane_alignment_resources.config` and individual REDUX sample allocations.
Its JSON retains every lane's size and preview allocation beside complete
library and paired-case plans. File size is a workload proxy, never a read count.
No full payload scan, hashing, lane filtering or scientific tool change.

## Alignment policy and evidence

Let G be compressed R1+R2 GiB for the actual task:

- CPUs:64, fixed for every alignment task.
- Memory:48+2G GiB, rounded up to16GiB. The48GiB base retains reference,
  sorting and safety allowances; workload growth gives larger tumour lanes
  more RAM than the recently measured normal lanes.
- Estimated pairs:compressed input bytes/120. This is explicitly a workload
  estimate, never a measured read count. The32normal lanes had142.96–148.15
  compressed bytes/pair; using120allows roughly19%more pairs than the lowest
  observed ratio. This is not a guarantee for every library type.
- Time:estimate/4million pairs per minute, multiplied by2, plus15minutes;
  round up to whole hours with a1hour floor. The24completed64CPU normal jobs
  sustained4.50million pairs/min pooled (individual sample groups4.27–4.59).
  Use the complete-task rate including startup/BAM conversion/sorting, not the
  faster BWA-only rate. Queue time is excluded from the measured throughput.
- Refuse planning/evaluation above768GiB or336h; never silently clamp an
  excessive request and rely on retries.

The new32-lane CPU/RAM/throughput assessment contains completed normal lanes
only: KQ0254 at16CPUs and KQ0322/KQ0323/KQ0400 at64. Peak batch RSS was
30.39–45.71GiB. Different inputs and hardware make this observational evidence,
not a controlled scaling benchmark. Tumour needs are not established by those
normal-only peaks. The source assessment/throughput files are retained in the
coordinator's private resource audit dated20261119T092747-0500.

Additional historical memory context:32completed64CPU lanes from KQ0205, KQ0004, KQ0005 and
KQ0012 took0.69–2.51h. Slurm batch RSS observations span67.68–144.52GiB;
these are step observations, not a guaranteed simultaneous whole-job maximum.
Average busy CPUs span41.57–44.90. KQ0205's four tumour lanes represented
roughly67GiB compressed per lane and reached144.52GiB RSS. The new RAM formula
would allocate192GiB at67GiB, retaining margin above that historical observation.
The formula is a conservative envelope, not a fitted regression or universally
validated minimum RAM. The32recent normal jobs and historical memory context
are different cohorts and are not pooled into a new throughput benchmark.

Current pending tumour previews (all64CPUs): KQ0254→128GiB/3h;
KQ0322/KQ0323→144GiB/4h; KQ0400→128–144GiB/4h, according to each pair's size.
KQ0254 normal4.78–5.02GiB would request64GiB/1h for a new alignment, but
its completed normal alignments are reused. REDUX remains whole-sample:
normal39.23GiB →64CPU/144GiB/16h; tumour274.25GiB →64CPU/224GiB/44h.
Native analytical arguments, reference, all lanes/mates and downstream scope
are unchanged. The existing coordinator's staged-native queue is reused.

## Validation and execution

Python regressions cover lane scope, adding an unrelated/larger lane without
inflating another, fixed64CPUs, estimated-count labelling, maximum refusal and
per-sample REDUX. Explicit completed-task CPU history preserves the original
command hashes for cached16CPU lanes; it never lowers CPU for a new task.
Resume packages supply `completed_alignment_cpus.json` only for completed
source tasks, and the entry passes it through `--completed-alignment-cpus`.
`tests/test_lane_resource_config.nf` is a real Nextflow metadata-only canary:
its synthetic sparse files have known sizes; tasks only echo resolved CPU,
memory and time. A test-only local executor capacity allows these directives
to be inspected without allocating the requested memory or running tools.
Run the canary head under Slurm, with4CPUs/8GiB and no real patient inputs.
Real-case run/restart IDs and canary outcomes are recorded in
`tests/results/lane-resource-sizing/README.md` and the coordinator's private
`phase2_native_20261118` queue. Old packages, work/cache and all input files
remain preserved; new packages record the new resource release and use resume.
