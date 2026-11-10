# Diagnose speed before changing allocations

Owner clarification,10 November2026: prioritize faster computation and investigate
missing tool settings; lowering CPU caps alone is not the solution. Earlier
proposed caps have been withdrawn from executable source. Native sizing remains
**input-v2** and BAM preparation remains **v0.3.0**. Neither running full-case
package nor its Slurm allocation is changed. No resource-tuning release is issued.

## Pinned source and observations

- nf-core/oncoanalyser3.0.0: `7c74c87a43749952b38c9a18915947570f0595a0`,
  `modules/local/redux/main.nf`, `modules/local/isofox/run/main.nf`,
  `modules/nf-core/samtools/sort/main.nf`, `conf/modules.config`.
- REDUX2.0.5: hmftools commit
  [`8ec3714a248c01d49a499a684082d82a95d45e14`](https://github.com/hartwigmedical/hmftools/tree/8ec3714a248c01d49a499a684082d82a95d45e14/redux).
- ISOFOX2.0.1: hmftools commit
  [`ecd6b49792f604c9d6859e2db5a8cee9c4fffba3`](https://github.com/hartwigmedical/hmftools/tree/ecd6b49792f604c9d6859e2db5a8cee9c4fffba3/isofox).

Inspected task commands contain the requested thread arguments. A missing
`-threads64` is not the explanation for the observed low use.

| Step | Average busy CPUs out of64 | Diagnosis and next action |
|---|---:|---|
| BWA-MEM2,8 completed DNA lanes | 41.57–44.01 | Substantial parallel use. Inspect sort/compression and storage phases before claiming alignment can saturate64 cores. |
| REDUX normal/tumour, completed FASTQ route | 16.94 /11.36 | Tumour log has127/128 partitions finished29minutes before the last. Final concatenation/indexing then takes about14minutes. Native partition/overlap settings are candidates. |
| RNA samtools sort, completed | 5.31 | Native command has64 threads but no explicit `-m`; the tool's memory-per-thread default applies. Benchmark sort buffers, compression and storage placement. Disk volume alone is not proof of I/O wait. |
| ISOFOX, completed | 11.14 | Source creates chromosome tasks in `Isofox.processBam`; chromosome1 finishes about13minutes after chromosome3 and7minutes after chromosome2 in this run. No finer task-size option was found in `IsofoxConfig`. Do not split BAMs and concatenate expression/fusion outputs naively. |
| AMBER/COBALT, completed | 19.19 /26.14 | Thread arguments are present. Profiling/segmentation phases need separate investigation; no lower cap or new setting is enabled. |
| SAGE germline, completed | 3.60 | Native `-germline -panel_only` bounds its workload; elapsed80seconds. No scientific scope expansion merely to consume CPUs. |
| CIDER DNA/RNA, completed | 2.49 /1.47 | Targeted analysis completes in142/85seconds. Measure startup versus analysis before pursuing more threading. |
| SAGE somatic, live snapshot | 61.65 | About96% CPU use; retain64 CPUs. |
| ESVEE, live snapshot | 23.01 | Multi-phase run incomplete; observed438.54/448GiB RAM. Do not reduce RAM or infer later-stage parallelism from this snapshot. |
| VirusBreakend, live snapshot | 4.30 | Log in initial viral classification/taxonomy phase. Inspect later phases before changing resources. |
| BAM mate preparation, completed | 3.78 |5h54m total: name-sort startup about55minutes, Picard traversal3h19m, remaining coordinate sort/index about1h40m. No thread argument is supplied to Picard's record traversal; allocating more CPUs does not itself parallelize it. Investigate streaming and sort buffers, then a separately proven equivalent alternative if needed. |

Figures come from completed parent Slurm TotalCPU/elapsed or labelled provisional
`sstat` batch accounting. MaxRSS is a batch observation, not a guaranteed whole-job
peak. This is one real case, not a general scaling calibration.

## REDUX controls found in the actual installed version

`ReduxConfig.registerConfig` registers `partition_ratio` (default2) and
`parallel_concat` (off by default). `ReduxApplication.createPartitionThreads`
uses the ratio times the worker count to build a shared partition queue. Finer
partitions might reduce a long tail; more partitions also increase overhead.
`parallel_concat` starts `FinalBamWriter` while partition work proceeds, reserving
one worker. This can overlap final writing but may increase I/O or spend time
waiting. Neither is guaranteed faster. They are absent from the inspected live
commands. The module accepts these through its existing `task.ext.args`; native
tool code does not need patching to test them.

## Isolated64-CPU pilot

Maintained scripts:

- `workflow/benchmark_redux_performance.py`
- `workflow/benchmark_redux_performance.sbatch`
- `tests/test_redux_performance.py`

The pilot replays an existing native REDUX task into new directories. Every
trial gets64 CPUs, the same120GiB heap inside192GiB allocated RAM, the same
container/reference/input BAM and scientific flags. Modes:

1. Native default partition ratio2 and no parallel concatenation.
2. `-partition_ratio8` (512 partitions at64 workers).
3. `-partition_ratio8 -parallel_concat` (63 workers,504 partitions).
4. Repeat baseline to expose cache/order effects.

The input is the existing4.87GB pre-REDUX normal BAM from the chr17/chr19 fixture;
no genomic transfer or new caller run is needed. Trial order shares one node.
`measurements.json` records elapsed time and child CPU accounting; final Slurm
accounting describes the entire pilot. Outputs and logs remain on scratch.

This pilot can identify overhead or promising settings. It cannot establish
improvement on the full case's problematic final partition. Before adopting a
setting, inspect output read/flag/tag and native QC/statistical equivalence,
then test a promising candidate against baseline on a representative full input.
Do not alter consensus, filters, read scope, BQR, jitter or reference to gain speed.

The next bounded experiments are sorting buffers/compression and overlapping
Picard output with coordinate sorting. No Picard replacement, named-pipe adapter,
node-local staging or production setting is implemented here. Node-local storage
requires enough capacity for the real spill volume; it is not a universal fix.

Execution, test evidence and pending results are in
[the verification record](../tests/results/resource-tuning/README.md).
