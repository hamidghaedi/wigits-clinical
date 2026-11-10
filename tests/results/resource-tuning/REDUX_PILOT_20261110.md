# REDUX 64-CPU partitioning pilot

10 November 2026. Native tool: REDUX 2.0.5, hmftools
`8ec3714a248c01d49a499a684082d82a95d45e14`; native task source:
nf-core/oncoanalyser 3.0.0 `7c74c87a43749952b38c9a18915947570f0595a0`.
Benchmark source: `workflow/benchmark_redux_performance.py` and its sibling
Slurm script. Both full-case workflows and production allocations are unchanged.

## Performance result

Slurm 12446613 completed 0:0 in 17m43s on frnt119 with 64 CPUs and 192GiB RAM.
All four sequential trials used the same 4.87GB pre-REDUX normal BAM from the
chr17/chr19 fixture, reference, 120GiB Java heap, tool/container and scientific
flags. Only the native partitioning/concatenation controls differed.

| Trial | Elapsed seconds | Average busy CPUs | CPU allocation use | Speed versus initial baseline |
|---|---:|---:|---:|---:|
| Default ratio2 | 378.55 | 5.49 | 8.57% | 1.00x |
| Ratio8 | 159.46 | 13.78 | 21.54% | 2.37x |
| Ratio8 with parallel concatenation | 160.70 | 14.34 | 22.41% | 2.36x |
| Repeated default ratio2 | 363.89 | 5.82 | 9.09% | 1.04x |

Per-trial CPU is the wrapper's child accounting, not separately scheduled Slurm
steps. Whole-job Slurm TotalCPU is 2h24m59s; batch MaxRSS is 30.20GiB. The ratio8
trial remains 2.28x faster against the repeated baseline. Parallel concatenation
adds no measured benefit here; a single trial cannot establish a small difference.

The restricted workload strongly favors finer partitions: default 128 partitions
have 7 partitions exceeding 1million reads, holding 99.12% of processed reads.
Ratio8 has 512 partitions with 24 over 1million, holding 98.91%. The parallel mode
has 504 partitions with 24 over 1million. Other partitions contain cross-chromosome
mates and are not literally empty. This explains how the fixture exposes idle
workers, but the 2.3x gain cannot be extrapolated to whole-genome data or the
full case's problematic final partition.

There is a more specific limitation in the pinned implementation:
`redux/src/main/java/com/hartwig/hmftools/redux/PartitionThread.java:123–141`
balances ordinary chromosomes by base length, excludes alt-contig lengths from
the balancing calculation and appends those contigs to the final partition.
This matches the full tumour log's final chrY/chrM/3344-region partition, which
finished 29minutes after 127/128 partition readers. A larger partition ratio does
not distribute that alt/decoy tail across workers, so it is not established as
the remedy for the actual full-case bottleneck.

## Equivalence findings

Successful process exit and faster execution do not establish equivalent results.
All trials report 91,847,076 reads and identical native unmapping statistics.
Jitter parameters and the decompressed microsatellite table are identical.
However, small differences are present:

| Comparison to default | Duplicate reads | Duplicate groups | ALIGNMENT_ONLY consensus groups | Changed BQR count rows | Changed recalibrated-quality rows | Largest absolute quality change |
|---|---:|---:|---:|---:|---:|---:|
| Default values | 7,907,227 | 7,099,769 | 6,776,431 | — | — | — |
| Ratio8 | 7,907,226 | 7,099,768 | 6,776,430 | 958/1536 | 210/1536 | 0.18 |
| Ratio8 plus parallel concatenation | 7,907,227 | 7,099,769 | 6,776,431 | 1016/1536 | 255/1536 | 0.15 |
| Repeated default | 7,907,227 | 7,099,769 | 6,776,431 | 0/1536 | 0/1536 | 0.00 |

The BQR table keys are identical across trials. Absolute sums of BQR count
differences are 115,887 and153,815, against the baseline 10,803,709,447 total
counts. These are real numerical differences, not simply row order or plot
metadata. The duplicate-frequency TSV also has some different high-count bin
labels even in the repeated baseline; do not treat every TSV byte difference
as a new effect of partitioning. No downstream calling comparison has been run,
so neither biological insignificance nor a variant-call impact is established.

Order-agnostic record/flag/quality/position/CIGAR/mate/all-auxiliary checks
completed under Slurm 12446759, using samtools 1.22.1: COMPLETED 0:0 in 57seconds,
16 CPUs/32GiB. The explicit
command avoids `checksum -a` because its format-roundtrip defaults sanitize and
exclude some tags. These are isolated experiment comparisons, not new routine
patient input gates. Results:

| Output | All alignment records, including secondary/supplementary | Compared to baseline |
|---|---:|---|
| Baseline | 98,322,812 | Reference |
| Ratio8 | 98,322,811 | One fewer output record; content checksums differ |
| Ratio8 plus parallel concatenation | 98,322,812 | All reported order-agnostic content checksums identical |
| Repeated baseline | 98,322,812 | All reported order-agnostic content checksums identical |

The checksum comparison includes QNAME, all flag bits, sequence, quality,
position, CIGAR/MAPQ, mate fields and all auxiliary tags, without sanitization.
It excludes headers and record order. Equal checksums are evidence for the
compared data, not a cryptographic proof of scientific equivalence; unequal
checksums establish a difference.

Slurm 12446767 localized the ratio8 difference: COMPLETED 0:0 in 3m18s with
24 CPUs/64GiB. The absent record is one REDUX-generated consensus mate, identified
by its consensus-name suffix and `CR` auxiliary tag. Its alignment crosses a new
partition boundary; the other mate starts on the other side and is retained.
Baseline has both records for that consensus name; ratio8 has one. This supports
a partition-boundary effect, rather than establishing that an original sequencer
read was dropped. Exact read name and coordinates remain private. Original
family handling and downstream impact are not yet characterized.

The initial localization helper used an invalid GNU sort option and failed to
propagate its background error despite a scheduler exit0 (job12446764). Its
empty output was rejected, not treated as an equivalence result. The corrected
helper uses `--parallel=6`, explicitly waits on both pipeline processes and writes
to a separate directory; all localization findings above come from12446767.

The combination mode retains baseline BAM content by these measurements but
still changes the BQR sidecar. That distinction matters because downstream
calling can consume the sidecar; no downstream impact has yet been measured.

**Disposition: promising speed improvement; no production adoption.** Resolve
and characterize record/BQR differences first, then confirm a promising setting
on a representative full input. No CPU cap is lowered and no native scientific
option is disabled to improve the performance percentage.

## Evidence location

```text
/global/scratch/hpc6049/wigits_runs/resource_tuning_dev_20261110/performance_pilot/
  execution.json
  redux.12446613.log
  redux_trials/measurements.json
  redux_trials/<mode>/command.sh
  redux_trials/<mode>/command.log
  redux_trials/<mode>/redux_KQ0205-gDNA/
  equivalence.sbatch
  equivalence.12446759.log
  equivalence/{bqr_comparison,partition_workload,checksums}.json
  equivalence/<mode>.checksum.tsv
  locate_record_difference_v2.sbatch
  read_difference.12446767.log
  equivalence/read_names_v2/
```

Private BAMs, per-read information and full accounting remain on scratch.

## Verification and remaining work

This increment changes documentation only; benchmark source was already tested
at development commit `bf9a704`. `git diff --check` passes. Both full-case
controllers remain untouched. The observed settings are not installed into a
future resource policy or scientific workflow.

Next: resolve the consensus-boundary behavior and BQR changes before considering
a production setting; a full-genome confirmation would then be needed because
the alt-contig tail and sparse fixture partition layout differ. No full-input
performance rerun or downstream caller comparison was submitted in this increment.
