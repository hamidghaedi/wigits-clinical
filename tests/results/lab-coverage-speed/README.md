# Coverage speed regression — November25

The earlier maintained DNA module uses Parabricks bammetrics4.6.0-1 with32 CPU
threads/64GiB and input-sized walltime. The historical shell worker requested
32CPU/250GiB. The new native QC route instead uses Picard3.1.0-SNAPSHOT,
4CPU/16GiB and maxForks4. Live tasks consume approximately one CPU and1–1.5GiB;
completed tumour coverage tasks took6h38–9h55. Increasing the Picard allocation
alone does not create parallel execution. Native Picard's exact source identity
was selected for standards measurement, but its performance cost was not tested
against the previous accelerated path before deployment.

The separate COVERAGE_BENCHMARK entry reuses the existing installed Parabricks
image on one completed normal and tumour. CPU-only,32 threads/64GiB, up to2 tasks;
walltime uses the old2min/GiB model with2h floor. No GPU, realignment or source
BAM modification. Original Picard jobs, completed metrics, final report sources
and future production defaults remain unchanged during the comparison.

Use the identical native BAM/index, HMF FASTA/FAI, existing complete-dictionary
primary-chromosome interval list, MAPQ20/baseQ20/cap250 and default unpaired=false.
Compare every reported metric and depth histogram numerically without tolerance,
silent averaging, or a new scientific threshold. Record differences; do not call
Parabricks execution native Picard. Source/interval/filter alignment is required;
a format-compatible table alone does not establish interchangeable measurements.

NVIDIA documents bammetrics as accelerated CollectWgsMetrics with a num-threads
option and the matched quality/cap defaults:
https://docs.nvidia.com/clara/parabricks/4.6.0/documentation/tooldocs/man_bammetrics.html
The project source is `dnaseq_nf/modules/parabricks_bammetrics.nf` and the older
coordinator `scripts/workers/bammetrics_worker.sh`; these are read-only references,
not imported DNA workflows or replacements for new report packages.

Offline suite125 tests,3 existing skips; adapter suite3 pass. Comparison tests
cover numerical formatting, real metric/histogram differences and missing data.
Shell syntax and diff checks pass. Full-size timing and numerical comparison are
pending execution; no speedup or equivalence is claimed before results exist.

Direct adversarial review: preserve full source scope and exact reference;
no sampling, truncated reads, extra raw-data checks or increased coverage cap.
The comparison is isolated and writes only new outputs. Do not adopt a backend
on vendor claims alone or change current metrics to manufacture agreement.
Changing the future measurement producer requires documented results and actual
source labels; the existing native-Picard production route stays bound meanwhile.
