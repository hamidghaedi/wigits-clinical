# Distributed examples and measured-count adapter

Owner-authorized experiment, 11 November 2026. Development remains on
`hartwig_adoptation`. This does not modify active native workflows or the
existing DeepSomatic inference/postprocessing jobs, and changes no default
caller, scientific filter or downstream feature.

## Example-generation experiment

Entry: `workflow/benchmark_deepsomatic.py`, component 0.1.0, with
`workflow/benchmark_deepsomatic.sbatch`. The completed baseline is the full
KQ0205 native-reference REDUX run: 64 CPUs, 64 shards, 19h47m42s Slurm elapsed,
104 GiB peak batch RSS. No chromosome restriction or downsampling is introduced.

The experiment replays the exact native make_examples_somatic arguments from
the completed baseline. Only output paths and native runtime TSV collection
change. Slurm starts two tasks on two distinct nodes, with 32 CPUs per task.
Rank 0 owns even shard IDs, rank 1 odd shard IDs; together they own exactly
0–63. This keeps 64 total concurrent workers and the original 64-shard output
namespace. Small-model output shards are retained alongside main-model examples.
No inference or postprocessing experiment is automatically submitted.

Allocation: 128 GiB per node, 40h initial limit. Time is twice the measured
same-input baseline rounded to four hours; memory per node already exceeds
the measured total single-node baseline peak. This does not assume linear
speedup. Two-node placement can improve memory/I/O availability; it cannot
halve a purely CPU-bound workload with unchanged total CPUs. Queue time,
hardware, filesystem cache/load and instrumentation differ from the baseline.
Record those limits rather than attributing every timing difference to topology.

The dependent comparison uses four CPUs and 16 GiB to stream-decompress and
hash all 64 main-model shards, 64 small-model shards and 64 example-info JSONs
against the baseline. Gzip timestamps are ignored. Exact decompressed equality
is strong replay evidence; unequal files need semantic inspection and are not
automatically evidence of a scientific difference. Original files are never
deleted. About 1.05 TB of additional compressed examples is expected from the
measured baseline. Local GPFS reports no enforced user scratch quota and ample
filesystem space; the generic diskusage_report fails to parse this site's quota
output, so direct mmlsquota was used.

## Adapter increment

Entries: `workflow/deepsomatic_adapter.py`, `run_deepsomatic_adapter.py` and
`deepsomatic_adapter.sbatch`, component 0.1.3. Explicit external config binds the VCF, sample
names, BAMs, reference, native SAGE output, resource directory and container
paths. It does not scan for inputs or choose tumour/normal identities.

1. Retain the complete original caller VCF. In component 0.1.1, explicitly
   exclude only `RefCall` and `GERMLINE` from the somatic annotation derivative;
   count these exclusions by original FILTER. Retain `NoCall` and every other
   class unchanged, with no VAF or quality threshold. Germline analysis remains
   the separate native branch. The pilot has 651526 total records but only
   15279 outside RefCall/GERMLINE, so annotating all reference/genotype records
   would add substantial irrelevant work. This scope never promotes a NoCall.
2. Add the normal genotype column using actual DeepSomatic NAD/NDP. Leave
   normal GT missing. Derive normal AF from NAD/NDP; missing measurements stay
   missing, zero depth does not become a clean normal genotype. Preserve all
   retained tumour fields and FILTER values. Add tumour AF from VAF.
3. Retain the original NAF vector without relying on its inconsistent native
   Number=R header. The derivative declares Number=. and explicitly says it
   remains in original allele order. It is not a projected normal AF field.
4. Use native bcftools `norm -m -any --old-rec-tag DS_ORIGINAL_ALLELE` for
   multiallelic projection. Do not perform reference realignment, rewrite
   failed/germline calls as PASS, or restrict the derivative to SAGE overlap.
5. Apply SAGE 5.0.2 ReferenceData/TierSelector and reference-only RepeatInfo
   search through existing helpers ported from DNA adapter commit `7093033`.
   Add the reference-centred trinucleotide as native VariantReadContext does.
   Original DNA source remains untouched. Native high-confidence, hotspot,
   driver-gene and Ensembl resources are the oncoanalyser 3.0 bundle.
6. Check that native annotation preserves alleles, FILTER, QUAL values and
   genotype fields. Audit exact allele matches against the existing native
   SAGE output, recording measured RC_INFO/RC_CNT availability, LPS presence,
   tier agreement and separate caller counts without copying classifications.
   Take only INFO additions from the Java output and keep the exact normalized
   source payload to avoid treating VCF serialization as new measurements.
   Multiple native contexts for one allele are recorded, not chosen or treated
   as an annotation-job failure. All requested allele rows remain in the audit.

This is an **intermediate annotation product**, not a released `sage_dir`.
No pipeline consumers are disabled. PAVE population/PON/effect annotation and
PURPLE/downstream execution remain native and await the context increment.
Microhomology in native SAGE is derived from selected read context, sometimes
augmented by reference context; it is not silently fabricated from flanks here.
Missing authentic RC_INFO, RC_CNT, phasing and contamination-recovery semantics
are explicit open work. No changed annotations are presented as validated calls.

## Native context probe

`workflow/probe_sage_context.py` runs unmodified SAGE 5.0.2 on up to 1000
PASS alleles from the completed pilot using `specific_var_file`. The native
CSV schema is `Chromosome,Position,Ref,Alt`. Existing BAMs, resources, candidate
filters and TINC flags are retained; only the queried scope, threads, heap and
output location change. This tests whether genuine evidence can be obtained;
it does not force alleles, lower thresholds or promise all will be emitted.
Native SAGE outputs remain separate from DeepSomatic outputs. Contamination
estimates and phase sets from a restricted probe are not sample-wide results.

Probe component 0.1.1 stages links to both samples' actual REDUX BQR, jitter
parameters and microsatellite tables in a single experiment-local directory,
passed via native `-jitter_bqr_dir`. SAGE 5.0.2 otherwise defaults that directory
to the first normal BAM's parent, which cannot locate the tumour sidecars when
published BAMs are in separate sample directories. No file is recomputed and
neither BQR nor jitter modelling is disabled. The original Nextflow task avoids
this by staging both samples together; this probe now supplies equivalent inputs.

All native branches must remain available for final integration. An intersection
of the two callers, synthetic read-context vectors, or switching off SAGE_APPEND
or NEO does not meet that requirement. The November 11 reviewed recommendations
supersede those shortcuts in the earlier draft wiring contract.

## Evidence and execution

Runtime packages and job IDs are recorded in
`tests/results/deepsomatic-distributed/EXECUTION.md`.

Source anchors: Google r1.10 `make_examples_core.processing_regions` and
`make_examples_options.runtime_by_region`; Hartwig tag `sage-v5.0.2`,
`SimpleVariant.loadSimpleVariants`, `SageConfig` specific_var_file handling,
`Candidates.candidates`, `VariantReadContext.trinucleotide`,
`VariantReadContextBuilder.setMaxRefRepeat`, and
`CandidateSerialisation.toContext` (MH/read-context semantics). The local
source review directory is recorded in the findings review.
