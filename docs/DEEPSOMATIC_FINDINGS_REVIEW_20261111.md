# Review of DeepSomatic performance and downstream integration proposals

11 November2026. Review of `7069cb5`; recommendations, not implementation or
approval of new scientific defaults. No running or queued job was changed and
no analytical job was submitted. The existing batch256 correction remains.

## Recommendation

Keep the native REDUX inputs and three caller phases. Prefer a locally measured
GPU type for future inference; test distributed example generation before
altering input reads. Continue the adapter work, but do not adopt the proposed
feature-disabling or synthetic read-context shortcuts. Full-feature replacement
is not demonstrated yet. Retain the native SAGE arm as the complete control
while resolving context for DeepSomatic-only alleles.

## Performance decisions

| Proposal | Recommendation | Evidence and limit |
|---|---|---|
| Separate CPU examples, GPU inference, CPU postprocessing | Keep | Avoids reserving a GPU throughout example generation. Already implemented. |
| Batch256 | Keep | Existing local memory/canary evidence; already explicit in both launchers. |
| Prefer A30 over untyped GPU request | Adopt for a future tested profile | Previous A30 throughput is better than the observed RTX6000 runs. Hardware comparisons are not controlled at identical batch/workload; compare queue plus execution time. A100 remains a benchmark candidate, not an assumed rate. |
| More CPU nodes for examples | First performance experiment | Use the existing global shard namespace across Slurm array jobs, each owning disjoint task IDs. Initially compare one64-core job with two32-core jobs at the same64 total shards; then assess scaling to more shards/nodes. Keep native callable scope, model and thresholds unchanged. One central set of examples/sidecars feeds native inference. Do not promise linear speedup. |
| Independently selected chromosome jobs | Defer behind ordinary shard distribution | Extra interval-boundary and output-gathering logic is unnecessary for the first experiment. If later used, preserve the native callable scope rather than asserting every FASTA contig is necessarily called, and check boundary behavior. |
| Two GPU inference jobs | Later bounded test | Both consume disjoint example shards; gather all main and small-model outputs exactly once before joint postprocessing. Stock r1.10 does not establish support for arbitrary comma-separated input specifications in these paths. Canonical shard staging or a tested gather is required. |
| Adapter trimming before alignment | Test only | Changes inputs to SAGE, ESVEE and other consumers, not just DeepSomatic. Measure caller accuracy/concordance as well as time. |
| Hard-clip all soft clips in a DeepSomatic-only BAM copy | Do not adopt as default | Removes genuine variant/indel context as well as possible adapter sequence. It is not equivalent to adapter-aware FASTQ trimming. A narrowly designed comparison could test it, but speed alone cannot justify it. |
| Disable local realignment or raise candidate/VAF thresholds | Do not adopt | Changes the scientific caller behavior; not a resource fix. |
| Memory/time reductions | Recalibrate after completion | Useful for scheduling, not evidence of faster computation. Preserve headroom across cases; current sample has not established a universal160GB/16h/2h policy. |
| Internal thread limits | Hypothesis only | Prior process inspection saw one active thread and288 sleeping threads per inspected worker. A large thread count alone does not demonstrate harmful oversubscription. |

Native `--runtime_by_region` produces a TSV of reads, candidates and phase
timings. Enable it in a future bounded benchmark rather than building a new
monitoring layer. The get-reads timing includes realignment in r1.10; it does
not by itself separate assembly cost from BAM retrieval. The strongest initial
benchmark holds BAM, reference, regions, model, batch and total shards constant.

### Corrections to the causal claims

- The reported84–86% versus8–14% soft clipping is from **three10kb windows**,
  not a whole-BAM rate. It is a useful lead, not proof of adapters or the full
  runtime mechanism. Source code confirms soft clips can trigger assembly-window
  selection, and that realignment is enabled by default. Profiling or a matched
  intervention is still needed to establish its contribution here.
- The caller already automatically enables `trim_reads_for_pileup` with the
  FFPE alt-aligned pileup. This does not replace adapter trimming and is not a
  missing speed switch to turn on again.
- Per-shard file I/O accounting is not physical GPFS traffic or a reliable
  completion fraction. Google partitions regions among task IDs, rather than
  deliberately assigning every genomic region to every shard. Indexed block
  rereading/overlap can incur overhead, but the claim that64 shards necessarily
  perform64 complete physical BAM reads is not established. Do not use the
  proposed77–85% progress estimate as a measured ETA.
- Pilot versus full run and old versus new input comparisons also change
  hardware, reference, regions or preprocessing. Neither proves a universal
  32-core optimum or a particular trimming speedup.

## Downstream integration decisions

The native tumour `sage_dir` input is a useful insertion boundary: retain
REDUX, AMBER, COBALT, ESVEE, SAGE_GERMLINE and native PAVE/PURPLE consumers.
Keep a separate experimental output directory/cache for the substituted arm;
do not overwrite the native control while testing the alternate somatic input.

| Data/function | Recommended handling |
|---|---|
| DeepSomatic calls, FILTER, QUAL and tumour counts | Preserve original VCF and source meanings. Make a separate adapter derivative; exclude reference-only records from somatic candidates explicitly. Do not relabel germline/failed calls as passing somatic calls. Let native PAVE apply its configured filtering; do not prematurely discard evidence merely because it might be filtered downstream. |
| Normal counts | The completed local pilot already has FORMAT `NDP`, `NAD`, `NAF`. Map these measured values into the normal sample column after checking full-run completeness and allele indexing. Preserve original fields. Do not automatically re-pileup all loci, or invent normal GT=0/0 from the somatic output. If a normal genotype is not supplied, keep it missing; absent counts are not zero. |
| TIER, repeats, trinucleotide context, microhomology | Reuse pinned native resources/implementations where semantics match. These are real annotations, not evidence of SAGE having called the allele. Missing TIER impairs ordinary PURPLE fitting, although hotspot-type exceptions exist in its source. |
| RNA support / SAGE_APPEND | Keep enabled. It reconstructs a measured variant read context from `RC_INFO`. Ordinary reference flanks are not a proven substitute. Intersecting with native SAGE provides context only for matching alleles and leaves DeepSomatic-only variants unresolved. |
| NEO | Keep enabled. Native `RC_CNT` categories are FULL, PARTIAL_CORE, CORE, REALIGNED, REF, TOTAL. `[AD_alt,0,0,0,AD_ref,DP]` invents those categories; matching the vector length or avoiding an exception is not validation. NEO1.3 sums the first four support categories and uses genotype DP for coverage. |
| Phasing | Preserve genuinely measured phase information with exact allele/representation mapping. Do not silently accept loss of LPS/MNV effects as a feature-preserving replacement. SAGE-only intersections cannot establish phase support for every new allele. |
| TINC | Preserve native contamination diagnostics as an independent result during validation. Do not copy SAGE recovered status to DeepSomatic calls or claim equivalent contamination recovery without an explicit tested policy. |
| LINX and other downstream outputs | Keep their upstream inputs/native rules, but compare results. The same ESVEE calls do not guarantee identical LINX output if PURPLE fits, CN or enriched SV inputs change. |

### What native SAGE can and cannot currently provide

SAGE5.0.2 `specific_var_file` defines regions around requested variants and
filters the candidates it discovers. `CandidateStage` still discovers candidates
from tumour reads; `Candidates.candidates` retains matching discovered alleles.
It is **not an established force-genotyping/annotation mode for arbitrary
DeepSomatic alleles**. Its append mode likewise expects existing SAGE read
context and does not discover missing tumour context.

The preferred development direction is DeepSomatic discovery plus genuine
native SAGE read-context/evidence collection, with the two caller's fields kept
distinct. First reuse already available SAGE context for exact matched alleles.
Then test whether the native restricted pass supplies complete context for the
remaining alleles; it may not. If it cannot, a minimal measured-context extension
or explicit consumer support is needed. Do not fill missing fields, disable
consumers, or restrict the final list to the intersection and call that complete.
Retaining a full native SAGE control for now adds cost but preserves the full
reference output while this integration gap is resolved.

## Next sequence of work, subject to implementation instruction

1. Let current full-genome execution finish and preserve its baseline outputs.
2. Benchmark existing shard distribution across nodes and suitable GPU placement;
   add native region timing to the benchmark. Assess sequence-preserving options
   before altering BAMs. Compare final allele/filter/count results, not just speed.
3. Map real DeepSomatic normal counts and native reference annotations into a
   separate adapter output; keep both input VCFs intact.
4. Resolve authentic context/phasing for DeepSomatic-only alleles. Demonstrate
   RNA SAGE_APPEND and NEO consumption without disabled branches or invented counts.
5. Run native PAVE/PURPLE and all consumers on the substituted arm and update
   the maintained KQ0205 comparison, including changed findings and failures.
   Only then declare the integration feature-complete and move to a single
   optional Nextflow caller branch. Native SAGE remains the control/default.

## Source checks

Pinned public source snapshots and retrieval manifest:
`/global/scratch/hpc6049/deepsomatic_wiring_review_20261111/`.
The local pilot VCF header was read directly; no whole-BAM scan was performed.

- Google r1.10 [make_examples_core.py](https://github.com/google/deepvariant/blob/r1.10/deepvariant/make_examples_core.py): processing_regions task partitioning; region_reads realignment; runtime reporting.
- Google r1.10 [make_examples_options.py](https://github.com/google/deepvariant/blob/r1.10/deepvariant/make_examples_options.py): realign_reads default, trim_reads_for_pileup and runtime_by_region.
- Google r1.10 [window_selector.py](https://github.com/google/deepvariant/blob/r1.10/deepvariant/realigner/window_selector.py): soft clips among assembly-window signals.
- Google r1.10 [call_variants.py](https://github.com/google/deepvariant/blob/r1.10/deepvariant/call_variants.py): prepare_dataset normalizes one sharded pattern; [postprocess_variants.py](https://github.com/google/deepvariant/blob/r1.10/deepvariant/postprocess_variants.py): get_cvo_paths and joint small/main-model input handling.
- SAGE5.0.2 [Candidates.java](https://github.com/hartwigmedical/hmftools/blob/sage-v5.0.2/sage/src/main/java/com/hartwig/hmftools/sage/candidate/Candidates.java): restrictedVariants filters discovered candidates; [CandidateSerialisation.java](https://github.com/hartwigmedical/hmftools/blob/sage-v5.0.2/sage/src/main/java/com/hartwig/hmftools/sage/vcf/CandidateSerialisation.java): toCandidate dereferences RC_INFO.
- NEO1.3 [SomaticVariant.java](https://github.com/hartwigmedical/hmftools/blob/neo-v1.3/neo/src/main/java/com/hartwig/hmftools/neo/score/SomaticVariant.java) and [VariantReadSupport.java](https://github.com/hartwigmedical/hmftools/blob/neo-v1.3/hmf-common/src/main/java/com/hartwig/hmftools/common/variant/VariantReadSupport.java): actual support categories and consumption.
- Native oncoanalyser3.0.0 local clone `7c74c87a43749952b38c9a18915947570f0595a0`: sage_calling input gating and modules/local/sage/somatic/main.nf TINC flags; completed B software manifest confirms NEO1.3, not an assumed NEO snapshot from another tool's tag.
