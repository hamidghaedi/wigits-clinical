# KQ0205 pipeline performance and output comparison

Snapshot: **11 November 2026, approximately09:09 EST**. This is a maintained
comparison of analytical outputs, not a new report or clinical clearance.
Inputs, jobs and outputs were inspected read-only. No new analytical jobs were
submitted. Refresh pending rows when their source jobs finish.

## What we can conclude now

The two current Hartwig routes give the same primary PURPLE purity/ploidy
(24%/2.90) and the same two QC warnings. The native-reference route therefore
has **not resolved these CN warnings**. This is a completed PURPLE result even
though its containing workflow is still waiting for VirusBreakend.

The short-variant sets delivered to/through PURPLE differ substantially in
size and exact allele overlap. That is a reason to investigate their membership,
not evidence that one caller is more accurate. The old-to-new comparisons change
caller, SV caller, preprocessing, PURPLE version and annotation/filtering together.
The two current SAGE routes also differ in reference and BAM preparation.
The cleanest planned short-caller comparison—SAGE versus DeepSomatic on the
**same native REDUX BAMs**—is still pending DeepSomatic completion.

No threshold for acceptable disagreement has been invented. Counts, effect
sizes, QC and evidence of changed findings should inform the eventual decision.

## Routes and exact evidence roots

Paths below are shorthand used in the tables; all are full-size KQ0205 inputs.
The earlier20Mb/chr17+chr19 pilot is excluded from whole-genome comparisons.

| Key | Route and versions | Evidence root |
|---|---|---|
| A | Existing DNA workflow; Google DeepSomatic1.10.0 FFPE_WGS; PURPLE4.3; GRCh38.d1.vd1 | `/global/scratch/hpc6049/mohcc_runs/gi_frontenac_20261103/KQ0205/` |
| B | Existing BAMs, mate preparation then native REDUX2.0.5; SAGE5.0.2/PAVE1.9/ESVEE2.0.1/PURPLE4.4; custom GRCh38.d1.vd1 | `/global/scratch/hpc6049/wigits_runs/KQ0205_full_bam_custom_20261109_r1/` |
| C | Original full DNA FASTQs, four lanes per sample, native HMF masked reference; BWA-MEM2 2.3 then REDUX; same native callers as B | `/global/scratch/hpc6049/wigits_runs/KQ0205_full_native_20261109_onco3.0.0/` |
| D | Google DeepSomatic1.10.0 FFPE_WGS on the completed REDUX BAMs from C; same HMF reference | `/global/scratch/hpc6049/wigits_runs/KQ0205_full_deepsomatic_native_20261110/` |

B/C use nf-core/oncoanalyser3.0.0, revision
`7c74c87a43749952b38c9a18915947570f0595a0`, Nextflow26.04.4. B's consolidated
`results/pipeline_info/nf_core_oncoanalyser_software_mqc_versions.yml` records
the versions above. C's consolidated version file is still incomplete; its
completed VCF headers independently confirm SAGE5.0.2/PAVE1.9/PURPLE4.4.
A's `launch/release.txt` records checkout commit
`8e80d456a7ea5de9073567249c5b3ac5ceeea6d6`; its Nextflow log identifies the
entry revision `da86c0e5d0` at `dnaseq_nf-gi-closure-20261019/main.nf`.
D uses frozen runner0.2.0 (`b7510cd`) with an owner-approved inference-command
batch256 correction; maintained future launchers are commit `d233320`.

## Completion versus availability

| Output/work | A | B | C | D |
|---|---|---|---|---|
| Caller VCF | Completed DeepSomatic | Completed SAGE | Completed SAGE | Pending; CPU examples running |
| PURPLE primary CN/purity | Available, warnings | Available, warnings | Available, warnings | Not yet run with D calls |
| SV | Existing GRIDSS/GRIPSS route | ESVEE completed | ESVEE completed | No separate SV change intended |
| RNA | Outside this DNA release | Not supplied in B samplesheet | STAR/ISOFOX completed | Outside caller-only scope |
| ORANGE | Not this historical reporting route | PDF/JSON present | Pending workflow dependencies | No independent report requested |
| Whole execution | Historical completed analysis | Resumed workflow succeeded | Resumed workflow running | Running |

B's resumed head12449093 finished0:0 in6h15m51s:23 new tasks,12 cached,
zero failed tasks. This is **resume elapsed time**, not the complete cost from
original BAMs. C's resumed head12449091 was6h30m elapsed; VirusBreakend12449525
was active after passing the prior count-parser failure. D example job12447879
was14h50m elapsed; inference12447880 and postprocess12447881 remain dependent.
The inference command still explicitly contains `--batch_size 256`.

## Performance: measured task execution, not queue time

Nextflow `realtime`, `%cpu` and `peak_rss` are shown. CPU equivalents are
`%cpu / 100`; they are averages across a task, not instantaneous utilization.
Trace memory units are reproduced as supplied. Cached rows retain the original
task measurement and do not represent repeated work during resume.

| Task | Route | Job | Runtime | CPU allocation | Average busy CPU equivalents | Peak RSS | Resume status |
|---|---|---:|---:|---:|---:|---:|---|
| DeepSomatic, all3 phases | A |12414711|11h35m45s|32|17.32|60.7GB|Historical completed |
| REDUX normal | B |12446408|28m20s|64|17.62|31.2GB|Reused |
| REDUX tumour | B |12446409|1h11m30s|64|15.33|34.1GB|Reused |
| SAGE somatic | B |12446634|3h37m46s|64|60.27|96.8GB|Reused |
| ESVEE | B |12449096|4h04m43s|64|52.07|344.5GB|Completed after targeted JVM workaround |
| PURPLE | B |12449991|38s|—|6.74|13.2GB|New completed |
| REDUX normal | C |12444775|31m36s|64|16.95|29.4GB|Reused |
| REDUX tumour | C |12446060|2h09m14s|64|11.36|46.9GB|Reused |
| SAGE somatic | C |12446324|3h26m03s|64|60.21|182GB|Reused |
| ESVEE | C |12446323|12h49m14s|64|2.68|441.9GB|Reused |
| PURPLE | C |12449535|38s|—|6.65|12.8GB|New completed |
| DeepSomatic examples only | D |12447879|14h50m and running|64|Approximately63|Approximately104GiB|Running; not a final measurement |

C's four normal DNA alignment tasks each took50m08s–52m39s; four tumour tasks
2h25m32s–2h30m39s, with64 allocated CPUs and approximately42–44 busy CPU
equivalents. These execute in parallel; do not sum them into elapsed runtime.
RNA STAR took1h08m20s, sorting53m07s, duplicate marking3h10m18s, ISOFOX40m12s.
RNA scope exists only in C, so overall B/C walltime is not a matched comparison.
B also incurred initial BAM preparation/reference building and failed attempts;
those must be included before claiming that either route saves total time.

A's native DeepSomatic phase times were6h27m examples,4h57m inference and10m35s
postprocessing. D has already exceeded the earlier example phase. Its latest
periodic counters cover64 shards, approximately50.7million candidates and43.1million
examples, versus A's approximately33.7million/27.1million final periodic counters.
These are **not final somatic calls or percent genome completed**. Hardware,
alignment/reference and candidate workload all changed; no single cause of the
runtime difference is isolated. See the linked [runtime investigation](../deepsomatic-substitution/RUNTIME_REVIEW_20261111.md).

Performance source files:

- A: `launch/nextflow_trace.txt` and `work/22/70ba719541a0567e75b335e070cb85/KQ0205.deepsomatic_ffpe.run.log`.
- B: `results/pipeline_info/execution_trace_2026-11-11_00-02-08.txt`.
- C: `results/pipeline_info/execution_trace_2026-11-11_02-37-34.txt` (still growing).
- D: `results/examples.log`, stage JSON records and `sacct`/`sstat` for the listed jobs.

## Available CN and biomarker outputs

These are native analytical values, **not an endorsement for report inclusion**.

| Native PURPLE field | A | B | C |
|---|---:|---:|---:|
| Purity |21%|24%|24%|
| Purity fit bounds |18–25%|20–28%|20–28%|
| Ploidy |2.86|2.90|2.90|
| Ploidy fit bounds |2.70–3.00|2.68–4.90|2.70–3.05|
| Somatic penalty |0.0000|0.0302|0.0308|
| CN segments |402|1,364|1,384|
| Unsupported CN segments |313|948|986|
| Deleted genes |648|649|649|
| TMB/Mb |5.6289|2.2284|2.0549|
| Native MSI label |MSS|MSS|MSS|
| Microsatellite indels/Mb |0.0000|0.0644|0.0560|
| Native SV tumour mutational burden |1|111|90|
| QC |WARN_DELETED_GENES; WARN_HIGH_COPY_NUMBER_NOISE|Same|Same|

Native MSS output does not supersede the project's purity-dependent reporting
policy. The old zero somatic penalty is an observed field, not proof that its
short-variant input supplied the same context used by SAGE/PAVE in B/C.

Sources: A `results/cnv/KQ0205/purple/purple/`; B/C respective
`results/KQ0205_full_bam_custom/purple/` and `results/KQ0205_full_native/purple/`:
`KQ0205-DNA.purple.purity.tsv` and `KQ0205-DNA.purple.qc`.

## Small-variant overlap: completed PURPLE-facing outputs only

This first inexpensive comparison reads the small completed
`KQ0205-DNA.purple.somatic.vcf.gz` files in those directories. **It is not a
comparison of all raw DeepSomatic/SAGE calls**. A's file contains onlyPASS
records; B/C also retain PON/PONGnomad-tagged records. The filtering/annotation
paths differ, and raw caller FILTER labels must not be treated as equivalent.

| File-level observation | A | B | C |
|---|---:|---:|---:|
| PASS VCF records |16,094|6,671|6,134|
| Unique exact PASS alleles |16,111|6,671|6,134|
| PONGnomad records |0 in this derivative|589|523|
| PON;PONGnomad records |0 in this derivative|437|386|
| PON records |0 in this derivative|396|304|

Method: FILTER exactly`PASS`; splitALT on commas; compare unique
`(CHROM, POS, REF, ALT)` keys as written. No normalization, liftover or
representation-equivalence inference was performed. Multi-ALT splitting accounts
for A's record/allele difference. Indel representation and reference masking can
make exact overlap underestimate biological agreement.

| Comparison | Shared exact alleles | Only first | Only second | Jaccard |
|---|---:|---:|---:|---:|
| A versus B |1,260|14,851|5,411|5.85%|
| A versus C |1,181|14,930|4,953|5.61%|
| B versus C |4,939|1,732|1,195|62.79%|

B/C therefore share4,939 exactPASS alleles, with C having537 fewer netPASS
alleles (8.05% fewer than B). That net reduction hides1,732 B-only and1,195
C-only alleles. The very low A-versus-native overlap needs a review of caller
and filter behavior before any accuracy claim. D is needed to separate the
caller effect from the alignment/reference effect more cleanly.

## Structural and RNA outputs

PURPLE SV VCF **records**, not unique biological events:

| FILTER | A | B | C |
|---|---:|---:|---:|
| PASS |2|227|205|
| INFERRED |352|1,092|1,131|
| PON |2|9|13|

SV breakends can occupy multiple records; INFERRED records are not directly
observed ESVEE events. Do not interpret these counts as a validated hundredfold
gain in SV sensitivity. GRIDSS/GRIPSS-to-ESVEE and PURPLE version/input changes
are confounded. Match events by type, breakpoints, orientation and support in a
later bounded comparison; literal VCF ID/row equality is not adequate.

B has a completed ORANGE5.0.1 PDF/JSON and LINX output. C has LINX, CHORD,
signature and RNA results available, but ORANGE is not yet published.
C's ISOFOX summary reports QC=PASS,645,219,066 total fragments and419,174,197
duplicate fragments. This is an RNA tool summary, not a laboratory acceptance
decision. B's samplesheet contains tumour/normal DNA only; there is no matched
RNA arm in B. Consequently missing B RNA is out of scope, not an analysis failure.

## Refresh and decision checklist

- [ ] Record terminal status, job times and final software versions for C and D;
  preserve failed attempts separately from successful execution.
- [ ] Capture D raw VCF/filter counts after postprocessing finishes. Larger raw
  caller VCF scans/normalization should be done under Slurm, not as login-node
  full-genome work. Keep these outside the above PURPLE-derivative comparison.
- [ ] Compare SAGE and DeepSomatic on C's identical BAMs/reference using explicit
  allele normalization and separate raw-caller/PASS/annotation-filtered scopes;
  report SNV and indel results separately, VAF distributions, normal support and
  known reviewed findings. Do not change either caller to manufacture agreement.
- [ ] Once the approved adapter has run D calls through PAVE/PURPLE with the same
  AMBER/COBALT/ESVEE inputs, compare purity/ploidy fits, penalty, CN gene/segment
  changes and QC. D's caller alone cannot establish a PURPLE outcome.
- [ ] Compare B/C CN segments by genomic overlap and gene copy numbers, and SVs
  as events with support. Check whether differences involve drivers or known
  reviewed findings rather than only global counts.
- [ ] Review completed ORANGE outputs side by side for changed findings, not
  just PDF presence; keep RNA-only differences separate.
- [ ] Summarize active compute, queue time, prep/reference setup and reused work
  separately before deciding future resource sizing or preferred entry route.

Update this document in place with a new dated snapshot and exact evidence;
do not substitute absent output with zero or declare the comparison complete
while the same-BAM DeepSomatic arm is pending.
