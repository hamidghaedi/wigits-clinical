# KQ0205 distributed examples and native annotation experiment

11 November 2026. Original native workflows and DeepSomatic jobs were not
modified, held, cancelled or restarted. New jobs use frozen source copies and
separate scratch output trees. Existing example preparation completed naturally
before this experiment; original GPU inference is running.

## Distributed examples

Runtime root:
`/global/scratch/hpc6049/wigits_runs/KQ0205_deepsomatic_distributed_20261111/`.

| Job | Purpose | State at inspection |
|---|---|---|
| 12447879 | Original 64-core baseline | COMPLETED, exit 0; 19h47m42s; 109137892 KiB peak batch RSS; frnt128 |
| 12452810 | Same 64 global shards, 2 nodes × 32 CPUs | RUNNING on frnt123 and frnt124 |
| 12452811 | Compare decompressed main/small-model payloads and example-info files | Pending afterok:12452810 |

`submission.json`, `comparison_submission.json`, `resources.json` and
`release.json` record exact commands, accepted allocations and frozen component
0.1.0 hashes. `results/plan.json` records the unchanged scientific command and
baseline. `results/workers/0.json` and `1.json` record nodes and disjoint task
ownership; per-shard logs show native task IDs with the denominator 64.
Runtime TSVs and comparison.json are pending completion. No speedup or output
equivalence is claimed yet. Original GPU/postprocessing jobs 12447880/12447881
and native-run jobs remain separate and untouched.

## Measured-count/native-annotation adapter

Runtime root:
`/global/scratch/hpc6049/wigits_runs/KQ0205_deepsomatic_adapter_20261111/`.

The original completed 20 Mb pilot contains 651526 records: 597564 RefCall,
38683 GERMLINE, 15096 NoCall and 183 PASS. All 651526 supplied NAD/NDP; 308
were multiallelic. These are raw caller classes, not all somatic mutations.
The explicit somatic derivative excludes only RefCall/GERMLINE, retains all
other original FILTERs and the complete original caller file, and does not
replace the independent native germline branch.

| Job | Purpose | State at last recorded inspection |
|---|---|---|
| 12452819 | Initial all-record annotation exercise, component 0.1.0 | Native counts, projection, tiers, repeats and TNC succeeded. Final exact-text preservation check stopped on Java VCF serialization. Outputs preserved under pilot/. |
| 12452826 | Explicit candidate derivative, component 0.1.1 | Stopped on equivalent missing-vector serialization. Frozen package and output preserved. |
| 12452827 | Native context probe dependent on 0.1.1 | Automatically cancelled because the prerequisite failed; no analysis ran. |
| 12452831 | Pilot annotation with source payload retained verbatim, component 0.1.2 | All annotation/preservation/compression steps succeeded; audit stopped on multiple native SAGE records for an exact allele. |
| 12452832, 12452833 | Dependent probe/full adapter for 0.1.2 | Automatically cancelled because the prerequisite failed; no analysis ran. |
| 12452839 | Pilot annotation/audit, component 0.1.3 | COMPLETED, exit 0; all eight steps passed; see pilot_v013/execution.json |
| 12452840 | Native SAGE restricted-context probe on the pilot PASS set | FAILED in 5 seconds before calling: default BQR directory was the normal BAM directory, missing the tumour sidecar there |
| 12452902 | Corrected probe, component 0.1.1 | COMPLETED, exit 0, 2m01s, 24047532 KiB peak batch RSS; context_probe_v014/ |
| 12452841 | Same adapter on the full DeepSomatic VCF | Pending afterok:12452839:12447881; does not change either prerequisite job |

The first check exposed equivalent QUAL formatting (0.3 versus 0.30) and
missing-vector serialization (.,. versus .). Component 0.1.2 takes only the
new INFO annotations from Java and retains the exact normalized caller QUAL,
FILTER and genotype payload. Original native output is retained for inspection;
no candidate or measurement is invented to satisfy a format check.

The native SAGE file also legitimately contains multiple records at some exact
alleles. Component 0.1.3 records every matching context in a JSON column and
marks these as ambiguous; it neither chooses one nor fails the whole annotation
job. Only a unique native match counts as immediately available context.

**Verified pilot result (0.1.3):** 15279 candidate source records became 15281
biallelic rows after native projection of two multiallelic records. All allele,
FILTER, QUAL and genotype payload values were preserved through annotation.
Native reference search found 2116 qualifying repeats and 13165 measured
no-qualifying-repeat rows. The audit found 1899 exact native SAGE allele matches,
1894 with unique available read context and five with multiple native records.
Among the 183 DeepSomatic PASS rows, 68 have exact unique native SAGE context;
115 have no exact match. All 183 remain in the annotation product. These are
representation-exact comparisons, not claims that all unmatched alleles are
biologically absent from SAGE. The restricted native probe is pending; no RNA,
NEO or full downstream equivalence is claimed.

Native region and repeat helpers were ported from DNA adapter commit 7093033,
using SAGE 5.0.2 and the current oncoanalyser 3.0 resource bundle. Source DNA
checkout is unchanged. Actual normal counts come from DeepSomatic, not a new
BAM scan. Native read-context availability is audited separately; no RC_INFO,
RC_CNT, phase sets or contamination-recovered flags are fabricated or inferred.
No downstream consumer is disabled and no experimental sage_dir is published.

`submission_v0.1.3.json` records current exact commands. The full annotation allocation
is four CPUs, 48 GiB and eight hours, based on 15279 pilot candidate records
scaled by reference span (about 2.46 million estimated candidates); this is a
conservative proxy, not a measured full-case rate. `full_resources.json`
records that basis. The bounded SAGE probe has eight CPUs, 32 GiB and two hours
for at most 1000 supplied PASS alleles; restricted TINC/phasing output must not
be treated as a sample-wide result.

## Local verification

- `python3 -m unittest discover -s tests`: 48 tests, OK, 3 skipped (existing optional integration checks).
- `python3 -m py_compile workflow/benchmark_deepsomatic.py workflow/deepsomatic_adapter.py workflow/run_deepsomatic_adapter.py workflow/probe_sage_context.py`: exit 0.
- `bash -n workflow/benchmark_deepsomatic.sbatch workflow/deepsomatic_adapter.sbatch workflow/probe_sage_context.sbatch`: exit 0.
- `git diff --check`: exit 0 after retaining clean line endings for the added draft-contract notice.

Follow-up: the BQR staging correction adds two regression tests. Complete suite:
50 tests, OK, 3 existing skips. Original BAMs, sidecars and running jobs unchanged.
`probe_sidecar_submission.json` and `context_package_v0.1.1/release.json` record
the exact new command and source fingerprints. Missing tumour BQR was a staging
failure, not evidence against support for the 115 unmatched PASS alleles.

**Completed native context probe:** all 183 DeepSomatic PASS alleles were
requested. SAGE returned genuine context for exactly the same 68 alleles as
the original native control: zero additional supported alleles, zero lost
previous matches, 115 still without matching SAGE context. This was verified
by allele-key sets, not matching totals alone. Evidence:
`context_probe_v014/probe_audit.tsv`, its JSON summary, and
`comparison_to_native_control.json`. Both samples' actual BQR tables loaded
successfully (1536 records each), and native SAGE/TINC completed. Restricted
TINC/phasing remains unsuitable as a full-sample result. These findings do not
establish absent BAM support or false-positive status for the unmatched calls.
The restricted discovery route cannot supply feature-complete replacement
context; no variant was removed or synthetic context emitted.

At the follow-up check, distributed job 12452810 was running at approximately
23 minutes, with both nodes progressing. Full inference 12447880 was running
at about 2h32m. All 64 completed example summaries total 63898627 main-model
examples; the inference log had reached 7206656 (11.28%) at about 0.126 seconds
per 100 examples. Linear extrapolation suggests about 20 hours more inference
at that rate, excluding postprocessing, and is not a completion guarantee.

Remaining: benchmark completion/timing/equality; pilot/full annotation results;
native context coverage of DeepSomatic-only alleles; genuine RNA append, NEO,
phasing and TINC handling before complete native downstream substitution.


## Native evidence export and consumer increment (11 November; 12 November UTC)

Maintained implementation and limits: `docs/DEEPSOMATIC_EVIDENCE_ADAPTER.md`.
Runtime root (all paths below relative to it):
`/global/scratch/hpc6049/wigits_runs/KQ0205_deepsomatic_adapter_20261111/evidence_export_v0.1.0/`.
Original caller inputs, native control and executing full-case jobs were untouched.

| Job | Result | Evidence |
|---|---|---|
| 12454341 | COMPLETED 0, 7s; isolated SAGE5.0.2 evidence overlay built | `build/jar-sha256.txt`, `submission.json` |
| 12454342 | COMPLETED 0, 11s; 180/183 PASS alleles matched | `probe/execution.json`, `probe/probe_audit.tsv` |
| 12454391 | COMPLETED 0, 2s; reference normalization changes no PASS allele | `normalize_submission.json`, `reference_normalized.vcf` |
| 12454441 | COMPLETED 0, 10s; native ordered candidate export | `probe_r2/execution.json`, `probe_r2/probe.sage.somatic.candidates.tsv.gz` |
| 12454446 | COMPLETED 0, 2m15s; 180 context roundtrips, 180 native RNA append, 180 native NEO reader, 180 PAVE records | `consumers/execution.json`, `consumer_submission.json` |
| 12454454 | COMPLETED 0, 6s; initial observed-context experiment, **not accepted as final evidence** | `seed_measurements/`; the chr19 haplotype-only context does not encode the insertion in its CIGAR |
| 12454470 | COMPLETED 0, 6s; strict exact-CIGAR seed and original-BAM tumour/normal measurement | `seed_measurements_exact/`, `seed_package_exact/` |
| 12454471 | COMPLETED 0, 2m15s; all four consumer checks pass for 181 alleles | `consumers_r2/execution.json`, `consumer_r2_submission.json` |

The overlay JAR SHA256 is
`3a639126a511e96793fb97856d76b9c62c6dc7295af5ae89e01ce3df802b814d`;
original installed JAR SHA256 is
`c8945ba7d30e77f35564d98849c01474f73b1bb7af201d9fb17d93df3ca79f10`.
The patch changes only retention at two native filtering gates, requires explicit
requested alleles, and forbids restricted-scope TINC. All native counts, quality
metrics and filters remain measured; DeepSomatic decisions are not replaced.

The 180 discovery contexts comprise 176 unique and four multiple-context alleles.
The four are resolved using the native ordered primary candidate (full/core
matches then event count). All eight raw context records remain in the evidence
VCF; counts are not summed. The complete master preserves all 15281 normalized
candidate rows and all183 PASS calls with unchanged original genotype payload,
FILTER and QUAL. `all_with_evidence.vcf.json` records exact preservation.

The additional exact-CIGAR allele is chr17:2539597 G>GT. Native append measures
7 full-support tumour fragments /89 depth and1 /81 in normal; original caller
counts remain alongside these values. These are measurements, not a new
classification or an assertion of a clean normal. The failed haplotype-only
chr19 seed is excluded from the accepted supplement and is covered by a
regression test requiring the exact event in the serialized context CIGAR.

Still unresolved at this checkpoint: chr17:3383908
AACAAATATTTATTGAATAAATATTGTAT>A and chr19:8679269 C>CATT. Both remain in the
complete caller file. Diagnostic jobs12454450/12454459 used native Google
realignment but changed partition boundaries. The first deletion probe did not
reproduce its candidate. Job12454480 preserves the original1000-base partition
grid and reproduces the exact deletion with91 native supporting reads, agreeing
with the original caller. This is a diagnostic-scope correction, not a change to
any original call. `partition_aligned_realigner/candidates.json` records the
candidate, its actual supporting read names and sample labels privately.

Small inspection jobs12454458 and12454469 failed before measurement (Google
BAM diagnostic pairing flags, and importing zipped packaged Python modules,
respectively). Corrected diagnostics12454467/12454468 completed0. No production
inputs or resource installations changed. Native diagnostic BAMs can overwrite
normal with tumour at a shared region path; they are never used for final
fragment counting. A separate observer of native post-normalization reads is
being tested to avoid that ambiguity and capture actual later CIGAR changes.

This remains an incomplete full-feature substitution. NEO testing is its native
count reader, not prediction; full PURPLE/LINX/ORANGE, phasing and TINC behavior
remain to be demonstrated. Full DeepSomatic inference12447880 continues,
postprocessing12447881 remains dependent. No claim of complete readiness or
whole-genome compatibility is made from the181-allele consumer exercise.


### Final context increment result

All115 initially unmatched PASS contexts are now represented. Slurm12454494
completed0 in7s, measuring all three observed indel contexts with native
SAGE_APPEND against the original tumour/normal REDUX BAMs. Native context
selection for the two realignment-derived contexts used unmodified
`AltContext.addReadContext/selectCandidates` (job12454489, completed0 in7s).
`native_primary_inspection/sample_identity.json` verifies their observed read
identities against Google's native candidate support and tumour sample labels.
No artificial read, changed quality or invented support category was used.

Slurm12454495 completed0 in2m16s for all183 context/RNA/NEO-reader/PAVE records.
The subsequent review found five differences between the interim reference-only
repeat annotations and actual native read-context annotations: three different
repeat counts and two native absences. Merge component0.1.2 now retains both in
the master and projects the actual native context annotations into the consumer
view. Native absence is preserved explicitly, not guessed from missing data.
The original caller fields and complete candidate derivative remain unchanged.

**Final Slurm12454502: COMPLETED, exit0, 2m16s**, component0.1.2:

- 183 native context roundtrips;
- 183 native RNA append records;
- 183 accepted by the installed NEO1.3 count reader;
- 183 native PAVE records;
- Exact input/RNA/PAVE allele multisets equal; all68 original native matches
  plus all115 formerly unmatched alleles retained;
- All15281 normalized candidate rows preserved in the master, including every
  original PASS/NoCall, QUAL and genotype value.

Authoritative receipts: `consumers_v012/execution.json`,
`consumers_v012/final_allele_audit.json`,
`consumers_v012/all_with_evidence.vcf.json` and
`consumer_v012_submission.json`. The final allele-multiset audit was performed
read-only after completion; the maintained consumer test also checks that
multiset in future runs. Native PAVE may add its own filters; it does not remove
or rewrite the original master calls.

Three PASS calls have zero native SAGE fragment support despite positive
DeepSomatic AD (one deletion, two SNVs). `zero_native_support.json` records all
three with exact original/native counts and native filters. They remain present.
The deletion has Google91 ALT reads /141 depth, versus native0 fragments /30
context depth. These different methods/units must not be forced into agreement,
interpreted as clean absence, or substituted for one another in the master.
NEO's consumer projection uses actual native fragment depth with its fragment
numerator; original caller depth remains available. No equivalent sensitivity
or full downstream scientific validation is claimed.

The optional native post-normalization observer completed as12454496 (0,29s):
14 separate sample/region BAMs; native286 candidates,239 examples,83 small-model
examples, matching the unobserved partition-aligned diagnostic totals. Its first
two setup attempts12454485/12454490 stopped on ZIP symlink extraction and
protobuf-path precedence, respectively. The final observer uses the image's
original packaged runtime and source computations; no environment installation
or source BAM modification occurred. This observer was not needed to invent
or replace the accepted seed measurements.

Final local verification:63 unit tests, OK,3 existing skips; Python compilation,
shell syntax and `git diff --check` exit0. Full-genome context scaling,
phasing/MNV and TINC semantics, full NEO prediction and substituted
PURPLE/LINX/ORANGE comparison remain open. `full_downstream_ready` is deliberately
false. This closes the pilot's missing-context problem, not the whole production
caller-substitution validation.
