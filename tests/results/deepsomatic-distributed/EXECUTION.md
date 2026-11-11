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
