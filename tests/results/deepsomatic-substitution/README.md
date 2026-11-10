# Native-reference DeepSomatic integration test

10 November2026. Source component0.1.0, commit `9d6ca38`, on
`hartwig_adoptation`. This is an isolated bounded caller test; neither live
native workflow nor the default caller is modified.

## Completed preparation

- Bound the completed tumour/normal REDUX BAMs produced by full native FASTQ
  realignment, their existing indexes and the exact HMF masked FASTA.
- Reused the installed official Google DeepSomatic1.10.0 GPU image and selected
  FFPE_WGS from the existing FFPE source assignment.
- Published source and an immutable experimental package; no BAM recovery,
  data transfer, environment installation or reference/caller substitution.
- `python3 -m unittest discover -s tests`:27 tests, exit0;3 existing real-tool
  tests skipped on login. New tests cover explicit model/sample/reference/region
  binding, CPU/partition limits, CUDA mask forwarding and raw/zero-call summaries.
- `bash -n workflow/test_deepsomatic.sbatch`, exit0;
  `git diff --check`, exit0. Scheduler test-only accepted16 CPUs/128GiB/one GPU/4h.

## Historical submission snapshot

Slurm **12446769**, running onfrnt153 when recorded. First scope is20Mb,
chr17:1–10000000 and chr19:1–10000000, read through the full BAM indexes.
The test uses16 CPUs,128GiB and one GPU; these are bounded-test resources,
not full-WGS preparation sizing. No GPU inference completion is claimed yet.

```text
/global/scratch/hpc6049/wigits_runs/deepsomatic_seam_20261110/
  case.json
  package_v0.1.0/{release.json,test_deepsomatic.py,test_deepsomatic.sbatch}
  caller.12446769.log
  caller_test_v0.1.0/{execution.json,runner.log,logs/,intermediate/}
  source/refs.json
```

`execution.json` contains the exact issued command, source-run binding and
Slurm ID. The raw VCF/sample/FILTER/FORMAT/INFO/version summary is added only
after a successful caller exit. The job is still generating examples; its CUDA
messages do not demonstrate successful GPU inference. Inference-device visibility
remains to be checked. Inspect device logs rather than treating allocation as proof.

## Completed pilot and full-run correction

Job12446769 completed0:0 in53m37s. GPU:0 inference is confirmed by TensorFlow
MirroredStrategy and successful predictions. Approximate stage times: examples
38m13s, inference12m13s, postprocessing1m53s. The indexed VCF has183 PASS calls
(136 SNVs,47 indel/complex alleles),38,683 GERMLINE,597,564 RefCall and15,096 NoCall
records. All records are within the requested20Mb. NDP/NAD/NAF normal evidence
is present despite the single tumour sample column.

The owner clarified that the requested execution was **whole genome**. The
regional test does not fulfil that request. Component0.2.0 now runs Google's
generated native commands as three dependent full-genome jobs, with no region
argument. Inputs are the completed full FASTQ-run REDUX BAMs, not the older
custom-reference BAM route. Jobs12447879/12447880/12447881 were submitted for
examples/inference/postprocessing. Source and runtime evidence are described in
the [substitution plan](../../../docs/DEEPSOMATIC_SUBSTITUTION.md).

Validation:31 unit tests passed (3 existing integration skips); bash syntax and
git diff checks passed. All three allocations passed scheduler test-only checks.
Full caller completion and downstream compatibility remain pending.

## Still pending

Whole-genome caller completion; paired normal/context
adaptation using current native resources; PAVE1.9/PURPLE4.4 comparison with
unchanged SV/CN inputs; exact SAGE side-output consumer resolution; eventual
optional workflow integration. Prior PURPLE4.3 experiments establish useful
groundwork, not completion of these3.0.0 tests. No equivalence, clinical
clearance, new report or production caller switch is implied.

See [the substitution plan](../../../docs/DEEPSOMATIC_SUBSTITUTION.md) for
pinned sources, supported injection boundary and explicit missing-data rules.
