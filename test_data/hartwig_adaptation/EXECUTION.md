# Preparation execution, 9 November 2026

Source commit: `9d883c669329939c16da2555e807b33c39578dbe` on
`hartwig_adoptation`. This is test-input preparation only.

| Library | Slurm job | CPUs | Memory | Time limit | Initial state |
|---|---:|---:|---:|---:|---|
| KQ0205 tumour DNA | 12440943 | 8 | 24 GiB | 12 h | RUNNING, frnt141 |
| KQ0205 blood-normal DNA | 12440944 | 8 | 24 GiB | 6 h | RUNNING, frnt142 |
| KQ0205 tumour RNA | 12440945 | 8 | 24 GiB | 6 h | RUNNING, frnt142 |

Submitted 17:15:53–17:15:56 EST. All three began running immediately. This
initial checkpoint is not evidence of completion or ready FASTQs.

Committed scripts were staged unchanged under:
```text
/global/scratch/hpc6049/hartwig_test_fastqs_20261109/KQ0205/package/
```

Exact submitted commands and job IDs:
```text
/global/scratch/hpc6049/hartwig_test_fastqs_20261109/KQ0205/preparation_jobs.json
```

Outputs and logs are beneath the same KQ0205 runtime root. Each library writes
`completion.json` only after gzip/FASTQ parsing, paired-name/count matching and
accounting for all selected primary BAM records succeed. Check all three
completion records and terminal Slurm accounting before calling the complete
test dataset ready. Preserve separately recorded singletons/other reads.

Validation before submission:

- `module load StdEnv/2023 samtools/1.22.1`.
- `python3 test_data/hartwig_adaptation/test_prepare_fastqs.py`: exit 0,
  one synthetic end-to-end regression passed. Eight primary records accounted
  for as three pairs, one singleton and one other read; OQ recovery and
  cross-chromosome/unmapped mates checked.
- `bash -n test_data/hartwig_adaptation/prepare_fastqs.sbatch`: exit 0.
- `git diff --cached --check`: exit 0.
- Slurm `--test-only` accepted both 12-hour and 6-hour resource requests.

Only small header/index inspections and the synthetic regression ran on the
login node. All real BAM extraction, conversion and FASTQ verification run in
the three Slurm allocations. No source BAMs, existing analyses or production
workflow settings were modified.

## Samplesheet and RNA completion checkpoint

The completed three-row native samplesheet is
`samplesheet.KQ0205.chr17_chr19.csv` beside this document. It passed the official
Oncoanalyser 3.0.0 JSON schema and the unmodified native `Utils.parseInput`
parser at `7c74c87a43749952b38c9a18915947570f0595a0`, using Nextflow 26.04.4.
The parser-only workflow resolved one group and all three intended sample/role
bindings; it launched no analytical processes. Validation artifacts are under
the runtime root's `samplesheet_validation/` directory; its `.nextflow.log`
is at the runtime root. This establishes metadata validity, not data readiness.

RNA job **12440945 COMPLETED, exit 0:0, elapsed 00:08:02**, with maximum RSS
908580 KiB. Its completion record reports **41,399,966 paired reads**, zero
singletons and zero other reads; all selected primary records were accounted
for. The R1/R2 compressed sizes are 3,022,368,548 and 3,121,712,533 bytes.
Tumour and normal DNA jobs were still RUNNING at this checkpoint. Do not launch
the complete DNA/RNA test until their completion records and successful Slurm
accounting are also available.

## Full FASTQ preparation completed

All three preparation jobs completed with exit **0:0**, and all three
per-library completion records passed paired-name/count and selected-record
accounting. Every library has zero singleton/other reads.

| Library | Paired reads | Elapsed | Maximum RSS |
|---|---:|---|---:|
| Tumour DNA | 102,114,544 | 00:24:05 | 2258736 KiB |
| Normal DNA | 45,474,999 | 00:16:40 | 1011748 KiB |
| Tumour RNA | 41,399,966 | 00:08:02 | 908580 KiB |

The data are ready for the chr17/chr19 execution test, not for whole-genome
scientific conclusions. Aggregate input receipt:
```text
/global/scratch/hpc6049/hartwig_test_fastqs_20261109/KQ0205/test_dataset_completion.json
```

## Metadata consumer verification

The earlier schema/parser pass was a structural check. Consumer-specific
findings and native limitations are now authoritative in
[`METADATA_CONSUMERS.md`](METADATA_CONSUMERS.md).

Slurm **12441043 COMPLETED, exit 0:0, elapsed 00:00:04, maximum RSS 132900 KiB**.
Probe source: `9ed2dab7aa1a9eb4783aebc9729552cc9b849b58`.
The installed NEO JAR was tested on one transcript from the staged public
reference. Native output:

```json
{"Carcinoid_match":false,"Small intestine_match":true,"Small_intestine_match":false,"unmatched_cancer_tpm":-1.0,"pan_cancer_tpm_available":true,"unmatched_equals_omitted":true}
```

The installed ORANGE JAR's `formConfiguredPrimaryTumorDoid` bytecode was
inspected: the DOID branch iterates an empty set; the display-text branch builds
a label with an empty DOID. This corroborates the tagged source review; it is
not a full ORANGE report execution.

Executed JAR SHA-256:

- NEO: `cd7e313d3cf6db24ee510887f5231f1464da628f5c2967ff79b1d4d3c8d9ec6c`.
- ORANGE: `2f64189b8680c62bec6e7c72cecc3cda51a9eb56950808f349cee7c6de666987`.

Evidence is in the runtime root's `terminology_review/native_binary_check/`;
exact launch and source commit are in `terminology_review/job_r2.json`.
Initial probe12441042 failed before Java execution because Slurm temporary
storage was not mounted in the container. The harness added an explicit bind
and the new job above passed. No native tool, reference or interpretation rule
was modified to obtain agreement.
