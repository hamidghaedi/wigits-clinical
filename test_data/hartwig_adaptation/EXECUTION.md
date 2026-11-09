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
