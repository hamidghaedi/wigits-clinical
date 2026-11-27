# Optional neo portability verification

27 November 2026. Release native-intake/v1.1.5.

Offline regressions cover nonempty mixed-origin records, contig/position order,
empty VCFs, undeclared contigs, retained annotations/genotypes/phasing/duplicates,
explicit container Python isolation, completed CH status reuse, failure precedence
and preserved/mismatched export handling.

| Verification | Outcome |
|---|---|
| `python -m unittest discover -s tests -q` |193 tests,3 existing skips; exit0. Includes executable recovery-head tests with synthetic module/receipt bindings.|
| `python -m pytest tests/test_neo_interface.py tests/test_neo_export.py workflow/tests/test_chip_curation.py -q` |54 passed, exit0.|
| `python -m unittest discover -s workflow/sage_compat/tests -q` |3 passed, exit0.|
| Tracked shell/Slurm/env `bash -n`; `git diff --check` |Exit0.|
| Real poisoned user-site and injected PYTHONPATH probe |Head12516154: old invocation fails deliberately, fixed container NumPy/pandas import and pVACseq CLI pass in both cases; synthetic interleaved4-record VCF indexes with records preserved.|
| Actual failed VEP VCF probe |12516158 COMPLETED0:0,7s; actual4-record annotation sorted/indexed with complete record/header retention; same container isolation probes pass.|
| Prediction-only Nextflow graph |Child12516155 COMPLETED0:0; exactly one NEO_PREDICT task under `PREDICT_EXISTING_NEO`, zero CH/preparation tasks. Stub output is labelled synthetic.|
| Full real reference prediction/review |Head12516154 COMPLETED0:0,19m47; child12516157 COMPLETED0:0,17m59. Original native/CH/export products reused.|
| Previously failing case replay using transferred real annotation/context/export files |Generic recovery head12516161 COMPLETED0:0,14m40; child12516162 completed in13m56. Exactly one prediction task; all66 selected/4 proximal records retained;63 candidate rows/66 columns and native NEO comparison written. No BAMs/FASTQs transferred or native/CH/export task rerun.|
| Original small CH/export preservation |All16 recorded files retain identical contents and mtimes. No bulk BAM/FASTQ hashing.|

Published tested application source: `095f016` and `50b6332`; their complete
`workflow/` trees are identical. Subsequent release documentation does not
change that tested code. Resources/images/models and scientific thresholds are
unchanged. Exact commands, accounting and receipts remain outside Git at
`/global/scratch/hpc6049/optional_neo_fix_20261127/`.

## Output comparison and limits

The reference re-prediction retains69/69 candidate identities, all66 columns,
all labels and every non-affinity field.51 rows are byte-identical as parsed
records;18 differ only in `ic50_mt`, `ic50_wt` or derived `fold_change`.
Maximum absolute differences are0.001nM,0.003nM and0.00001213 respectively;
maximum relative difference is0.000453%. This is consistent with numerical
inference variation, whose precise source was not isolated. The exact full-row
comparison therefore exits1; it is retained alongside the detailed comparison,
not relabelled as bitwise agreement. No threshold or selection rule was tuned.

Both the earlier successful reference and new replays expose the pre-existing
incomplete pVACseq proximal-phasing interface (gap `OP-NEO-PHASING`). Its missing-
anchor warnings and empty proximal-context tables are not cleared by a zero
execution exit. No clinical/vaccine validity or universal portability is claimed.
The actual Fir retry is prepared for the owner, not submitted by these tests.

## v1.2.0 common-launcher extension

The owner then requested explicit entry points in existing reusable scripts.
Published tested source `7461da8` adds `--entry native|neo-predict|neo-review`
to the existing launcher and generalises its internal optional Slurm head.
No case-specific operator script is added. The prediction/export/model/review
implementation files are byte-identical to v1.1.5.

- Primary suite:200 tests,3 existing skips; exit0. Native default behavior,
  optional entry dispatch, receipt separation, repeat-submission handling,
  dependencies/walltime overrides and required review inputs are exercised.
- Existing54 feature tests and3 adapter tests pass; shell/diff checks pass.
- The existing launcher submitted real review-only head12516192 (COMPLETED0:0,
  1m14) and child12516193 (COMPLETED0:0,24.9s). The trace contains only
  `REVIEW_EXISTING_NEO:NEO_REVIEW`.63/63 candidate rows match all66 columns
  exactly, with no repeated model inference/native/CH/export work.
- Exact command, submission receipt, trace and full-row comparison are retained
  in the same private evidence directory. Source/image/model/threshold bindings
  and the separately recorded phasing limitation remain unchanged.
