# Optional extension integration verification

26 November2026; integration branch `codex/optional-chip-neo` from published
v0.2.0 plus `chip/reporting-layer` at `b618a90`. Neither published main nor
issued execution packages have changed.

| Check | Outcome |
|---|---|
| `python -m unittest discover -s tests -q` |157 tests,3 existing skips; all executed tests pass.|
| `python -m unittest discover -s workflow/sage_compat/tests -q` |3 tests pass.|
| `python -m pytest tests/test_neo_interface.py tests/test_neo_export.py workflow/tests/test_chip_curation.py -q` |54 tests pass.|
| `bash -n` over all tracked shell/Slurm sources |35 files pass.|
| Installed pVACtools CLI inspection on a compute node |Slurm12514452 COMPLETED, exit0:0,23s;1CPU/4GiB/10min allocation. Help for run and five downstream commands, plus installed argument-parser source, retained outside Git.|
| Nextflow branch/option wiring |Slurm12514464 COMPLETED0:0 in35s, source `f0f9f4d`; both, CH-only, neo-only and neither pass under explicit stub execution. This is graph/interface validation, not scientific validation.|
| Real blood caller/curation |Corrected source `2fb4cdf`, child12514477 COMPLETED0:0 in2m22s; source BAMs unchanged, no staged BAM links published.|
| Real neoantigen enrichment |Head12514475 COMPLETED0:0 in29m06; child12514479 COMPLETED0:0 in27m30. Original native outputs were reused.|
| Review reuse / output fidelity |12514505 COMPLETED0:0 in1m38; all69 standalone variant/allele/peptide rows match in every column, without rerunning models.|
| Fir deployment / real-case extension comparison |Pending explicit site profile, assets and runtime validation; not claimed.|

Private technical evidence:
`/global/scratch/hpc6049/wigits_optional_extensions_20261126/`.
The login node cannot execute this container (`Failed to create user namespace`);
the compute-node inspection succeeded without an environment or package change.

New offline checks exercise independent options; blood/RNA eligibility from lab
metadata; absent RNA support; mismatched sample bindings; unavailable resources;
input-dependent allowances; failed/unpublished task receipts; preservation of
native completion; and command bindings to the completed caller/annotation code.

An inherited export test required some populated expression whenever the trimmed
fixture contained any ISOFOX genes, although none matched its selected variants.
It now checks every variant against the actual source join. Separate synthetic
matching records exercise both zero and positive TPM in the tumour genotype.
The exporter and source fixture are unchanged; no expression is fabricated for
real cases and no assertion is weakened to accept lost matching values.

The merged feature's vendor whitelist files remain byte-preserved, including
their upstream line endings. The source-history merge contained inherited CRLF
warnings; the subsequent integration increment uses LF and passes its own
`git diff --check`. No scientific thresholds were tuned for these tests.

## Real publication defect found and corrected

The first real canary12514465 exposed `publishDir(copy)` following the CH
script's staging links to original REDUX inputs. CH itself completed (child
12514466,2m26s) and export completed (12514467,4s), but copying those inputs
would unnecessarily duplicate large BAMs. Only that new canary head and its
prediction child12514468 were cancelled; unrelated workflows were untouched.
The partial53,635,212,500byte replica was removed from that canary's published
`blood_ch/stage/` after cancellation. Originals remain at their source paths.
The private cleanup receipt lists each removed replica; no source file was
deleted or rehashed.

Fix `2fb4cdf` removes only the task-created staging symlinks after successful
CH execution and before emitting its output directory. Unexpected regular
files in that input staging directory are preserved and stop publication.
Two regressions exercise both paths. The corrected real CH task publishes
only analysis products/provenance; no `stage/` or BAM copy is present.
Prediction reruns in a separate test directory; no native alignment or calling
is repeated. Both the original and corrected source checkouts are retained.

## Preserve the recorded aggregation, not accidental default drift

The successful full test reproduced all69 candidate identities/peptides,
measurements and final review labels.11 auxiliary pVACseq labels differed.
The old standalone aggregation metrics record RNA VAF0.25, depth10 and
allele-expression cutoff2.5, whereas the automatic aggregate inherits the
discovery call's0.0/5/0.0 settings. No prediction or measured evidence differed.

`674fb5d` restores explicit review aggregation, preserves the automatic aggregate
as separate evidence and adds a review-only Nextflow entry consuming completed
predictions. An initial reuse test12514499 stopped before writing because its
input folder name collided with its output folder. The source-preservation
guard caught it; `84cd94b` gives the staged inputs distinct names and adds a
regression. Test12514502 then completed0:0 in1m44, revealing that pVACtools
converts `--expn-val` to a2.5-fold allele-expression cutoff. Passing the derived
2.5 as the argument produced6.25 and five auxiliary label differences.

`2818716` uses the original CLI default1.0, yielding the recorded2.5 cutoff.
The final keyed comparison has69/69 identical candidates and zero changed
columns. This restores the documented standalone invocation; parameters were
not fitted to classifications. All models, annotation inputs and source
outputs are preserved, and prediction was not repeated for either review test.

The GitHub connector refused draft PR creation with HTTP403 (`Resource not
accessible by integration`). The feature branch is pushed; main is unchanged.
The previous direct-main exception applied only to v0.2.0 consolidation.
