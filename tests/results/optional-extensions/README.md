# Optional extension integration verification

26 November2026; integration branch `codex/optional-chip-neo` from published
v0.2.0 plus `chip/reporting-layer` at `b618a90`. Neither published main nor
issued execution packages have changed.

| Check | Outcome |
|---|---|
| `python -m unittest discover -s tests -q` |154 tests,3 existing skips; all executed tests pass.|
| `python -m unittest discover -s workflow/sage_compat/tests -q` |3 tests pass.|
| `python -m pytest tests/test_neo_interface.py tests/test_neo_export.py workflow/tests/test_chip_curation.py -q` |54 tests pass.|
| `bash -n` over workflow shell/Slurm sources |29 files pass.|
| Installed pVACtools CLI inspection on a compute node |Slurm12514452 COMPLETED, exit0:0,23s;1CPU/4GiB/10min allocation. Help for run and five downstream commands, plus installed argument-parser source, retained outside Git.|
| Nextflow branch/option wiring |Slurm12514464 COMPLETED0:0 in35s, source `f0f9f4d`; both, CH-only, neo-only and neither pass under explicit stub execution. This is graph/interface validation, not scientific validation.|
| Real blood caller/curation |Corrected source `2fb4cdf`, child12514477 COMPLETED0:0 in2m22s; source BAMs unchanged, no staged BAM links published.|
| Real neoantigen enrichment |Head12514475 / child12514479 in progress at this checkpoint, with the native export and local predictions running. End-to-end success is not yet claimed.|
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
