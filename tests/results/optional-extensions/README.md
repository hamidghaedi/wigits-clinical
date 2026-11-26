# Optional extension integration verification

26 November2026; integration branch `codex/optional-chip-neo` from published
v0.2.0 plus `chip/reporting-layer` at `b618a90`. Neither published main nor
issued execution packages have changed.

| Check | Outcome |
|---|---|
| `python -m unittest discover -s tests -q` |152 tests pass,3 existing skips.|
| `python -m unittest discover -s workflow/sage_compat/tests -q` |3 tests pass.|
| `python -m pytest tests/test_neo_interface.py tests/test_neo_export.py workflow/tests/test_chip_curation.py -q` |54 tests pass.|
| `bash -n` over workflow shell/Slurm sources |29 files pass.|
| Installed pVACtools CLI inspection on a compute node |Slurm12514452 COMPLETED, exit0:0,23s;1CPU/4GiB/10min allocation. Help for run and five downstream commands, plus installed argument-parser source, retained outside Git.|
| End-to-end optional Nextflow execution |Pending source publication and HPC canary; not claimed.|
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
