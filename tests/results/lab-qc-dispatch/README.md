# Early run-level laboratory QC dispatch

Owner request November24,2026. Component lab-qc/v0.2.0 adds a metadata-only
observer to the existing BAM-ready LAB_RUN_QC entry. Metric computation and
current v0.1.2 jobs remain unchanged. No analysis-head/caller dependency.

## Verification

All commands use the existing Python interpreter, with PYTHONPATH unset and
PYTHONNOUSERSITE=1. No environment or package change.

| Command/check | Outcome |
|---|---|
| `python -m unittest discover -s tests` | exit0;118 tests,3 existing skips |
| `python -m unittest discover -s workflow/sage_compat/tests` | exit0;3 tests |
| `bash -n workflow/run_laboratory_qc.sbatch` | exit0 |
| `git diff --check` | exit0 |
| Real source metadata-only manifest construction |2 sequencing runs;4+10 libraries;32+80 lane pairs; exact producing reference; task-local versions; no final combined version file read |

The12 dispatcher tests cover submission while SAGE is RUNNING, waiting for the
whole laboratory run, incomplete published copies, failed BAM tasks, baseline
exclusion, existing-job adoption, ambiguous submission, queue-pointer enrollment,
missing lanes, single delivery, failed-head visibility and use of the existing
stager with compact QC archival. Source metadata was read from the completed
producer traces and task directories; no bulk FASTQ/BAM scans or new metric jobs.
Real manifests and validation detail remain outside Git.

## Direct adversarial review

Blockers: none remaining for the configured native DNA queues. Missing library,
lane, exact reference or producer metadata remains a visible pending reason.

Warnings: full-size metric performance remains under measurement in the current
v0.1.2 jobs. New sequencing-run reports wait for every expected DNA library;
RNA-only queues, repeated chunk labels and combined multi-run BAM assignment
need explicit mapping. Failed/ambiguous submissions require operator recovery
rather than automatic duplicate submission. No universal runtime validation.

Notes: save intent before sbatch; global flock and per-run state make ticks
idempotent. Current jobs are externally managed and never resubmitted. Existing
analyses form the activation baseline. The tracker collector stays read-only;
the coordinator owns submission/delivery and the existing five-minute sync.
Source versions come from completed producer tasks. QC publication clears no
identity/scientific hold. Full-size metric completion is not claimed here.
