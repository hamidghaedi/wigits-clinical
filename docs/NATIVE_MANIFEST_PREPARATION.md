# Reusable native manifest preparation

Component `native-intake/v1.0.0`, on the WiGiTSv0.3.0 base. Operator instructions
are the standalone [Fir guide](FIR_EXECUTION.md); they are sample-independent.
Source is `workflow/prepare_native_manifest.py`, Python standard library only.

Input: complete original laboratory CSVs (repeat `--lab-manifest` across runs),
case ID, phase where absent in the laboratory source, shared site JSON and
optional exact library selectors / analysis ID / CH / RNA options. Excel sources
are exported to CSV. No patient-specific script, guide or hand-built samplesheet
is part of the supported process.

The explicit mappings are documented in the module docstring. The selected
samplesheet preserves native library IDs, recorded lanes, both mates and every
observed chunk. Acquisition run IDs come from laboratory paths; canonical
worker files supply filenames/sizes. Optional RunInfo XML supplies flowcell;
otherwise only one FASTQ header is read. No bulk read count, hashing, BAM scan,
transfer, scheduler submission or scientific filtering occurs during preparation.
Missing/ambiguous bindings stop preparation. Specimen facts govern eligibility;
fresh-frozen is the routine native route, while the new FFPE decision remains open.

Generated runtime inputs live at `file_manifests/<analysis_id>/`: samplesheet,
site/options, selected laboratory binding, input plan and source/release provenance.
They are data and receipts, not custom programs. The launcher reuses these inputs
and the existing native/optional entries. It never generates case-specific code.
Options selected in preparation are carried by the configuration into submission.
The historical explicit `--blood-ch` launcher assertion remains compatible.

Full DNA acquisition-run membership is separately retained under
`file_manifests/lab_qc/laboratory/<run>/<source_digest>/manifest.csv` with the
original source paths and hashes. Original laboratory inputs are never rewritten.
The shared registry points to each full run source; a patient-only CSV cannot
shrink existing scope. Normal and tumour in different runs retain separate scopes.
Dispatcher2.3.0 accepts these prepared-package roots in addition to existing queues.
Its BAM-ready events and the prior once-per-complete-run barrier are unchanged.
Multi-chunk FASTQs get unique child labels but still one summary row per library.
QC submits using the shared site environment, never the last patient configuration.

The site JSON is prepared once per release/site and reused across samples. It
binds `project_root`, `scratch_root`, exact workflow/source versions and commits,
reference root, extension resource paths and the existing QC dispatcher/stager.
Source fields point to the maintained BAM-event overlay, prepared once using
`prepare_lab_qc_event_source.py`; the analytical modules remain unchanged.
No environment, reference, model or package is installed by the preparer.

Validation covers two synthetic paired cases with different arguments, optional
RNA/CH, missing lanes/mates, multi-chunk input, ambiguous libraries, native route
eligibility, RunInfo without raw-header reads, full-run retention and refusal to
shrink it. A cross-case QC test shows no submission for a partial run and exactly
one submission after the other expected case is represented. Existing real-data
native/CH/neo and laboratory-QC evidence remains separately recorded; target-site
runtime requires an owner-launched execution.

## v1.1.0 — original-source discovery

The normal CLI is now `prepare_native_manifest.py <patient_or_sample_id>`.
The shared site profile is found automatically beside the installed scripts.
`laboratory_sources` supplies complete original CSV globs and source phase
bindings; explicit CSV/phase flags remain optional compatibility overrides.
All original sources are scanned before tumour/normal/RNA selection. Unrelated
incomplete RNA paths do not block another patient's valid pair. If several DNA
libraries exist, only a uniquely staged complete candidate can be selected
automatically; otherwise the ambiguity remains explicit.

Case-token normalization reproduces the tracker laboratory/Seqera rule.
Separator matching reproduces the safe Phase1 standardizer behavior; suffix
letters and numbered aliquots remain distinct. The broad historical pass4
(GDNA→DNA and1DNA→DNA) is deliberately excluded: previous processed manifests
contain one such potentially unsafe association. Exact source-confirmed aliases
are loaded from the maintained registry, with its evidence retained. Source
rows/IDs, actual storage names and matching method stay together.

`case_manifest.json` is the generated workflow input SSOT: declared and runnable
libraries, source-row metadata, phase, assay mode, options, run directory,
actual file bindings and version/resource provenance. CSV/options/site files are
projections from the same resolution, never independently edited masters. The
launcher preserves this record in the execution package. RNA states distinguish
NOT_DECLARED, NOT_STAGED, NOT_SELECTED, AMBIGUOUS and READY.

## v1.1.1 — execution-site bindings

Shared profiles may explicitly set `site`, `fastq_root`, `run_root`,
`manifest_root` (including `{phase}`), `qc_registry_root` and `launcher`.
Fir defaults remain unchanged. This lets Frontenac reuse its existing staged
inputs and canonical run tree without copying data or hand-editing a case
samplesheet. Optional `context_sources` locates approved clinical/runtime
sidecars by their recorded case ID. Their values/provenance are preserved;
conflicting contexts are not guessed. Original laboratory files still own
biological identities. The generated case manifest and sidecars are retained
with the actual execution package.

## v1.1.2 — Narval site binding

The canonical worker layout and original-source intake are shared with Narval.
A shared `config/native_site.json` is preferred when present; the existing
`config/native_fir.json` fallback preserves issued Fir installations. Narval
selects `site: narval`, a site-specific profile identity and the same worker
executor, task-local resource policies and release ceilings. Roots, account and
module versions are configuration inputs; no scientific module changes.
Target module/hardware/runtime verification remains pending until executed on
Narval. Staging is not runtime validation. Existing Fir packages remain pinned.
The generic owner guide is [NATIVE_EXECUTION.md](NATIVE_EXECUTION.md).

Narval site validation:179 primary offline tests pass (3 existing skips),3
adapter tests pass, shell syntax and `git diff --check` pass. These tests cover
generated WGTS/CH/neo inputs and shared root bindings; they do not establish
Narval module availability, scheduler execution or performance.

## v1.1.3 — configurable controller walltime

An optional site `controller_walltime` is retained as `WIGITS_CONTROLLER_TIME`
in generated `site.env`. The maintained launcher accepts `--controller-time`,
using CLI, prepared site, then96h default precedence. The native head's direct
Slurm fallback is96h. The default covers the longest of9 observed successful
native controller occupancies (71.8h) with24.2h margin, including child queue
waits. No task CPU, memory, walltime, inputs or analytical code changes.
Current KQ0061 remains running with its existing allocation by owner direction.

Validation for v1.1.3:181 primary tests pass (3 existing skips),3 adapter
tests pass, shell syntax and diff checks pass. Stubbed submissions verify the
default, shared-site value and explicit CLI precedence, while preserving the
prepared package and duplicate-submission behavior.
