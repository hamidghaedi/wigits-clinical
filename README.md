# wigits-clinical

An open implementation of a clinical reporting layer on top of the Hartwig WiGiTS pipeline. The aim is to take the research-grade output that `nf-core/oncoanalyser` already produces and add back the pieces that turn it into a signed clinical report, which Hartwig keeps in a private repository.

## Why

Hartwig's WiGiTS tools (PURPLE, LINX, SAGE, ORANGE and the rest) are open source and well maintained, and oncoanalyser runs them end to end. What stops at the ORANGE JSON is everything clinical. Evidence matching against knowledgebases, on-label and off-label reasoning, germline consent handling, QC-fail reports, the clinical summary sentences, and the final PDF all live in `hartwigmedical/oncoact`, which is not public. This repository rebuilds that layer in the open, using what can still be read from hmftools history and from Hartwig's published regulatory documents.

## What is here now

The maintained integration branch is **`main`**. The source version is **v0.2.0**,
consolidating manifest-complete laboratory QC and task-local RNA resource sizing
on the v0.1.0 baseline. Production uses a published numbered release or exact
commit; see [v0.2.0 notes](docs/releases/v0.2.0.md) for validation and release state.
The previous `hartwig_adoptation` history is retained. New work uses short-lived
branches and pull requests; see [CONTRIBUTING.md](CONTRIBUTING.md).

The implemented workflow launches nf-core/oncoanalyser 3.0.0's native Hartwig
analytical stack on Frontenac. The default installed source is
`3.0.0-qlmp.2.2`, with the documented ORANGE source-value corrections; it is
not described as byte-identical upstream source. It provides native original-
FASTQ execution, input-adjusted resources, read-only QC and efficiency reports.
Existing-BAM intake, DeepSomatic substitution and the QC-failed renderer are
separately selected routes/experiments. Their source being present does not
change the native default. The clinical layers described below remain proposed.

For new routine work, the owner selected native Oncoanalyser for fresh-frozen
tumours. FFPE routing remains an open decision. See
[the execution contract](docs/FULL_NATIVE_EXECUTION.md) for that policy and
the coordinator-owned native analysis/ORANGE tracker integration.

The earlier DNA-branch adaptation inputs and metadata work are now maintained
here; see [the migration record](docs/migrations/DNA_HARTWIG_ADOPTATION_20261110.md).

| Path | What it is |
|---|---|
| `docs/hartwig_wgts_dataflow.md` | How a default WiGiTS WGTS run is put together. Stages, data sources, the dependency graph, and how ORANGE assembles its report. |
| `docs/hartwig_clinical_layer.md` | What is known about SERVE, PROTECT, ROSE and the patient reporter from public history and regulatory documents. |
| `docs/djerba_wgts_dataflow.md` | The same analysis for OICR's Djerba, kept here as the comparison point that motivated this work. |
| `docs/gap_register.md` | The running list of features that are missing or private, with what is known and what we intend to do about each. |
| `docs/directory_layout.md` | **Canonical directory-layout contract.** Where every kind of file lives on the HPC side (git / shared project storage / scratch), what already exists where, and the machine-readable `workflow/conf/frontenac_paths.env`. |
| `docs/neoantigen_extension_proposal.md` | Neoantigen discovery and review extension on top of native oncoanalyser outputs (NEO, LILAC, SAGE_APPEND, ISOFOX): what is already covered, what is added, phased plan. Gap G28. |
| `docs/neoantigen_interface.md` | **Canonical oncoanalyser-to-pVACtools interface.** Every field the neoantigen layer consumes, with source file, meaning and conversion rule, verified against the KQ0004 and subject_a runs. |
| `workflow/` | Native Oncoanalyser launch profiles, optional BAM mate preparation, input-size resource planning, transfers and read-only resource reports. |
| `test_data/hartwig_adaptation/` | Migrated chromosome-test input preparation, samplesheet, metadata probes and native-run evidence. |
| `docs/FULL_NATIVE_EXECUTION.md` | Full original-FASTQ entry using the native HMF reference and all DNA lanes. |
| `docs/NATIVE_BAM_ENTRY.md` | Existing-BAM entry with explicit custom reference and optional mate preparation. |
| `docs/RESOURCE_REPORTING.md` | CSV/HTML efficiency reporting and the non-blocking future-launch completion hook. |
| `docs/LABORATORY_QC.md` | Event-driven, manifest-complete laboratory run QC: 32-thread Parabricks coverage, Picard inserts, QC-only fastp and Phase1 presentation. |
| `docs/RNA_RESOURCE_SIZING.md` | Input-v5 sizing from each RNA task's actual FASTQ/BAM size; DNA resource rules retained. |
| `docs/RESOURCE_TUNING.md` | Full-case CPU diagnosis and isolated64-CPU performance experiments before changing defaults. |
| `docs/DEEPSOMATIC_SUBSTITUTION.md` | Isolated Google caller test and the minimal somatic-VCF substitution boundary; native defaults unchanged. |
| `docs/DEEPSOMATIC_DISTRIBUTED_TEST.md` | Two-node example-generation benchmark, measured-count/native-annotation adapter and genuine SAGE context probe. |
| `CONTRIBUTING.md` | How work moves between this repository and the HPC side where testing happens. |
| `VERSIONS.md` | Current repository release, upstream/source/resource bindings and component tags. |
| `docs/RELEASING.md` | How reviewed changes reach main, and how an upstream Oncoanalyser update becomes a new tested release. |

## Intended shape

Three parts, in the order they will be built.

1. **Evidence layer.** A SERVE-compatible actionability directory built from knowledgebases we are licensed to use, and a PROTECT reimplementation that reads the ORANGE JSON.
2. **Interpretation layer.** A ROSE-style rules table and a germline consent model, both versioned in this repository.
3. **Reporting layer.** A reporter that reads the ORANGE JSON plus the two layers above and writes a PDF and a JSON record, including the QC-fail report variants.

These proposed clinical layers sit downstream of oncoanalyser. They are not
implemented or enabled by this release. Existing ORANGE source corrections and
explicit experimental patches are documented separately in `workflow/`.

## Status

Native execution and the listed operational tools are implemented and have
recorded tests/runs. This is a research workflow release, not clinical assay
validation. See [the v0.1.0 release record](docs/releases/v0.1.0.md) for exact
scope, checks and limitations, and [the gap register](docs/gap_register.md) for
unfinished work. Report templates remain owned by the reporting repositories.
