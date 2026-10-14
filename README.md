# wigits-clinical

An open implementation of a clinical reporting layer on top of the Hartwig WiGiTS pipeline. The aim is to take the research-grade output that `nf-core/oncoanalyser` already produces and add back the pieces that turn it into a signed clinical report, which Hartwig keeps in a private repository.

## Why

Hartwig's WiGiTS tools (PURPLE, LINX, SAGE, ORANGE and the rest) are open source and well maintained, and oncoanalyser runs them end to end. What stops at the ORANGE JSON is everything clinical. Evidence matching against knowledgebases, on-label and off-label reasoning, germline consent handling, QC-fail reports, the clinical summary sentences, and the final PDF all live in `hartwigmedical/oncoact`, which is not public. This repository rebuilds that layer in the open, using what can still be read from hmftools history and from Hartwig's published regulatory documents.

## What is here now

The project starts as documentation. Nothing runs yet.

| Path | What it is |
|---|---|
| `docs/hartwig_wgts_dataflow.md` | How a default WiGiTS WGTS run is put together. Stages, data sources, the dependency graph, and how ORANGE assembles its report. |
| `docs/hartwig_clinical_layer.md` | What is known about SERVE, PROTECT, ROSE and the patient reporter from public history and regulatory documents. |
| `docs/djerba_wgts_dataflow.md` | The same analysis for OICR's Djerba, kept here as the comparison point that motivated this work. |
| `docs/gap_register.md` | The running list of features that are missing or private, with what is known and what we intend to do about each. |
| `workflow/` | Where the Nextflow layer will live. Empty apart from a note on intended structure. |
| `CONTRIBUTING.md` | How work moves between this repository and the HPC side where testing happens. |

## Intended shape

Three parts, in the order they will be built.

1. **Evidence layer.** A SERVE-compatible actionability directory built from knowledgebases we are licensed to use, and a PROTECT reimplementation that reads the ORANGE JSON.
2. **Interpretation layer.** A ROSE-style rules table and a germline consent model, both versioned in this repository.
3. **Reporting layer.** A reporter that reads the ORANGE JSON plus the two layers above and writes a PDF and a JSON record, including the QC-fail report variants.

All three sit downstream of oncoanalyser and never modify it. A thin Nextflow wrapper in `workflow/` will run oncoanalyser and then the three layers.

## Status

Documentation only. See `docs/gap_register.md` for what is open.
