# Release bindings

Repository-wide release: **v0.1.0** (`VERSION`). This is the first consolidated
release, incorporating the published `hartwig_adoptation` history through
`7a8f207`. Component versions below have independent namespaces; v0.1.0 does
not rename or replace their immutable tags.

| Component | Bound version / identity |
|---|---|
| Upstream workflow | nf-core/oncoanalyser **3.0.0**, commit `7c74c87a43749952b38c9a18915947570f0595a0` |
| Default installed source | **3.0.0-qlmp.2.2**, commit `05b8ec3ac953248d1f025528f03fe647c6f8ce30`; native analysis/eligibility plus [ORANGE source-value mapping](docs/ORANGE_SOURCE_VALUES.md) |
| HMF reference bundle | **GRCh38, 3.0.0--8**, with the release's GRCh38_hmf genome/index configuration |
| Site Nextflow / Apptainer modules | **26.04.4 / 1.3.5**, as selected in `workflow/conf/frontenac_paths.env` |
| Lane/sample resource policy | `input-resources/v0.6.0`, **input-v4**; BWA64CPUs, lane-local RAM/walltime, sample-local REDUX |
| RNA optical-duplicate comparison policy | `rna-duplicates/v0.1.0`; ordinary duplicate marking retained |
| Read-only allocation/queue reporting | `resource-report/v0.2.1` |
| Native somatic output review | `native-qc/v0.1.0` |
| Bounded indexed SNP support | `snp-support/v0.1.0` |
| Laboratory sequencing-run QC | `lab-qc/v0.3.0`; completed-BAM output event dispatch (no polling) and separate Nextflow entry, Picard image3.1.0 / fastp1.0.1 / MultiQC1.32; execution versions remain task-recorded |
| Explicit runtime recovery | `runtime-recovery/v0.3.0`; select only for the diagnosed failure and original attempt |

The analytical component/container versions are bound by that upstream source's
native module configuration, not independently upgraded in this wrapper:
REDUX2.0.5, SAGE5.0.2, PAVE1.9, AMBER4.3, COBALT3.0, ESVEE2.0.1,
PURPLE4.4, LINX2.3.1 and ORANGE5.0.1. These bindings were also observed in
the completed development-cluster source-version records. A run's emitted
software versions, container identities and explicit overrides remain its own
source of truth; do not infer a historical run's tools from this current table.

Upstream currency checked November23,2026: GitHub's official latest stable
release API returns **3.0.0**, published September1,2026. This checks release
metadata, not universal compatibility or analytical performance.
[Official release](https://github.com/nf-core/oncoanalyser/releases/tag/3.0.0),
[reference-resource instructions](https://nf-co.re/oncoanalyser/3.0.0/docs/usage/#reference-data-urls).

DeepSomatic, native-context export, custom-reference BAM intake and the older
QC-failure eligibility renderer remain explicitly chosen routes/experiments.
Their component/model/source versions belong to their own issued packages;
none replaces the native SAGE caller by default. The general new-FFPE route
remains undecided; existing batch-specific authorizations retain their scope.

Paths come from the existing site config and directory contract. No reference
bundle, container, installed pipeline or executing package is rebuilt or
updated by publishing this version.
