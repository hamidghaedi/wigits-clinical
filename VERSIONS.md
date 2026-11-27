# Release bindings

Repository version: **v0.3.0** (`VERSION`), adding the tested independently
optional blood CH and RNA-informed neoantigen branches to v0.2.0. Native
Oncoanalyser/source/reference bindings remain unchanged. See
[release notes](docs/releases/v0.3.0.md) and [Fir execution](docs/FIR_EXECUTION.md).

| Component | Bound version / identity |
|---|---|
| Maintained native intake and optional recovery | **native-intake/v1.2.0**; existing launcher selects native, neo-predict or neo-review; retains v1.1.5 record-preserving VEP sorting and explicit container Python isolation |
| Upstream workflow | nf-core/oncoanalyser **3.0.0**, commit `7c74c87a43749952b38c9a18915947570f0595a0` |
| Default installed source | **3.0.0-qlmp.2.2**, commit `05b8ec3ac953248d1f025528f03fe647c6f8ce30`; native analysis/eligibility plus [ORANGE source-value mapping](docs/ORANGE_SOURCE_VALUES.md) |
| HMF reference bundle | **GRCh38, 3.0.0--8**, with the release's GRCh38_hmf genome/index configuration |
| Site Nextflow / Apptainer modules | **26.04.4 / 1.3.5**, as selected in `workflow/conf/frontenac_paths.env` |
| DNA lane/sample resource policy | **input-v4**, from `input-resources/v0.6.0`; BWA64CPUs, lane-local RAM/walltime, sample-local REDUX |
| RNA resource policy | **input-v5**, `85043737ebda623a2dc1e6f0b2a37e95cb60cb38`; STAR, sort, duplicate marking and ISOFOX use their actual task input sizes; WGS-only packages omit RNA selectors |
| RNA optical-duplicate comparison policy | `rna-duplicates/v0.1.0`; ordinary duplicate marking retained |
| Read-only allocation/queue reporting | `resource-report/v0.2.1` |
| Native somatic output review | `native-qc/v0.1.0` |
| Bounded indexed SNP support | `snp-support/v0.1.0` |
| Laboratory sequencing-run QC | `lab-qc/v0.4.1`; completed-BAM events (no polling),32-thread Parabricks4.6.0-1 coverage / Picard insert size / fastp1.0.1 / MultiQC1.32; Phase1 presentation |
| Explicit runtime recovery | `runtime-recovery/v0.3.0`; select only for the diagnosed failure and original attempt |

The analytical component/container versions are bound by that upstream source's
native module configuration, not independently upgraded in this wrapper:
REDUX2.0.5, SAGE5.0.2, PAVE1.9, AMBER4.3, COBALT3.0, ESVEE2.0.1,
PURPLE4.4, LINX2.3.1 and ORANGE5.0.1. These bindings were also observed in
the completed development-cluster source-version records. A run's emitted
software versions, container identities and explicit overrides remain its own
source of truth; do not infer a historical run's tools from this current table.

Upstream currency checked November26,2026: GitHub's official latest stable
release API returns **3.0.0**, published September1,2026. This checks release
metadata, not universal compatibility or analytical performance.
[Official release](https://github.com/nf-core/oncoanalyser/releases/tag/3.0.0),
[reference-resource instructions](https://nf-co.re/oncoanalyser/3.0.0/docs/usage/#reference-data-urls).

DeepSomatic, native-context export, custom-reference BAM intake and the older
QC-failure eligibility renderer remain explicitly chosen routes/experiments.
Their component/model/source versions belong to their own issued packages;
none replaces the native SAGE caller by default. The general new-FFPE route
remains undecided; existing batch-specific authorizations retain their scope.

Site launch settings default to Frontenac. Fir has an explicitly selected
profile and owner-operated launcher; its first full target-site run remains
a runtime validation, separate from the completed Frontenac feature tests. Paths come from the
existing site config and directory contract. No reference
bundle, container, installed pipeline or executing package is rebuilt or
updated by publishing this version.

---

# Pinned versions — neoantigen discovery layer

Everything the layer depends on, with where it comes from and how it was verified. A result produced with
a different combination is a different analysis and is not pooled with these.

Recorded 2026-11-23 against case KQ0004 on Frontenac (CAC, Queen's University).

## Upstream pipeline

| Component | Version | Source |
|---|---|---|
| nf-core/oncoanalyser | 3.0.0 (commit `7c74c87a43749952b38c9a18915947570f0595a0`) | nf-core |
| HMF reference bundle | `hmf_pipeline_resources.38_v3.0.0--8` | HMF, staged under `data/references/oncoanalyser/3.0.0/` with `SHA256SUMS` |
| Reference genome | `GRCh38_masked_exclusions_alts_hlas.fasta` (HMF masked GRCh38) | HMF bundle — used for calling, LILAC, VEP (`--use_given_ref`) and IGV review alike |
| Nextflow | 26.04.4 | module |
| Apptainer | 1.3.5 | module |

### WiGiTS / hmftools components this layer reads

SAGE 5.0.2 · PAVE 1.9 · PURPLE 4.4 · LINX 2.3.1 · LILAC 2.0 · ISOFOX 2.0.1 · NEO 1.3 · ORANGE 5.0.1
(versions as shipped in oncoanalyser 3.0.0; the full suite is recorded in the run's `pipeline_info/`).

## Annotation

| Component | Version | Note |
|---|---|---|
| Ensembl VEP | 115.2 (container `ensembl-vep_release_115.2.sif`) | |
| VEP cache | 115, **merged** (Ensembl + RefSeq), GRCh38 | the merged cache is BAM-edited, so `--use_given_ref` is **mandatory** or VEP rejects the HMF masked reference |
| VEP plugins | `Wildtype`, `Frameshift` (pVACtools-supplied) | installed to `data/knowledgebases/vep_plugins/pvactools/` |
| Required VEP flags | `--symbol --terms SO --tsl --biotype --transcript_version --hgvs --flag_pick` | pVACseq rejects input lacking `--tsl`/`--biotype` |

## Prediction

| Component | Version | Note |
|---|---|---|
| pVACtools | 5.4.0 (`pvactools-5.4.0.sif`, 12 GB) | pvacseq used; pvacfuse/pvacsplice not used in this version |
| MHCflurry | 2.0.6 | both `MHCflurry` (affinity) and `MHCflurryEL` (eluted-ligand) |
| MHCnuggets | 2.4.1 | class I (`MHCnuggetsI`) |
| TensorFlow | 2.15.1 | MHCflurry backend, inside the pVACtools image |
| htslib (`bgzip`, `tabix`) | 1.22.1 (`samtools-1.22.1--h96c455f_0.sif`) | **neither** the pVACtools nor the VEP image ships `bgzip`; it is borrowed from the samtools image |

**Not available — the NetMHC family.** NetMHCpan, NetMHCpanEL and NetMHCIIpan are listed by pVACtools but
require an IEDB standalone licence that is not held for this project. The ensemble is therefore
MHCflurry + MHCflurryEL + MHCnuggetsI. The IEDB **web** API is not used either: patient-derived peptide
sequences are not sent to an external service. This is the single largest methodological limitation of
the layer, because the NetMHC family is the most extensively benchmarked, and it is stated in the review
SOP and in any report the layer produces.

### Known defect in pVACtools 5.4.0

The 5.4.0 container is missing its own `supporting_files/cancer_census_hotspot_gene_list.tsv`, so
`pvacseq run` aborts in post-processing **after** the predictions are complete. The layer therefore runs
`binding_filter`, `coverage_filter`, `transcript_support_level_filter`, `top_score_filter` and
`generate_aggregated_report` as separate calls on the completed `all_epitopes.tsv`. Recorded here so the
workaround is not mistaken for a stylistic choice.

## Knowledgebases

| Resource | Version | Use |
|---|---|---|
| Cancer Gene Census | v104 (768 symbols) | `gene_of_interest` flag on candidates; not a filter |
| OncoTree / CKB / CIViC etc. | see `docs/evidence/` | evidence layer, not used by this layer |

## Layer code

| Component | Identity |
|---|---|
| Interface contract | `docs/neoantigen_interface.md` |
| Tier rules | `docs/neoantigen_extension_proposal.md` §3.4, printed by `neo_candidates.py --explain` |
| `rules_version` stamped into outputs | `docs/neoantigen_extension_proposal.md@2026-11-23` |
| Branch | `neo/discovery-layer` |

## Platform

Frontenac (CAC, Queen's University): Slurm, accounts `def-hpcg6049_cpu` / `def-hpcg6049_gpu`; compute nodes
have outbound internet; node-local scratch `$SLURM_TMPDIR` (462 GB) is used for `TMPDIR` and
`APPTAINER_TMPDIR`. Containers are shared from
`containers/NXF_SINGULARITY_CACHEDIR/`.
