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
