# Directory layout contract

**Status: canonical.** This document is the single authority for where things live on the HPC side of
`wigits-clinical`. Every script, config and document in this repository refers to these locations and to
nothing else. If a path is not listed here, it is not a sanctioned location; add it here first, then use it.

Machine-readable form: [`workflow/conf/frontenac_paths.env`](../workflow/conf/frontenac_paths.env) — one
variable per root below. Scripts `source` that file; they do not hard-code paths. The two files are kept
in step: a change to one is a change to both, in the same commit.

Cluster: Frontenac (CAC, Queen's). Account `hpc6049`, group `hpcg6049`. Verified 2026-10-15.

## 1. Principles

1. **Three tiers, three filesystems.** Versioned text lives in this git repository. Shared, immutable,
   large assets (containers, reference data, knowledgebases) live once on project storage and are shared
   with the sibling pipelines. Mutable, large, regenerable run data lives on scratch.
2. **Reuse the sibling layout.** `wigits-clinical` sits beside `somatic_pipeline_sm` and `dnaseq_nf` under
   the same `somatic_pipeline/` parent and shares their container cache and reference root. We do not
   create a second copy of anything they already hold.
3. **Nothing permanent lives in a home directory or in `/tmp`.** Home is per-user and small; `/tmp` is
   node-local. The one current exception (the Nextflow pipeline clone under `~/.nextflow`) is listed as a
   known gap in §5.
4. **Every durable directory carries a `README.md`** saying what it is, which version it is, where it
   came from, and a checksum file where the contents are downloaded.
5. **Immutable once released.** Reference bundles and released run outputs are never edited in place; a
   new version gets a new directory.

## 2. The roots

| Variable (`frontenac_paths.env`) | Path | Tier | Purpose |
|---|---|---|---|
| `WIGITS_REPO` | `/global/project/hpcg6049/somatic_pipeline/wigits-clinical` | git | this repository: code, configs, docs, samplesheets, test summaries |
| `PIPELINE_ROOT` | `/global/project/hpcg6049/somatic_pipeline` | project | shared parent of all somatic-pipeline work (sibling repos live here) |
| `CONTAINER_DIR` | `${PIPELINE_ROOT}/containers/NXF_SINGULARITY_CACHEDIR` | project, shared | **existing.** All Apptainer images, for every pipeline. oncoanalyser's images land here as `quay.io-*.img` |
| `REF_ROOT` | `${PIPELINE_ROOT}/data/references` | project, shared | **existing.** All reference data, one subdirectory per tool/bundle version |
| `ONCO_REF_DIR` | `${REF_ROOT}/oncoanalyser/<oncoanalyser-version>` | project, shared | **new.** The complete reference set oncoanalyser `<version>` pins (genome, indexes, HMF resources). See §3 |
| `KB_ROOT` | `${PIPELINE_ROOT}/data/knowledgebases/<source>/<version>` | project, shared | **new, reserved.** Licensed/clinical knowledgebases for the evidence layer (CKB, CIViC, OncoKB, ...). Never in git |
| `RUN_ROOT` | `/global/scratch/hpc6049/wigits_runs` | scratch | **new.** One directory per run (§4). Work dirs, intermediate files, Slurm logs. Regenerable |
| `RELEASE_ROOT` | `${PIPELINE_ROOT}/wigits_releases` | project | **new.** Durable copy of what a finished run produced that we keep: ORANGE JSON/PDF, clinical-layer outputs, provenance. No BAM/CRAM |
| `SOFTWARE_ROOT` | `${PIPELINE_ROOT}/software` | project | **new, reserved.** Pinned checkouts of external pipelines (§5) |
| `ONCO_SOURCE_DIR` | `${SOFTWARE_ROOT}/oncoanalyser/3.0.0-qlmp.1` | project | Isolated, tagged local QC-report fix; upstream asset cache stays unchanged |

`ONCO_SOURCE_VERSION=3.0.0-qlmp.1` and `ONCO_REFERENCE_VERSION=3.0.0` deliberately
distinguish the patched reporting interface from the unchanged native reference
bundle. The launcher defaults to unmodified native 3.0.0. The patched source
is experimental and requires an explicit `ONCO_SOURCE` override; see the workflow runbook.

Agent scratch (`/global/home/hpc6049/.claude-science-scratch`, `/global/scratch/hpc6049/.claude-science-scratch`)
is harness-internal: per-job wrappers and logs for jobs an assistant submitted. Nothing we rely on lives there.

## 3. Reference data — `ONCO_REF_DIR`

oncoanalyser pins one reference set per release; we store exactly that set, keyed by the oncoanalyser
version, so a `-r <version>` run has one unambiguous reference directory.

```
data/references/oncoanalyser/
  3.0.0/
    README.md                         source URLs, download date, how it was staged
    SHA256SUMS                        of every file below
    reference_data.config             Nextflow config pointing oncoanalyser at these files
    genomes/GRCh38_hmf/
      GRCh38_masked_exclusions_alts_hlas.fasta (+ .fai, .dict, .img)
      bwa-mem2_index-2.2.1/
      gridss_index-2.13.2/
      star_index-gencode_38-2.7.3a/
    hmf_pipeline_resources.38_v3.0.0--8/
      common/  dna/  rna/  misc/
```

Rules: populated once by `--mode prepare_reference` (or verified manual download), then read-only
(`chmod -R a-w`). Never symlinked from a run `work/` directory. The older
`${REF_ROOT}/hmf_pipeline_resources.38_v2.3.0--2/` belongs to `dnaseq_nf` (AMBER/COBALT/PURPLE 4.x) and
is not an oncoanalyser 3.0.0 input; it stays where it is.

Current state (2026-10-15): the 3.0.0 set exists only inside the smoke-test work dir
(`/global/scratch/hpc6049/oncoanalyser/test_3.0.0/work/`). Staging it into `ONCO_REF_DIR` is the next
step; until then that work dir must not be deleted.

## 4. Runs — `RUN_ROOT` and `RELEASE_ROOT`

```
/global/scratch/hpc6049/wigits_runs/
  <run_id>/
    samplesheet.csv          byte-identical copy of the committed samplesheet
    launch.sh                the exact command used (copied by the launcher)
    work/                    Nextflow work dir (-work-dir)
    results/                 oncoanalyser --outdir, then clinical-layer outputs
    logs/                    Slurm .out of the head job, .nextflow.log copies
```

`<run_id>` = `<label>_<YYYYMMDD>_onco<version>`, e.g. `test_20261014_onco3.0.0`,
`KQ0315_20261101_onco3.0.0`. `<label>` is a case id, cohort name or `test`.

Scratch is volatile. Anything a run produced that must outlive it is copied to
`${RELEASE_ROOT}/<run_id>/` together with `samplesheet.csv`, `launch.sh`, `pipeline_info/` and a
`RELEASE.md` naming the oncoanalyser version, pipeline commit, reference directory and container cache
state. BAM/CRAM are not released here; alignments stay on scratch or in the existing sibling archives.

The smoke test currently at `/global/scratch/hpc6049/oncoanalyser/test_3.0.0/` predates this contract; it
is the legacy location of `test_20261014_onco3.0.0` and is retired once §3 staging has consumed it.

## 5. Software — `SOFTWARE_ROOT`

Today the oncoanalyser clone lives in `~/.nextflow/assets/` because `nextflow pull` put it there. That is
a user home directory and violates principle 3. Target state:

```
software/
  nf-core-oncoanalyser/
    3.0.0/        git checkout of tag 3.0.0 (commit 7c74c87a43)
```

with the launcher calling `nextflow run ${SOFTWARE_ROOT}/nf-core-oncoanalyser/3.0.0 ...` so a run never
depends on network access to GitHub or on one user's home. Nextflow itself and Apptainer come from
cluster modules (`nextflow/26.04.4`, `apptainer/1.3.5`), recorded in `VERSIONS.md` once it exists.

## 6. In this repository

```
wigits-clinical/
  docs/                      design documents, this contract, gap register
  workflow/
    conf/
      frontenac.config             Slurm + Apptainer site config
      frontenac_limits.config      production resource ceilings
      frontenac_paths.env          THE paths (this contract, machine-readable)
      reference_data.config        copy of ONCO_REF_DIR/<ver>/reference_data.config once staged
    run_oncoanalyser_test_frontenac.sh
    main.nf, modules/, assets/     clinical layers (not yet written)
  samplesheets/
    <run_id>.csv                   oncoanalyser samplesheet for each run, as submitted
    <run_id>.csv.sha256
  tests/results/<name>/            summaries only; outputs stay on the HPC (CONTRIBUTING.md)
```

Samplesheets are the run manifests of this project and are versioned here, mirroring
`dnaseq_nf/run_manifests/`. They contain sample identifiers and HPC file paths only — no clinical fields.

## 7. What already exists where (inventory, 2026-10-15)

| Thing | Where it is now | Contract location | Action |
|---|---|---|---|
| oncoanalyser 3.0.0 code | `~/.nextflow/assets/.repos/nf-core/oncoanalyser/clones/7c74c87a43.../` | `SOFTWARE_ROOT/nf-core-oncoanalyser/3.0.0` | clone to target; update launcher |
| 29 tool images (14 GB) | `CONTAINER_DIR/quay.io-*.img` | same | none |
| GRCh38 reference set for 3.0.0 (~68 GB) | smoke-test `work/` | `ONCO_REF_DIR/3.0.0` | stage (`prepare_reference`), checksum, write `reference_data.config` |
| HMF resources v2.3.0--2 | `REF_ROOT/hmf_pipeline_resources.38_v2.3.0--2` | same (dnaseq_nf asset) | none |
| smoke-test outputs (157 MB) + work (123 GB) | `/global/scratch/hpc6049/oncoanalyser/test_3.0.0` | `RUN_ROOT/test_20261014_onco3.0.0`; summary already in `tests/results/` | move after staging; release ORANGE outputs to `RELEASE_ROOT` |
| site configs, launcher | `WIGITS_REPO/workflow/` | same | source `frontenac_paths.env` instead of hard-coding |

## 8. Open questions

- Scratch retention policy on Frontenac is not documented in the sibling project; until CAC confirms it,
  treat `RUN_ROOT` as deletable at any time and release promptly.
- Whether `RELEASE_ROOT` should hold full per-case oncoanalyser `results/` (a few GB per case without
  alignments) or only the ORANGE outputs and clinical-layer reports. Default until decided: full `results/`
  minus alignments.
