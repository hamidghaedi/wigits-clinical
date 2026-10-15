# nf-core/oncoanalyser 3.0.0 smoke test on Frontenac (CAC)

Date: 2026-10-14. Outcome: **pass** — `Pipeline completed successfully`, 42 tasks COMPLETED, 4 CACHED, 0 errors.

## Command

```bash
sbatch workflow/run_oncoanalyser_test_frontenac.sh
# which runs, inside a 2 CPU / 8 GB Slurm job (the Nextflow head):
nextflow run nf-core/oncoanalyser -r 3.0.0 -profile test,apptainer \
    -c workflow/conf/frontenac.config \
    -work-dir /global/scratch/hpc6049/oncoanalyser/test_3.0.0/work \
    --outdir   /global/scratch/hpc6049/oncoanalyser/test_3.0.0/results -resume
```

Every pipeline task was submitted as its own Slurm job (`--account=def-hpcg6049_cpu`, no partition).

## Versions

| Component | Version |
|---|---|
| nf-core/oncoanalyser | 3.0.0 (commit 7c74c87a43) |
| Nextflow | 26.04.4 (module `nextflow/26.04.4`; pipeline requires >= 25.10.4) |
| Apptainer | 1.3.5 (module `apptainer/1.3.5`) |
| HMF pipeline resources | `hmf_pipeline_resources.38_v3.0.0--8` (staged by the pipeline from data.oncoanalyser.com) |
| Genome | GRCh38_hmf (`GRCh38_masked_exclusions_alts_hlas.fasta`, bwa-mem2 2.2.1 index, STAR 2.7.3a index, GRIDSS 2.13.2 index) |
| Tool images | see `software_versions.yml` in this directory |

## Inputs

nf-core test samplesheet `fastq_eval.subject_a.wgts.tndna_trna.minimal.csv` (tumour/normal DNA + tumour RNA FASTQ, ~3 MB each). Mode `wgts`.

## Outputs (left on the HPC)

`/global/scratch/hpc6049/oncoanalyser/test_3.0.0/results/` — 157 MB; `subject_a/orange/subject_a.tumor.orange.{pdf,json}`,
`multiqc_report.html`, per-tool directories (purple, linx, sage, pave, esvee, lilac, cuppa, chord, sigs, neo, virusbreakend, ...).
Work dir 123 GB (includes ~50 GB of staged reference bundles).

## Issues found and fixed on the way

1. `/cvmfs/soft.computecanada.ca/config/profile/bash.sh` is not `set -u` safe — source it before `set -u`.
2. `~/.bashrc` exports `TMPDIR=/global/scratch/hpc6049/tmp`; Nextflow forwards TMPDIR into the container but binds only the
   work dir, so `sambamba sort` failed with "Read-only file system". Fixed in `workflow/conf/frontenac.config`
   (`TMPDIR=$SLURM_TMPDIR` exported and bind-mounted).
3. Image pulls (`mksquashfs`) run on the head node; with build temp on GPFS the GRIDSS image took >35 min on 2 CPUs.
   Launcher now uses node-local `APPTAINER_TMPDIR` and 4 CPUs. Images are cached in the shared
   `containers/NXF_SINGULARITY_CACHEDIR`, so later runs skip the pulls.
4. A `-c` site config outranks `-profile test`; resource ceilings therefore live in `workflow/conf/frontenac_limits.config`
   (production only) so the test profile's 4 CPU / 30 GB / 1 h clamps survive.

## Files here

- `execution_trace.txt` — Nextflow trace (status, realtime, peak RSS per task)
- `software_versions.yml` — tool versions collected by the pipeline
- `params.json` — resolved pipeline parameters
