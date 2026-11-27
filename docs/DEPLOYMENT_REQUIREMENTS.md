# Deploying the native workflow on another cluster

Reviewed 27 November 2026 against native-intake/v1.1.4 (`295c0f1`).
The technical repairs below are now published in **native-intake/v1.1.5**;
see [release and recovery commands](releases/native-intake-v1.1.5.md) and the
[verification receipt](../tests/results/optional-portability/README.md).
Current reusable operation uses the same launcher's explicit
[v1.2.0 entry points](releases/native-intake-v1.2.0.md), retaining those fixes.
Frontenac's successful executions are the reference deployment. Reuse their
exact source, images, resources and selected options; put cluster differences
in the shared site configuration. A successful source-site run does not establish
that the target has the same Python environment, mounts or predictor dependencies.

## Requirements to record once per deployment

| Layer | Requirement | Evidence before routine use |
|---|---|---|
| Controller/cluster | Existing Bash, Python, Nextflow with its Java runtime, Apptainer and authenticated Slurm access | Record actual executable/module versions and account/partition from the site binding; run a small compute-node task through that binding. Fir has executed with Nextflow26.04.4/Apptainer1.3.5; these are observed versions, not newly established minimum versions. |
| Host laboratory-QC helper | Exact configured Python can import PyYAML; existing helper imports must work in its noninteractive job environment | Fir/Narval binding is `/home/ghaedi/miniforge3/bin/python`, with owner-installed PyYAML6.0.3. Check that interpreter, rather than installing into an unrelated shell or user site. |
| Analytical containers | The release's exact images and their required runtime libraries, accessible on compute nodes | Reuse the recorded image inventory/transfer receipts. Exercise the selected tools in containers on the target; a file landing or CLI version alone does not test a predictor. |
| References and predictor data | Pinned native genome/index/bundle, selected VEP cache/plugins, HLA/predictor resources and models | Reuse the deployment inventory. Run the selected predictor on a tiny supported input using the staged models/cache; do not assume a successful Python import proves prediction works offline. |
| Paths and isolation | Canonical project/scratch trees, input symlink targets visible through container binds, writable output/work/cache | Exercise the actual configured binds and environment. Record Python/library origins when investigating a container import failure. Avoid inherited host or user-site packages supplying container dependencies. |
| Optional features | Native analysis, blood CH and additional RNA/pVACseq prediction have separate readiness/status | Test each selected feature. Preserve successful native/CH products when an optional prediction fails; failure never means a negative biological result. |

These are deployment checks, not new per-patient gates. No repeated BAM/FASTQ
scans, hashes, read counting or polling are added. Existing receipts are reused.
Run compute tests under Slurm; preserve running jobs and frozen packages.

For the currently configured host QC interpreter, the lightweight import check is:

```bash
env -u PYTHONPATH PYTHONNOUSERSITE=1 /home/ghaedi/miniforge3/bin/python -c \
  'import sys, yaml; print(sys.executable); print(yaml.__version__)'
```

Run it in the same job environment as the helper. Fir's successful BAM-ready
replays already demonstrate this path works. Narval's Conda package record and
configuration have been observed; an actual target import/runtime remains pending.
The existing dispatcher's `--bam-ready` replays recover missed events; neither
alignment nor completed metrics need repeating.

## Known issues: keep the categories separate

| Observation on Fir | Category and current disposition |
|---|---|
| QC dispatcher could not import `yaml` | Host Python dependency/binding. Corrected interpreter configuration; both missed notifications replayed successfully. |
| Additional pVACseq prediction reports NumPy1.24.2 failing to load `libflexiblas.so.3` with Python3.11 at `/usr/local/bin/python` | Host user-site contamination is observed: the full Fir traceback imports pandas from `/home/ghaedi/.local/lib/python3.11/site-packages/`. Repaired in v1.1.5 with explicit container Python isolation. Real poisoned-user-site/PYTHONPATH probes pass; no host numerical-library installation is required by the repair. Fir's own recovery remains owner-run. |
| Additional neo prediction fails tabix with `Chromosome blocks not continuous` on the proximal VEP VCF | Adapter data-ordering defect, not an installation prerequisite. Repaired in v1.1.5 by sorting task-owned VEP output before compression/indexing, including recovery from existing exports. Real replay retains all66 selected and4 proximal input records; no hand-edited production output is used. |

The two inspected Fir cases completed native WGTS (48/48 tasks each), including
ORANGE and native NEO, and completed blood CH. Their additional RNA/pVACseq
prediction failed at the two different points above. Earlier Frontenac extension
tests, including the 69-row comparison, passed on their tested inputs/environment.
That evidence is retained; it does not resolve the later target/input failures.
Raw case logs remain in the private execution records, outside Git.

### Why the successful reference missed these failures

The Frontenac reference proximal VCF contained **zero records**. The failing
Fir export contains **three somatic records followed by one germline record**;
its chromosome blocks are `chr3, chr6, chr10, chr6`. The exporter source is
byte-identical across these executions. Appending separately ordered call sets
does not make a combined ordered VCF. The corresponding VEP output retains
the problem. The other inspected Fir case also has zero proximal records and
passes indexing before encountering the separate Python import failure.

Compute probe12516094 on Frontenac completed0:0 in4s: tabix rejected a copy of
the failing VCF (exit1), then indexed a sorted copy (exit0). All four complete
record lines, including sample values, and the header were preserved. The same
probe imported the source pVACtools image's NumPy1.26.4 successfully, using its
bundled OpenBLAS libraries. The original attempt12516093 stopped before container
execution because the cluster module setup was sourced under `set -u`; the
corrected probe temporarily relaxed nounset only for module initialization.

This identifies two missing test conditions: nonempty mixed-origin proximal
records and an existing host Python user site. `--cleanenv` clears ordinary
environment inheritance but does not disable Python's user-site lookup in a
mounted home directory. Repairs must cover these conditions without dropping
variants, changing scientific thresholds or reinstalling host numerical libraries.
That was the diagnosis checkpoint. The v1.1.5 verification now demonstrates
both repairs and a complete replay of the affected case's optional work on
Frontenac. No Fir retry or original production output change is claimed.
The separate pVACseq proximal-phasing interface remains incomplete (gap
`OP-NEO-PHASING`); technical completion does not validate phased peptide sequences.

## Containerising the remaining software: proposed next increment

Most native analytical software already runs in task containers. Keep that
Nextflow design: an image per tool/task, with shared external reference data.
Nextflow and Slurm submission remain a small host layer. A single image around
the whole controller would still need scheduler access and host bind/runtime
integration, while introducing nested-container concerns. It would not fix the
unsorted-VCF defect.

The practical next increment is to package remaining Python helpers with their
dependencies (including PyYAML), pin those images and build recipes, and make
container environment isolation explicit. The diagnosis above now supports
preserving the existing pVACtools image and excluding host user-site imports.
Confirm that repair on Fir using the image's own packages before restarting
prediction. Keep predictor models/cache pinned outside images where needed.

Test the resulting bundle once on a compute node at each site with small native
and selected-extension inputs. Record exact source/image/resource bindings and
outcomes; a different site's installation should then require only its site
configuration. This is a proposal, not an implemented helper container, newly
validated deployment or authorization to change current jobs/install packages.

Implementation inspected: `workflow/run_full_native.sbatch`,
`workflow/conf/fir.config`, `workflow/dispatch_laboratory_qc.py` and
`workflow/extension_tasks.py`. See also [optional execution](OPTIONAL_EXTENSIONS.md).

Primary references: [Nextflow task containers](https://docs.seqera.io/nextflow/container),
[Nextflow with Apptainer](https://docs.seqera.io/nextflow/container/apptainer),
[Apptainer environment handling](https://apptainer.org/docs/user/1.3/environment_and_metadata.html)
and [FlexiBLAS](https://www.mpi-magdeburg.mpg.de/projects/flexiblas).
