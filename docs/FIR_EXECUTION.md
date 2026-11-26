# Running the workflow on Fir

Use the same two commands for every case: prepare the inputs from the original
laboratory manifests, then submit the workflow. Change the arguments, not the
scripts. The shared configuration is `mohcc/config/native_fir.json`.

Before starting, the FASTQs must be in:

```text
/scratch/ghaedi/mohcc/fastq/<original_run_id>/<laboratory_sample_id>/
```

Use the **complete original laboratory CSVs**, including the matched normal's
manifest when it was sequenced in a different run. Export an Excel manifest to
CSV first. Do not make a patient-only copy of a run manifest.

## 1. Prepare the workflow manifest

Set these values for the case you want to analyse:

```bash
mohcc=/home/ghaedi/projects/def-ghaedi/ghaedi/mohcc
case_id=PATIENT_ID
phase=Phase_1
analysis_id="${case_id}-native-$(date +%Y%m%d)-r1"

python3 "$mohcc/scripts/prepare_native_manifest.py" \
  --lab-manifest "$mohcc/file_manifests/laboratory_inputs/TUMOUR_RUN.csv" \
  --lab-manifest "$mohcc/file_manifests/laboratory_inputs/NORMAL_RUN.csv" \
  --case-id "$case_id" --phase "$phase" --analysis-id "$analysis_id" \
  --config "$mohcc/config/native_fir.json" \
  --blood-ch
```

If tumour and normal share one manifest, supply it once. Add another
`--lab-manifest` for RNA when it comes from another run.

- Omit `--blood-ch` when you do not want blood CH analysis. It requires a normal
  recorded as blood in the laboratory manifest.
- Add `--rna-neoantigens` for the optional RNA-informed neoantigen analysis.
  Recorded tumour RNA and its FASTQs must be available.
- Use `--dna-only` if you explicitly want to omit RNA.
- If there are several tumour or normal libraries, select them with
  `--tumour-sample` / `--normal-sample` using the exact laboratory sample IDs.

The command writes `file_manifests/<analysis_id>/samplesheet.csv`, `site.env`
and the supporting metadata. It includes every recorded lane, both mates and
all observed FASTQ chunks. Missing inputs or ambiguous laboratory bindings stop
preparation. It does not transfer data, count all reads or submit any jobs.
The routine route is native Hartwig for fresh-frozen tumours; the new FFPE
route still requires a separate analytical decision.

## 2. Submit the workflow

```bash
bash "$mohcc/scripts/run_native_fir.sh" \
  --manifest "$mohcc/file_manifests/$analysis_id/samplesheet.csv" \
  --config "$mohcc/file_manifests/$analysis_id/site.env"
```

The configuration carries the options selected in step1. The launcher prints
the Slurm job ID and output directory. Check progress with `squeue -u "$USER"`;
head logs are in the printed run directory's `logs/` folder.

Outputs are under
`/scratch/ghaedi/mohcc/analysis/WiGiTS/<version>/<analysis_id>/`:

| Location | Contents |
|---|---|
| `results/<case_id>/` | Native analysis, alignments and ORANGE outputs |
| `extensions/results/<case_id>/` | Selected blood CH / RNA neoantigen outputs |
| `extensions/state/status.json` | Separate optional-analysis status |
| `results/pipeline_info/resource_usage/` | CPU, memory, runtime and queue reports |

After correcting a failed job, repeat **step2** with `--resume`. Keep the same
manifest, configuration and work directory. Do not repeat preparation or change
an existing execution's inputs; use a new analysis ID for a different analysis.

## Laboratory run QC

The preparer separately registers the **whole DNA sequencing-run membership**
from each original laboratory manifest. A selected patient's samplesheet never
replaces that full-run list. Tumour and normal from different runs belong to
separate run-level QC reports.

Each completed native BAMTOOLS task sends a BAM-ready event. When every expected
DNA library and lane is represented, the existing laboratory QC workflow starts
once for that run, without waiting for variant calling. Missing libraries keep
it pending; no timer polls for readiness and no partial run is labelled complete.
Run state and full registered manifests are in `file_manifests/lab_qc/`.
Successful reports are retained under project `<phase>/qc_reports/<run_id>/`;
central return and Laptop_sync publication retain their own delivery receipts.
