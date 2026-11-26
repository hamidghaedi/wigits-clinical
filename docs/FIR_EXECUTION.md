# Running a patient on Fir

Use the same tools for every patient. The preparer scans the **original lab
manifest files** configured for this project and finds the tumour, matched
normal and any declared RNA. You do not choose run manifests or assemble IDs.

## 1. Prepare

```bash
mohcc=/home/ghaedi/projects/def-ghaedi/ghaedi/mohcc
python3 "$mohcc/scripts/prepare_native_manifest.py" PATIENT_ID
```

Use the patient ID (for example `KQ0061`) or a full laboratory sample ID.
The tool prints the resolved IDs, RNA availability, prepared directory and the
exact command to submit the workflow.

Optional analysis choices: add `--blood-ch` for blood CH or
`--rna-neoantigens` for RNA-informed neoantigens. These do not change how the
sample identities are discovered. Use `--dna-only` to explicitly omit RNA.

The generated `case_manifest.json` is the workflow input SSOT. It retains the
original lab values, identity mappings, runs, specimen details, kits/dates,
file paths, options and source/version evidence. The native samplesheet and
runtime configuration are generated from those resolved bindings. Do not edit
these files independently; preparation never modifies the original lab files.

Known Phase1 separator/zero-padding differences and sourced aliases are applied
without conflating blood DNA, tumour DNA, RNA or numbered aliquots. An ambiguous
binding is reported rather than guessed. All recorded lanes and both mates are
required for DNA. Declared RNA without staged files remains recorded as
`NOT_STAGED`; it is not relabelled as absent from the laboratory record.

FASTQs stay in the canonical location:

```text
/scratch/ghaedi/mohcc/fastq/<original_run_id>/<sample_directory>/
```

## 2. Run

Copy the submission command printed by preparation. Its form is:

```bash
bash "$mohcc/scripts/run_native_fir.sh" \
  --manifest /printed/prepared/directory/samplesheet.csv \
  --config /printed/prepared/directory/site.env
```

It prints the Slurm head job ID and run directory. Check progress with
`squeue -u "$USER"`; head logs are in the printed run's `logs/` folder.
Native results are in `results/<patient_id>/`; optional CH/neo outputs are in
`extensions/results/<patient_id>/`. After correcting a failed job, repeat the
same submission command with `--resume`. A different analysis needs a new
`--analysis-id`, not edits inside an existing execution.

## Laboratory QC

QC uses **every expected DNA library from the full original run manifest**,
independently of the patient selected above. Each BAM-ready event updates that
run; the final expected library triggers its run-level MultiQC once, before
variant calling finishes. No polling or partial-run completion is used.
Full-run membership/state is in `file_manifests/lab_qc/`; successful reports
are retained under project `<phase>/qc_reports/<run_id>/` before central return.

The shared configuration is `config/native_fir.json`. Its `laboratory_sources`
points to the original lab directories and their phase bindings; the maintained
alias registry contains source-confirmed corrections. New original lab files
go in those source directories, so later preparation discovers them automatically.
Routine native analysis uses fresh-frozen tumours; new FFPE analysis still needs
its separately agreed analytical route. No software or resources are installed
by either command.
