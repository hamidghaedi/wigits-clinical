# Running a patient on an execution worker

Use the same reusable tools on Fir and Narval. Set `mohcc` to the persistent
project root on that worker (Narval: `/project/def-ghaedi/$USER/mohcc`).
FASTQs remain on scratch in `fastq/<original_run>/<laboratory_sample>/`.

```bash
mohcc=/project/def-ghaedi/$USER/mohcc
python3 "$mohcc/scripts/prepare_native_manifest.py" PATIENT_ID
```

The preparer finds tumour, normal and declared RNA from the original laboratory
CSVs in `file_manifests/laboratory_inputs/`. It prints resolved IDs, RNA status
and the exact submission command. For a full DNA/RNA analysis, RNA must be
`READY`; add `--blood-ch --rna-neoantigens` to enable those eligible branches.
Neither option changes the native caller. New FFPE analysis retains its separate
route decision. Missing/ambiguous identities, lanes or mates remain explicit.

Run the command printed by preparation. On Narval its form is:

```bash
bash "$mohcc/scripts/run_native_workflow.sh" \
  --manifest /printed/prepared/directory/samplesheet.csv \
  --config /printed/prepared/directory/site.env
```

The launcher prints the Slurm head ID and scratch run directory. Results,
work, logs and state stay in `analysis/WiGiTS/<version>/<analysis_id>/`;
source stays in project `releases/` and shared assets in project `resources/`.
Check with `squeue -u "$USER"`. Use the same command with `--resume` after a
corrected failure; do not change the frozen package or rerun completed work.

`case_manifest.json` is the generated input SSOT; its CSV/site/options files
are projections. Do not edit them independently. All expected DNA libraries in
the complete original sequencing manifest govern laboratory QC. BAMTOOLS-ready
events trigger one run-level report after the final library, without polling or
waiting for variant calling. The per-run registry is `file_manifests/lab_qc/`.
Neither command installs software, copies inputs or deletes sources. Owner-run
target execution is required before claiming a new site's runtime validated.

## Controller walltime

New submissions default to96hours (`4-00:00:00`), based on9 completed native
controllers at11.7–71.8hours, including downstream job queue waits. The head's
own pre-start queue is excluded. Recovery-only heads and failed/current attempts
are not used to choose the default. Historical task caching remains represented.
The extra24.2hours is a starting margin, not a guarantee on a new cluster or for
full optional CH/neo runs. Individual task resources are unchanged.

Set `controller_walltime` in the shared site JSON before preparation, or append
`--controller-time 3-00:00:00` to the printed submission command for an explicit
per-submission override. CLI wins over the prepared site value, then the96h
fallback. Slurm validates the duration. The selected request is recorded in
`state/head_walltime.txt`. Already submitted jobs keep their allocation; changing
the shared setting or installing a new release does not update them. The native
head's direct `sbatch` fallback is also96h and accepts normal `sbatch --time`.
