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
