# Running a patient on an execution worker

Before routine use on a new worker, follow the [deployment requirements](DEPLOYMENT_REQUIREMENTS.md).
Transferred assets and configured paths alone do not establish target runtime readiness.

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

For completed native/CH results with failed additional neo prediction, use
the same launcher with `--entry neo-predict --run-dir ORIGINAL_RUN --workflow-dir
PINNED_WORKFLOW`. `--entry neo-review --existing-predictions DIR` starts after
model inference. See [entry-point commands](releases/native-intake-v1.3.0.md).
Original packages/exports are retained and recovery results are separate.

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

## Slurm dependencies and two case chains

Append `--dependency afterany:PREVIOUS_HEAD_JOB_ID` to a prepared case’s printed
submission command. The launcher passes this verbatim as one sbatch argument;
Slurm validates job IDs and the dependency expression. It immediately prints the
new queued head ID; it does not wait on a login node or submit a polling job.
The selected expression is recorded in `state/head_dependency.txt`.

Use `afterany` for independent cases: the next head becomes eligible when the
previous head terminates, including failure/cancellation. `afterok` requires
success and can leave successors blocked if the predecessor fails. Native
compound expressions are accepted unchanged. CLI overrides an inherited
SBATCH_DEPENDENCY; the variable is then removed before head/nested submission.
Dependencies are per-submission, never stored as a shared site default. A
scheduler rejection propagates without recording a new successful submission.

To queue two independent chains, start one behind each of two existing heads.
Keep each chain’s most recently submitted head ID and make its next case depend
on that ID. Under normal end-to-end completion this permits at most two case
controllers, including native analysis and selected optional branches. It does
not limit the number of process jobs within each case. Scheduler availability
can leave fewer than two running. Static chains cannot share an idle slot when
one chain finishes earlier; they are not a dynamic two-slot work pool. A killed
or timed-out controller can leave child jobs; those need an explicit review
before treating that case as finished. No automatic cancellation/retry is added.
Running heads and their allocations remain unchanged.

Example (after preparing the next case):

```bash
bash "$mohcc/scripts/run_native_workflow.sh" \
  --manifest /prepared/samplesheet.csv --config /prepared/site.env \
  --dependency afterany:PREVIOUS_HEAD_JOB_ID
```

On Narval the maintained script is `run_native_workflow.sh` with the same options.

## Corrected phased export from completed native results

Use the same launcher with `--entry neo-prepare --run-dir ORIGINAL_RUN
--workflow-dir PINNED_WORKFLOW` to build new derivative neo inputs. Once that
entry completes, use `--entry neo-predict --existing-export NEW_EXPORT_DIRECTORY`
with the same original run and pinned workflow. Export-only writes
`entry_status.json` with prediction `NOT_RUN`; it does not declare model completion.
No native or CH task is in either recovery graph. See the
[phase interface](neoantigen_interface.md#5-phasing--sage-phase-sets--germline).
Do not point an active attempt to new source or overwrite its exports.
