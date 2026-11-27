# Running a patient on Fir

For a new installation, use the [deployment requirements](DEPLOYMENT_REQUIREMENTS.md)
once to record host, container and selected-feature readiness. The current Fir
native/CH successes and unresolved additional neo failures are distinguished there.

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

New installations use shared `config/native_site.json`; the preparer falls back
to `config/native_fir.json` on existing installations. Its `laboratory_sources`
points to the original lab directories and their phase bindings; the maintained
alias registry contains source-confirmed corrections. New original lab files
go in those source directories, so later preparation discovers them automatically.
Routine native analysis uses fresh-frozen tumours; new FFPE analysis still needs
its separately agreed analytical route. No software or resources are installed
by either command.

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
bash "$mohcc/scripts/run_native_fir.sh" \
  --manifest /prepared/samplesheet.csv --config /prepared/site.env \
  --dependency afterany:PREVIOUS_HEAD_JOB_ID
```

On Narval the maintained script is `run_native_workflow.sh` with the same options.
