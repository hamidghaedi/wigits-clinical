# Run a prepared case on Fir

Use the manifests supplied in `mohcc/file_manifests/<analysis_id>/`.
The CSV lists every tumour and normal lane and both mates. `site.env` selects
the pinned WiGiTS release, native reference, account and output directory;
`blood_ch.json` binds the recorded blood normal and post-native CH analysis.
These are inputs to the launcher, not files to reconstruct by hand.

From a Fir login shell:

```bash
bash /home/ghaedi/projects/def-ghaedi/ghaedi/mohcc/scripts/run_native_fir.sh \
  --manifest /path/to/supplied/samplesheet.csv \
  --config /path/to/supplied/site.env \
  --blood-ch
```

It prints the Slurm head job ID and run directory. The head launches native
Oncoanalyser tasks, followed by the selected blood CH task. Omit `--blood-ch`
for a new native-only case. RNA-informed neoantigens remain independently
selected in the prepared configuration; no RNA is assumed from a case name.
The launcher never transfers FASTQs or copies reference/container payloads.

```bash
squeue -u "$USER"
tail -f /path/printed/by/the/launcher/logs/head-<job_id>.out
```

Native results: `<run>/results/<case>/`; native trace/report:
`<run>/results/pipeline_info/`; optional CH:
`<run>/extensions/results/<case>/blood_ch/`; separate optional status:
`<run>/extensions/state/status.json`. The native completion receipt is retained
if a later optional task fails. Resource CSV/HTML is under native
`results/pipeline_info/resource_usage/`; it observes jobs without changing them.

If the head fails, inspect its error and task logs. After addressing the cause,
repeat the **same** command with `--resume`. It reuses the original package and
Nextflow work/cache, and refuses a second head while the prior one is active.
Do not change the manifest or release inside an existing execution.

The module bindings default to Nextflow26.04.4 and Apptainer1.3.5, the existing
tested versions. These must be available on Fir; no installation or upgrade is
performed. The first owner-operated Fir run is the target-site runtime check.
DNA BWA stays64CPUs with lane-sized RAM/walltime. Fir's final production ceiling
is64CPUs/768GiB/7days. The head requests4CPUs/8GiB/6days23hours. These are
limits, not predicted durations. Laboratory QC retains its BAM-ready event
entry and site selection; full laboratory run membership and delivery need the
corresponding issued QC configuration, not a claim that one case is a run.
