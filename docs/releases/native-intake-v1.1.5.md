# native-intake/v1.1.5 — optional neo portability fix

27 November 2026. Based on published v1.1.4 and the recorded Fir diagnosis.
Native callers, CH calling, references, containers, predictors and scientific
selection/ranking thresholds are unchanged. Source is published as an immutable
component on the maintained `codex/reusable-native-intake` branch; this is not a
main promotion or a replacement of already issued packages.

Changes:

- Sort task-owned selected/proximal VEP output by header contig order and
  numeric position before bgzip/tabix. Complete record lines, headers,
  duplicate/multiallelic records, sample values and phase fields are retained.
  Existing unsorted exports can be reused without editing them.
- Explicit container `--env PYTHONNOUSERSITE=1 --env PYTHONPATH=` excludes host
  Python user-site packages and overrides injected Python paths, retaining the
  image's own NumPy/pandas and numerical libraries.
- `PREDICT_EXISTING_NEO` starts with a completed export. It cannot invoke native
  analysis, BLOOD_CH or NEO_PREPARE. It checks the export's case/sample binding
  and writes prediction/review products into a separate recovery directory.
- Status collection reads completed CH status from the original results without
  copying or altering it. Current recovery failure takes precedence over prior
  completion. The original failed attempt remains evidence.

## Resume a failed optional prediction

Use the script from the newly installed pinned release, passing the **original
native run root**. This is the same command for any eligible case:

```bash
workflow=/path/to/project/mohcc/releases/WiGiTS/native-intake-v1.1.5/workflow
run=/path/to/scratch/mohcc/analysis/WiGiTS/ORIGINAL_VERSION/ANALYSIS_ID
sbatch --account=SITE_ACCOUNT \
  --output="$run/logs/neo-recovery-%j.out" \
  --error="$run/logs/neo-recovery-%j.err" \
  --export="ALL,WIGITS_RECOVERY_WORKFLOW=$workflow" \
  "$workflow/resume_optional_neo.sbatch" "$run"
```

The head defaults to4CPU/8GiB/24h, including child queue waits; override its
walltime with native `sbatch --time` if needed. Prediction retains the completed
export's input-adjusted task allowance. The head reuses the original site/module
binding, completed native receipt, original option/resource configuration and
`extensions/results/<case>/neo_export/`. New work/results/status go to
`extensions/recovery/native-intake-v1.1.5/`. It records the new workflow path;
the versioned installed source identifies its exact commit. Do not submit a
second recovery head while one is active. No automatic retry or cancellation.

Original native/ORANGE results, original CH results and export files stay in
place. The recovery status references existing CH status directly. Native-only
tracker/report completion does not imply additional prediction success; use
the recovery's explicit status for the additional feature. Historical failed
status files are preserved rather than overwritten.

Verification receipts are recorded in
[`tests/results/optional-portability/README.md`](../../tests/results/optional-portability/README.md).
