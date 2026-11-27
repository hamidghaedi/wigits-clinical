# native-intake/v1.2.0 — one launcher, explicit entry points

27 November 2026. This extends the existing `run_native_fir.sh` launcher rather
than adding case-specific helpers. Narval uses the same source under its existing
`run_native_workflow.sh` name. The sorting/isolation fixes from v1.1.5 are retained;
its tag and installed source remain unchanged.

| `--entry` | Inputs | Work started |
|---|---|---|
| `native` (default) | Prepared samplesheet and site config from the original-laboratory preparer | Existing native workflow and its selected extensions; existing `--resume` behavior unchanged |
| `neo-predict` | Original native run root and a pinned workflow directory | Reuse completed export; VEP, sorting/indexing, models and downstream review. No native/CH/export task |
| `neo-review` | Same, plus explicit completed prediction directory | Aggregation and candidate/native-NEO comparison only. No VEP, prediction, native/CH/export task |

Normal new analysis remains the preparer's printed command. For optional recovery:

```bash
launcher=/path/to/project/mohcc/scripts/run_native_fir.sh
workflow=/path/to/project/mohcc/releases/WiGiTS/native-intake-v1.2.0/workflow
run=/path/to/scratch/mohcc/analysis/WiGiTS/ORIGINAL_VERSION/ANALYSIS_ID

bash "$launcher" --entry neo-predict --run-dir "$run" --workflow-dir "$workflow"

# If models already completed and only downstream review needs repeating:
bash "$launcher" --entry neo-review --run-dir "$run" --workflow-dir "$workflow" \
  --existing-predictions /path/to/completed/rna_neoantigens
```

Both optional entries reuse the original package's site/options/export and
write to `extensions/recovery/native-intake-v1.2.0/<entry>/`, preserving original
native/CH/source products. They use a shared internal Slurm head, not another
operator launcher. `--workflow-dir` explicitly selects the repaired immutable
code without altering the original package's source binding. No entry is inferred
from filenames, and arbitrary internal processes are not exposed as entry points.

`--controller-time` and native `--dependency` syntax work for all entries.
Optional heads default to24h; full native heads retain96h/site/CLI behavior.
Submission dependency variables are cleared before child jobs. Entry/head/time/
dependency receipts are saved separately in the chosen recovery directory.
Repeated submission requires `--resume` after the previous head ends; no automatic
retry, cancellation, pool or new case discovery is added.

Existing v1.1.5 prediction-only sbatch commands remain valid against that pinned
release. New operations use this common launcher. The separate phasing limitation
recorded in v1.1.5 remains; no scientific thresholds or models change here.

Verification: [entry and portability tests](../../tests/results/optional-portability/README.md).
200 unit tests (3 existing skips),54 feature tests and3 adapter tests pass.
Real review-only execution through the launcher completes0:0 and retains63/63
candidate rows exactly, without repeating inference. The underlying model/export
implementation is unchanged from the full v1.1.5 runtime tests.
