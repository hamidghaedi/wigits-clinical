# native-intake/v1.3.0 — supported tumour phase for pVACseq

27 November 2026. Adds a derivative SAGE-LPS adapter and an export-only entry to
the existing launcher. Main candidate selection, measured DNA/RNA/expression,
prediction models, thresholds, native NEO and CH are unchanged. Old packages and
running jobs retain their source; this release does not update KQ0061's run.

All selected somatic anchors appear in proximal context. Only unambiguous,
positively supported tumour somatic LPS becomes FORMAT/HP. Germline LPS IDs do
not establish somatic phase. Ambiguous/absent evidence stays unresolved; actual
tumour homozygotes use pVACtools' existing rule. Source fields, support, omitted
context and limitations are retained in a per-allele JSON ledger. Actual converter
uptake is recorded separately after prediction. See the [interface](../neoantigen_interface.md#5-phasing--sage-phase-sets--germline).

Use the shared launcher, with paths from your actual source/attempt:

```bash
launcher=/path/to/project/mohcc/scripts/run_native_fir.sh
workflow=/path/to/project/mohcc/releases/WiGiTS/native-intake-v1.3.0/workflow
run=/path/to/scratch/mohcc/analysis/WiGiTS/ORIGINAL_VERSION/ANALYSIS_ID
case=PATIENT_ID

bash "$launcher" --entry neo-prepare --run-dir "$run" --workflow-dir "$workflow"

# After neo-prepare completes successfully:
export_dir="$run/extensions/recovery/native-intake-v1.3.0/neo-prepare/results/$case/neo_export"
bash "$launcher" --entry neo-predict --run-dir "$run" --workflow-dir "$workflow" \
  --existing-export "$export_dir"
```

`native` and `neo-review` retain their previous meanings. `--existing-export`
also works for review, while prediction without that flag still selects the
original export explicitly. Optional entries write isolated versioned attempt
directories, default24h, support native Slurm `--dependency` and configurable
`--controller-time`, and never invoke native alignment/calling or CH. Repeating
an ended attempt requires `--resume`. Export-only status means prepared inputs,
not completed prediction. Reusing old prediction scores after changing peptide
context is scientifically invalid; predict the new sequences first.

Runtime validation is recorded in `tests/results/neo-phasing/README.md` when
complete. The source release must not be interpreted as scientific sign-out or
as support for pVACtools' unsupported non-missense/ambiguous proximal contexts.
