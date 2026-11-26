# Laboratory QC and RNA resource consolidation

November26,2026. Source baseline v0.1.0 / `81998b3`; input branch tips
`75b72b1` (laboratory QC) and `8504373` (RNA input sizing). Both histories are
retained through ordinary merges on `codex/consolidate-lab-qc-rna`.

## Checks on combined source

| Command/check | Result |
|---|---|
|`python -m unittest discover -s tests`|Exit0;136 tests, OK,3 existing integration skips.|
|`python -m unittest discover -s workflow/sage_compat/tests`|Exit0;3 tests, OK.|
|`bash -n` for every tracked shell/Slurm source|Exit0;31 files.|
|`git diff --check`|Exit0.|
|Main, lab and RNA source ancestry|All three branch tips are ancestors of the consolidated branch.|

The merge's executable source matches the previously tested clean combined
tree `876934accfd4b62e560568e2b25bd54918c3bbbf`; subsequent changes update
repository version and documentation. No clinical/research pipeline argument,
caller or upstream container was changed while resolving integration.

Existing component runtime evidence is preserved in
[laboratory QC](../lab-qc-parabricks/README.md),
[BAM-ready events](../lab-qc-events/README.md) and
[RNA resource sizing](../../../docs/RNA_RESOURCE_SIZING.md).
The full-size lab run rendered exact4/10-library reports with the source
component; RNA metadata sizing and the completed recovery remain component
evidence. These are not new executions of this combined repository version.

The skipped tests retain their real external-runtime prerequisites. This
consolidation did not submit a new analytical/Slurm job, transfer input data or
repoint a running package. Fir site bindings and a target runtime canary remain
separate preparation work before new Fir execution.
