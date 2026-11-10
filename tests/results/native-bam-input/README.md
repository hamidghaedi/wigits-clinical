# Native BAM entry development checks

November9,2026. Upstream nf-core/oncoanalyser3.0.0 at
`7c74c87a43749952b38c9a18915947570f0595a0`; Nextflow26.04.4,
Apptainer1.3.5. No analytical source modification or real BAM-analysis
submission in this increment.

| Check | Outcome |
|---|---|
| `test_full_input_preparation.py` | 6 tests pass; synthetic recovery plus count/bytes resource inputs, multi-lane aggregation, repeated-path refusal and envelope refusal |
| `test_frontenac_launcher.py` | 1 test/four configuration discovery paths pass |
| `test_two_hop_transfer.py` | 1 test passes; saved task adoption and hop ordering |
| `bash -n` native and BAM launchers | Exit0 |
| `git diff --check` | Exit0 |
| Native configuration flattening | Exit0; exact d1.vd1 paths; CPU ceiling64; REDUX64/208GiB/36h; somatic SAGE64/320GiB/60h; ESVEE64/432GiB/80h |
| Native BAM `-preview` | Samplesheet accepted and graph compiled; log printed success; trace contains header only, zero tasks. **Not a clean end-to-end pass:** Nextflow then waited in `OutputDsl.getOutput`. Preview JVM terminated by SIGTERM; enclosing launcher exited1. |

Invocation from the repository root:

```bash
WIGITS_WORKFLOW_DIR="$PWD/workflow" \
  bash workflow/run_bam_oncoanalyser_frontenac.sbatch \
  KQ0205_bam_entry_preview_20261109 \
  samplesheets/KQ0205_full_dna_bam_custom.csv -preview -ansi-log false
```

Private evidence: `/global/scratch/hpc6049/wigits_runs/KQ0205_bam_entry_preview_20261109/`,
including generated `state/input_resources.{json,config}`, `state/resolved.config`,
`state/preview_assessment.json`, native `.nextflow.log` and the zero-task trace.
The resolved FASTA is the original GRCh38.d1.vd1 reference; missing custom-genome
indexes are to be built by native reference preparation in a real Slurm run.

Remaining: actual BAM-format/reference compatibility checks, native missing-index
generation, full BAM execution and comparison with the preserved full native
FASTQ baseline. No production equivalence, completed comparison, report output
or caller/template change is claimed. Do not patch upstream to make preview
finalization pass.
