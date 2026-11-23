# Existing native run: verified completion at owner handover

Checked 9 November 2026, America/New_York. The owner asked to take over
submission/monitoring from Cloud Science. Inspection found that the native run
had already been submitted and was finishing; no duplicate was submitted.
No separate handoff comment was found in the shared coordination folders;
its location has been requested. This record documents observed execution,
not an unverified agent shutdown or a read of an unavailable chat comment.

## Execution

| Item | Verified value |
|---|---|
| Run ID | `KQ0205_chr17_19_20261109_onco3.0.0` |
| Controller job | **12441051**, COMPLETED, exit **0:0** |
| Start / end | 9 November 17:57:10 / 20:36:41 EST |
| Elapsed | 02:39:31 |
| Nextflow session | `91935fca-964a-47bb-9aa7-1b30ec403185` |
| Native pipeline | Oncoanalyser **3.0.0**, `7c74c87a43749952b38c9a18915947570f0595a0` |
| Mode / genome | `wgts` / `GRCh38_hmf` |
| Include/exclude overrides | None in captured parameters; native stage selection |
| Trace outcome | **39 COMPLETED, exit 0; zero failures, retries or cached tasks** |
| Samplesheet SHA-256 | `42922dba7dbd791b74e9cd6c5e267d47bcae071606948e4cdbe30a76b0e27f1b` |

The run's samplesheet is byte-identical to the prepared KQ0205 sheet under
`hartwig_test_fastqs_20261109/KQ0205/`. All three libraries were used: tumour
DNA, normal DNA and tumour RNA.

Captured launch uses the existing `wigits-clinical` launcher. That repository
was clean at inspection, HEAD `ae5175b3d0ae9f63eaa740b77e1a620bf933feaf`;
this current checkout observation is not a claim of a separately archived
wrapper commit at submission.

```bash
bash /global/project/hpcg6049/somatic_pipeline/wigits-clinical/workflow/run_oncoanalyser_frontenac.sh \
  KQ0205_chr17_19_20261109_onco3.0.0 \
  /global/project/hpcg6049/somatic_pipeline/wigits-clinical/samplesheets/KQ0205_chr17_19_20261109_onco3.0.0.csv
```

## Evidence and outputs

Runtime:
```text
/global/scratch/hpc6049/wigits_runs/KQ0205_chr17_19_20261109_onco3.0.0/
```

- `launch.sh`: captured native Nextflow command and launch metadata.
- `.nextflow.log`: terminal “Pipeline completed successfully”, execution
  statistics and session identity.
- `results/pipeline_info/execution_trace_2026-11-09_17-57-21.txt`: all 39 tasks.
- `results/pipeline_info/params_2026-11-09_17-57-34.json`: actual run parameters.
- `results/KQ0205_chr17_chr19/`: native analytical outputs.
- `results/multiqc_report.html`: generated MultiQC report.
- Controller stdout:
  `/global/scratch/hpc6049/.claude-science-scratch/.claude-science/jobs/36564270-59ce-4cd8-9f94-45218ad13b3b/slurm-12441051.out`.

The trace includes DNA/RNA alignment, REDUX, ISOFOX, AMBER/COBALT, SAGE,
PAVE, ESVEE, PURPLE, LINX, CHORD, signature fitting, CUPPA, NEO and MultiQC.
This establishes execution on the deliberately selected chromosome dataset;
it is not whole-genome scientific or patient-report clearance.

## ORANGE did not run

No ORANGE task appears in the completed trace and no ORANGE PDF is published.
The native `orange_reporting` runnable condition requires a PURPLE circos plot;
`purple/plot/KQ0205-DNA.circos.png` is absent, although the PURPLE somatic VCF is
present. Therefore successful workflow exit must not be described as completed
ORANGE report generation. No native gate or analytical result was modified to
force a report from this chromosome-only test.

The existing whole-genome ESVEE comparison is separate and was preserved.
