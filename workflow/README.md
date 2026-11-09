# workflow

Site configuration for Frontenac is in `conf/` (`frontenac.config`, `frontenac_limits.config`, `frontenac_paths.env`); the smoke-test launcher is `run_oncoanalyser_test_frontenac.sh`. All HPC paths come from `conf/frontenac_paths.env`, the machine-readable form of [`docs/directory_layout.md`](../docs/directory_layout.md).

Intended layout once work starts.

```
workflow/
  main.nf              thin wrapper that runs oncoanalyser then the clinical layers
  nextflow.config      profiles for the HPC side
  modules/
    serve/             build a SERVE-compatible actionability directory
    protect/           evidence matching against the ORANGE JSON
    rose/              conclusion sentences from the rules table
    reporter/          PDF and JSON report, including QC-fail variants
  assets/
    rules/             ROSE rules table, germline reporting gene table
    templates/         report templates
```

The base pipeline is `nf-core/oncoanalyser` in `wgts` mode. This wrapper does not modify it; it consumes the ORANGE JSON and the per-tool directories that oncoanalyser writes.
