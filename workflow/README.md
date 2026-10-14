# workflow

Nothing runs from here yet.

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
