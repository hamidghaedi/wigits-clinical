# Native somatic-call review

`workflow/audit_native_somatic_calls.py` v0.1.0 is a read-only comparison of
existing SAGE, PURPLE and ORANGE outputs. Run streaming VCF work through
`workflow/audit_native_somatic_calls.sbatch`, one array task per case. Original
results, reports and running analyses remain untouched. This is diagnostic QC,
not a validated false-positive rate, caller replacement or reporting policy.

The explicit private JSON manifest is a list of objects with `case_id`,
`orange`, `sage`, `purple`, and source metadata such as preservation and run
directory. Paths are supplied by the operator, never discovered across storage.
Case JSON and internal compressed PASS-site tables go to a separate QC directory.

```bash
sbatch --account=<account> --array=0-<last-index> \
  workflow/audit_native_somatic_calls.sbatch <python> \
  <released-script> <manifest.json> <qc-output-directory>
```

Compare all-record/filter counts, exact PASS records, SNV/non-SNV composition,
native VAF bins, normal support, qualities, fragment/read strand distributions,
unique fragment-coordinate counts and native trinucleotide context. Track the
selected report alleles back to both source VCFs, preserving the `reported`
flag and native caller counts. ORANGE catalog rows are not necessarily all
reportable alleles; do not interpret a whole-genome PASS count as a driver count.

Diagnostic thresholds are documented in the module docstring and are exploratory
review counters, not rejection criteria. No patient label, report prose or VAF
selection threshold changes. Missing support fields remain missing; AF is not
silently replaced with ALT/DP. Multiallelic rows are counted as records and not
expanded into allegedly independent mutations. Strand fractions are not strand
counts. These summaries are caller-derived, not an independent BAM review.

High C>T/G>A alone does not prove FFPE damage, and high burden alone does not
prove artefact. Examine native filters, low-VAF context, normal evidence,
laboratory preservation, signatures and CN/contamination QC together. Report
native PURPLE TMB unchanged; the review's record counts use different
definitions and cannot replace the native denominator. Never treat a workflow
exit 0 or an ORANGE PDF as scientific QC clearance.

Method sources: the executed VCF FORMAT/FILTER descriptions; pinned native
Oncoanalyser source/command and its published HMF component versions;
[SAGE documentation](https://github.com/hartwigmedical/hmftools/blob/master/sage/README.md)
for read-context and quality filters; and
[PURPLE documentation](https://github.com/hartwigmedical/hmftools/blob/master/purple/README.md)
for native burden/QC definitions. Live documentation is explanatory: do not
substitute current-master defaults for the actual executed source settings.
