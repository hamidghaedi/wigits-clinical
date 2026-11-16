# ORANGE source-value mapping

Owner direction of November 16, 2026 permits correcting missing TINC handling
and deletion RANGE values while preserving finalized report prose.

The maintained component is `workflow/overlays/orange_source_values/`, version
0.1.0. It uses the installed ORANGE 5.0.1 algorithms, data model and PDF writer.
Its only mutation of the native ORANGE record is somatic DEL `geneRange`.
Chapter implementations and report text come from the original ORANGE JAR.
SomaReport's finalized templates retain their existing approved text and layout.

## TINC

`##tincLevel=` in the PURPLE somatic VCF establishes an available TINC estimate.
A measured zero remains available with value zero. An absent header produces
`status=unavailable`, `value=null`. Invalid or conflicting estimates also produce
an unavailable value, with the reason recorded. Comparison to QC allows its
four-decimal rounding.

PURPLE's original `TincLevel` string remains in `native_qc_value` for source
inspection. It must not populate a measured TINC display when the header is
absent. The separate AMBER contamination metric and native QC codes retain
their source values. Missing TINC does not establish a negative result or clear
a QC warning. The source-aware observation is written in `*.source_values.json`.
ORANGE 5.0.1 has no numeric TINC report field. No TINC section is added.

SomaReport independently applies the same rule in `wigits_case.tinc_record`,
recorded inside `wigits_fit.json`. It uses its bound `purple_somatic_vcf` input.

## Deletion extent

Join each somatic DEL on exact gene, transcript and chromosome to
`*.purple.cnv.gene.tsv`. Intersect that row's start/end with
`*.purple.cnv.somatic.tsv`. Both tables use one-based inclusive coordinates.
Chromosome names permit the explicit `chr` prefix equivalence.

Use PURPLE 4.4's existing DEL criterion, copy number strictly below 0.5.
`FULL` requires low-copy segments covering the complete annotated gene interval.
Multiple adjacent segments can jointly provide full coverage. `PARTIAL` requires
complete observed coverage of the gene interval with both deleted and retained
bases. A missing/ambiguous transcript, coverage gap, overlapping segments or lack
of a qualifying deletion gives `N/A`, with an internal reason. No exons or LOH
are inferred. AMP and other driver types keep their native fields.

The native record before mapping is retained as `*.orange.source.json`. The
receipt records the gene/transcript, interval, supporting segments and old/new
range. Analytical PURPLE files retain their source values.

Pinned method sources are
[PURPLE 4.4 DeletionDrivers](https://github.com/hartwigmedical/hmftools/blob/purple-v4.4/purple/src/main/java/com/hartwig/hmftools/purple/drivers/DeletionDrivers.java)
and
[ORANGE 5.0.1 PurpleInterpreter](https://github.com/hartwigmedical/hmftools/blob/orange-v5.0.1/orange/src/main/java/com/hartwig/hmftools/orange/algo/purple/PurpleInterpreter.java).

## Released execution

`workflow/prepare_oncoanalyser_source_values.py <upstream-checkout> <new-destination>`
builds source version `3.0.0-qlmp.2.2` from upstream
`7c74c87a43749952b38c9a18915947570f0595a0`. The released local source is
`05b8ec3ac953248d1f025528f03fe647c6f8ce30`. Native analysis, report eligibility
and reference resources are retained. The earlier QC-failure eligibility
experiment is outside this source release.

New site packages select this release through `ONCO_SOURCE_DIR`; an explicit
`ONCO_SOURCE` override remains available. For completed cases use native
`--processes_manual orange` with all
completed source directories. Bind any saved native sampling date explicitly
when re-rendering. The original case run and reports remain available.

Unit coverage is `tests/SourceValuesTest.java`, compiled against the installed
ORANGE JAR under Slurm. The first compute-node test, job 12476532, passed all
19 assertions. Python repository tests passed 68 tests with three skips using
the existing Python runtime. The `-S` invocation omits installed PyYAML and
cannot run the unrelated YAML intake tests.

Real-source TINC observations are retained under
`/global/scratch/hpc6049/report_source_mapping_20261116/`. KQ0004 and the
FFPE DeepSomatic substitution lack the TINC header and now have unavailable
estimates. The native FFPE control retains its available estimate 0.025.
The KQ0004 report-only execution uses its complete previous analytical outputs.
Job 12476538 completed exit 0 in 1m45s. Its single ORANGE task 12476544
completed exit 0 in 10s. All native fields are identical to the original report
except the two CDKN2A transcript ranges and CDKN2B, corrected from PARTIAL to
FULL. PARD6G remains PARTIAL. PDF text comparison confirms the same three
substitutions, and page 2 was visually inspected. The report has 11 pages.

Review artifacts and comparison receipts are retained at
`/global/project/hpcg6049/somatic_pipeline/wigits_releases/KQ0004_source_values_20261116/`.
Finalized SomaReport template assets and patient prose have no edits.

The first launch 12476533 ended before computation because a supplemental config
path was relative to the launch directory. The next attempt 12476534 rendered
correct values but its versions-file terminator caused aggregation to fail.
The maintained source builder corrects that shell terminator in 2.2. Source
2.1 and its failed execution remain available as evidence.
