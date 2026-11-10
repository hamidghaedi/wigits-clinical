# DNA adaptation branch migration

Canonical development location: **hamidghaedi/wigits-clinical,
`hartwig_adoptation`**. Local source:
`/global/project/hpcg6049/somatic_pipeline/wigits-clinical`.

The owner requested migration of the adaptation work from the similarly named
DNA branch before retiring that branch. The source's eight November9 commits,
`9d883c6` through `b3440c9`, change only11 files in
`test_data/hartwig_adaptation/`. They were imported as eight commits, retaining
their original authors and technical messages. The destination's no-attribution-
trailer convention applies. The exact source-to-imported commit map, file blobs
and modes are in [the migration record](DNA_HARTWIG_ADOPTATION_20261110.json).

All11 files and their modes were identical at imported commit`5cc854d`.
The only subsequent migration addition to that imported directory is a README
notice identifying its new repository. Runtime paths and historical source
commit references in the handover/evidence remain unchanged.

This transfers chromosome-test FASTQ preparation, the prepared native input
sheet, NEO/ORANGE metadata/terminology probes, extraction execution evidence
and the native39-task handover. It does not import the older DNA pipeline,
its historical reporting implementation or unrelated changes inherited by
that branch. The shared source base`54e2232` is retained by other existing
DNA branches/remote refs, including`origin/feat/purple-qc-routing-20261018`.

The newer native BAM launch, original full FASTQ launch, resource sizing and
efficiency reporting already live on this destination branch. Earlier
`codex/*` work branches are fully contained in its history. Existing numbered
release tags and frozen run packages preserve their original commits.

## Retirement boundary

Slurm migration test12444964 completed exit0 in7s: all22 existing tests and
the migrated synthetic chromosome-extraction test passed (23 total). This
used tiny synthetic reads only; the real test datasets and active full runs
were not regenerated, relocated or interrupted.

The source DNA branch and checkout have not been deleted; the owner will
perform that step. After the imported branch is pushed and migration checks
pass, the adaptation-specific work no longer depends on retaining that branch
name. Existing scratch inputs/results and frozen packages must remain.
Deleting a branch is different from deleting its checkout or runtime directories.

If new commits arrive on the source after`b3440c9`, migrate those separately
before deleting it. No source-side files or running jobs were changed here.
