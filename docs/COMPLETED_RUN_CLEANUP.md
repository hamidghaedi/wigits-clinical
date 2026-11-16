# Completed-run bulk cleanup

Owner direction of November 16, 2026 authorizes releasing large local files
before retrieving and processing the next glioma case from Fir nearline.

`workflow/cleanup_completed_run.py`, version 0.1.0, requires the exact completed
head job, run directory, source input plan, durable release directory and cleanup
plan path. It records the explicit eligible paths and file identities, then
`--apply` removes those regular files with a per-file journal.

Published results stay in the original run. Each source library must retain an
independent final BAM/CRAM there. Non-alignment results are also copied to the
durable release directory. Manifests, source metadata, native counts, commands,
logs, final alignments and their indexes remain available.

Eligible files are the source plan's restored FASTQs and local archive copies,
whose source nearline custody and successful restoration are already recorded,
plus regular intermediate BAM/CRAM/FASTQ files of at least 1GiB under this run's
`work/`. Directory symlinks are never traversed. Results and shared references
are outside this deletion scope. Content is not rehashed or re-read for cleanup.

Every sequential case finishes native analysis and results preservation, then
bulk cleanup, before the next case's payload transfer starts. A failed run
retains its files for recovery. The controller records pending and completed
states separately.

The initial KQ0004 cleanup was Slurm 12476557. Its five focused synthetic tests
cover terminal-state handling, exact scope, symlinks, final-library preservation,
retry behavior and changed-file detection. Detailed runtime/source receipts
remain with the coordinator's sequential glioma execution package.
