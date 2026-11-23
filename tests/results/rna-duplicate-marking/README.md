# RNA optical-comparison regression

21 November 2026, existing Frontenac GATK 4.6.1.0 container; no installation.

- Slurm 12506204: COMPLETED, exit 0:0, elapsed 24 seconds, 1 CPU / 4 GiB.
- Baseline MarkDuplicates: 7.10 seconds; RNA setting: 5.86 seconds.
- 15 read pairs, 12 duplicate pairs, 80% duplication in both modes.
- First 11 SAM fields identical for all 30 records, including duplicate flags.
- Optical duplicate pairs: 12 in baseline; deliberately unassessed in RNA
  mode (the native auxiliary counter is written as 0).
- Unit suite: 83 tests pass, 3 existing integration skips. Bash syntax checks
  cover both edited launch/recovery scripts and the new test batch file.

Runtime files and real-case recovery inputs remain outside Git. The small
test proves the stated flag/interface invariants for its fixture; full-size
RNA recovery completion is a separate recorded execution outcome.
