# Staged native input release

Release `input-resources/v0.4.0`, November18,2026.

Owner relaxed the coordinator's independent demultiplexing-count prerequisite
for new native Oncoanalyser launches. All laboratory-bound lanes and both mates
remain required. Native tool inputs and analytical code do not change.
`workflow/conf/frontenac_paths.env` now supports `WIGITS_RUN_ROOT`, retaining its
historical default when unset. Resource coefficients and the full native entry
are unchanged. The coordinator owns staged-input queue preparation.

Validation with the existing Anaconda Python (PYTHONPATH unset):

- `python -m unittest discover -s tests`:75tests,3existing integration skips,
  exit0. This includes the explicit-root/historical-default regression.
- `bash -n workflow/run_full_native.sbatch workflow/run_oncoanalyser_frontenac.sh`:
  exit0.
- `git diff --check`:exit0.

An initial invocation using `python -S` could not import the existing PyYAML
dependency; rerunning without `-S` used installed packages and passed. No install
or environment modification was made. Real-case launch/accounting evidence stays
in the coordinator's private `phase2_native_20261118` queue and run packages.
