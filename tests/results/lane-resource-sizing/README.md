# Lane-local resource release validation

Release `input-resources/v0.5.0`, November19,2026. See
[policy and coefficients](../../../docs/LANE_RESOURCE_SIZING.md).

- `env -u PYTHONPATH python -m unittest discover -s tests`:77tests,
 3existing integration skips, exit0. No environment/package changes.
- `bash -n tests/run_lane_resource_canary.sbatch`:exit0.
- `git diff --check`:exit0.
- Metadata-only real Nextflow canary Slurm12487190:COMPLETED, exit0:0.
  Head4CPUs/8GiB, Nextflow26.04.4; native analytical tools never executed.

Canary command:

```bash
sbatch --account=def-hpcg6049_cpu --output=<runtime>/lane_canary_%j.log \
  tests/run_lane_resource_canary.sbatch <repository> <new-canary-directory>
```

The canary generated eight5GiB normal lane pairs and eight34GiB tumour lane
pairs as sparse synthetic files (no patient data, no FASTQ reading), generated
the actual resource config, then ran four echo-only local processes within the
Slurm head. Configured local-executor capacity is synthetic so task directive
resolution can be inspected without reserving/requesting those resources on
the cluster. Actual work uses only the small head allocation.

| Native task metadata | Resolved CPUs | Memory GiB | Time h |
|---|---:|---:|---:|
| Normal lane5GiB |16|96|4|
| Tumour lane34GiB |64|192|6|
| Complete normal sample40GiB |64|144|16|
| Complete tumour sample272GiB |64|224|44|

All four match the Python preview exactly. Synthetic sparse inputs were removed
by the canary after PASS; metadata, outputs and logs remain under the external
Phase2 queue's `canary_input_v3_r1/`. Warnings about absent paired-caller process
selectors are expected in this two-process canary. It does not validate these
memory/time limits universally or establish analytical equivalence/performance.

Live release currency check confirms native oncoanalyser3.0.0/HMF3.0.0--8
remain current. The native analytical source, original FASTQs, references,
caller settings and completed work remain unchanged. Cancellation/restart
receipts and actual new task allocations stay in the coordinator's private
`phase2_native_20261118` queue. No original execution package is overwritten.
