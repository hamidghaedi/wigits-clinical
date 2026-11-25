# Laboratory QC on completed-BAM events

Owner correction November24,2026: replace five-minute QC submission with a
completed-BAM event. v0.3.0 adds only an operational native output subscription;
native analytical modules and current v0.1.2 QC jobs remain unchanged.

Offline verification:122 tests pass,3 existing skips;3 SAGE-compat adapter tests
pass. Event tests show that the last expected BAM starts QC even without final
caller results, published copies or final trace/version files. Earlier events
wait for the remaining libraries; repeated events do not submit another job.
A failed native task cannot declare readiness. QC success triggers delivery once.
Isolated source preparation preserves the original source and reuses its binding.
Shell syntax and diff checks pass. No environment, package or plugin changes.

A small Slurm notification canary is provided by
`workflow/test_lab_qc_events.sbatch`. It checks the prepared native source's help
entry and runs a synthetic completed-output subscription alongside a20-second
caller. The sink records receipt time; acceptance requires that notification
arrives before that caller finishes. It performs no patient analysis or metric
computation. Actual execution evidence is recorded after the canary completes.

Direct adversarial review: no timer invokes the dispatcher; CLI requires an
explicit BAM or completion event. Native output topics carry successful task
paths, avoiding asynchronous publishing races. A blocking state lock preserves
concurrent final notifications. Whole-run laboratory scope, original sources,
scientific holds and one-time submission intent remain intact. Operational
callback failures log warnings without reclassifying analytical success; missing
source bindings remain visible errors. Full-size QC performance remains under
measurement, separate from notification correctness. Manual/remote/RNA-only
routes without this registered source mapping are not claimed as covered.


## Slurm and activation result

Slurm12510626 completed exit0 in47seconds on November24 local time.
Native source help with its real reference configuration exited0. The actual
Nextflow topic notification arrived at04:01:33.677UTC November25, before the
synthetic caller finished at04:01:53UTC. Both source parsing and event timing
passed. Candidate12510595 failed because the help invocation omitted the native
genome configuration;12510599 failed on the synthetic output's literal metadata.
The harness corrections preserve these failed receipts and change no analytical
code. Accepted source is the published lab-qc/v0.3.0 tag.

Real-task metadata replay:14 notifications produced exactly2 captured run-level
submission commands,0 actual production jobs. Current v0.1.2 QC and its delivery
job remain unchanged. Coordinator ops/lab-qc-events/v1.0.0 /97ccb23 removes the
periodic dispatcher and freezes the output hook/config in new native packages.
Private configuration switched to the accepted event runtime at04:02:46UTC;
source/base/hook bindings and validation receipts remain outside Git. No plugin,
new environment, main merge or running-job modification. Ordinary report copying
continues on the existing archive schedule and never triggers BAM readiness.
