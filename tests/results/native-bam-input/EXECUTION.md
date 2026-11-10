# KQ0205 full paired-DNA BAM execution

Owner-authorized submission on November9,2026 (America/New_York).

- Slurm controller: **12441614**, observed RUNNING on `frnt076`.
- Site source: `740b501`, branch `codex/native-bam-input`, release
  `bam-entry/v0.1.0`; resource policy `input-v2`.
- Scientific source: unmodified nf-core/oncoanalyser3.0.0,
  `7c74c87a43749952b38c9a18915947570f0595a0`.
- Run ID: `KQ0205_full_bam_custom_20261109_r1`.
- Scope: existing full tumour and matched blood-normal DNA BAMs and adjacent
  BAI files; custom GRCh38.d1.vd1 FASTA/FAI/dictionary, full native REDUX and
  downstream callers. The prepared DNA-only sheet is used unchanged. No
  BAM-to-FASTQ conversion, ORANGE patch, manual stage selection or RNA row.

All runtime details are external under:
`/global/scratch/hpc6049/wigits_runs/KQ0205_full_bam_custom_20261109_r1/`

The frozen read-only `package/` contains scripts/configuration, samplesheet and
`release.json` with source hashes, input sizes and native source SHA.
`state/execution.json` contains the exact sbatch command and submission time.
`state/input_resources.{config,json}` records measured BAM sizes and the
initial allocations. The head log is `logs/native_bam_head.12441614.log`.

| Process | CPUs | Memory (GiB) | Initial walltime (hours) |
|---|---:|---:|---:|
| REDUX | 64 | 208 | 36 |
| SAGE somatic | 64 | 320 | 60 |
| SAGE germline | 64 | 240 | 28 |
| ESVEE | 64 | 432 | 80 |
| AMBER | 64 | 208 | 20 |
| COBALT | 64 | 208 | 20 |
| VIRUSBreakend | 64 | 256 | 48 |

The BAM sizes are a declared workload proxy, not measured primary-read counts
or a validated runtime model. Reference-index creation uses native tool-specific
resources. Head allocation is4 CPUs/8GiB/13days; genomic tasks run as separate
Slurm jobs. The generated configuration is loaded last to override the former
32-core ceiling. Live scheduler admission previously accepted64 CPUs/432GiB/80h
on the eligible CPU nodes.

At this checkpoint, Nextflow26.04.4 and Apptainer1.3.5 had initialized with the
custom reference and the intended input/resource configuration. Analysis and
report generation are **pending**, not established by head startup. Preserve
the native QC/report gates and the complete trace/accounting for later review.

The separate original-FASTQ relay controller12441547 remains active; one of16
files had landed on Frontenac. Its state and payloads are unchanged by this
BAM execution. No unrelated running analysis or transfer was cancelled.

## Subsequent status check: missing mate-CIGAR tags

The actual tumour and normal REDUX jobs12441619 and12441618 each requested
**64 CPUs,208GiB,36hours**, as confirmed in `.command.run` and Slurm accounting.
Both failed after about5seconds with exit1 and the native message
`required mate CIGAR not set, exiting`. This is the input-format prerequisite
described in upstream `docs/usage.md`, not resource exhaustion. No SAGE/PURPLE
analysis has started. The requirement should have been checked before submission.

The head remains running while native reference-index job12441617 completes
(6 CPUs/36GiB/8h; observed about40minutes elapsed). Preserve this reference work.
The native documented repair is Picard `FixMateInformation` with mate-CIGAR
generation. Prepare new copies and retain the original BAMs before resuming;
do not remove the native requirement or invent MC tags from a read's own CIGAR.
Repair and analytical resumption remain pending at this checkpoint.

The independent FASTQ controller remains ACTIVE with no recorded exception:
2/16 files (71,819,382,338 bytes) verified on Frontenac. File3, tumour lane2 R1,
has11,777,081,344 of35,301,263,198 bytes transferred on hop2 with statusOK,
zero faults and no fatal error. Transfers remain unblocked.

## Systematic intake and automatic resume, November10

Release `bam-entry/v0.2.0`, source `3d06abd7d0ec9c95f855c20bf885c9e88d14a231`,
implements the mandatory [intake contract](../../../docs/BAM_INTAKE_CONTRACT.md).
Slurm test12441700 completed exit0 with24 tests passing using actual
samtools1.22.1 and GATK4.6.1.0, including deliberate record/sequence/quality/
annotation/alignment mutations, raw-input bypass refusal and failure gating.
Development-only trials caught and corrected binary stdout contamination by
the GATK wrapper before any real preparation; GATK now writes a normal BAM file.

The frozen package is under the same execution root at
`package_intake_v0.2.0/`, with source hashes and submission provenance in
`release.json`. Original `package/`, inputs, work/cache and logs remain unchanged.

| Job | Role | Allocation | State observed |
|---|---|---|---|
|12441707|Intake/resume controller|4 CPUs,8GiB,13days|RUNNING|
|12441708|Full tumour DNA intake|64 CPUs,192GiB,88h initial limit|RUNNING|
|12441709|Full matched-normal DNA intake|64 CPUs,192GiB,64h initial limit|RUNNING|

These are conservative input-volume limits, not runtime predictions. Every
record is inspected and conservation-checked; new prepared BAMs/indexes will
be under `inputs/prepared_bams_v0.2.0/<sample>/`. The source BAMs are not edited.
Receipts and dispatch state are stored beside the preparation outputs;
`state/intake_resume_execution.json` names the controller and exact command.

After both intake jobs pass and their prepared-input/reference bindings are
verified, the controller waits for the prior controller12441614 to become
terminal, then invokes the unchanged native pipeline with `-resume` in the
existing execution root. This preserves and reuses successful reference tasks
without simultaneous access to the same Nextflow cache. REDUX/callers have not
resumed yet at this checkpoint; they are gated on intake completion.

The independent FASTQ controller12441547 remains ACTIVE, with3/16 original
DNA FASTQs landed. No transfer job or reference job was cancelled.

## Follow-up: full scans active, reference cached, no new blocker

After about50minutes, both intake jobs12441708/12441709 remain RUNNING with
64CPUs/192GiB allocated. Their logs are at full source inspection/conservation
checksumming; neither has published PASS or begun the native resume. Live Slurm
accounting and an allocation-local process/I/O probe confirmed active processing.
The normal decoder had consumed about36.9GB of its172.25GB source; the current
checksums/audit files are finalized at end-of-stream, so empty intermediate
checksum files and the INSPECTING receipt are not themselves a stall.

Reference-index task12441617 completed exit0 in1h28m26s. Its successful work
directory/hash82/d8ffef and original cache remain in place; the original head
12441614 then ended FAILED because of the earlier REDUX error. The new controller
12441707 remains waiting for both verified intake outputs, with no new error.

FASTQ relay12441547 remains ACTIVE:5/16 files (178,776,539,315 bytes) have landed.
File6, tumour lane3 R2, is transferring on hop2 with8,514,699,264 of36,534,722,254
bytes copied, statusOK, zero faults and no fatal error. No new blocker or source
modification was found; no active job was restarted or cancelled.

## Subsequent checkpoint: FASTQ landing complete; intake progressing

FASTQ controller12441547 completed exit0 in5h45m33s. All16 original DNA FASTQs,
437,332,334,327 bytes total, are verified locally through the recorded two-hop
transfers. No active transfer remains and no exception was recorded.

At7h16m of intake runtime, both workers remain RUNNING with64CPUs/192GiB.
The normal source scan inspected1,894,142,243 records with no reference/RG,
supplementary-clipping, missing-sequence or missing-quality blocker. Missing MC
was confirmed on1,892,098,646 eligible records. Picard preparation and indexing
have finished; the new97,387,227,810-byte BAM and9,436,680-byte BAI are undergoing
full post-preparation conservation checks. This is **not yet a PASS**.

The tumour remains in full source inspection; live process I/O showed about
212.5GB of its242.6GB source consumed, with sustained CPU/read activity.
Neither intake receipt has passed, so controller12441707 has correctly not
resumed native analysis. Original BAMs and completed reference work remain
untouched; no new blocker or resubmission was identified.
