# Imported BAM preparation: measured costs and optimization assessment

**Superseded operational recommendation:** the owner subsequently directed
removal of both audits, not adoption of the buffering candidate. BAM entry
v0.3.0 follows the [current contract](../../../docs/BAM_INTAKE_CONTRACT.md).
The measurements below remain historical evidence only.

Date: November10,2026; stage times below are America/Toronto (EST).
Production intake source: `bam-entry/v0.2.0`, `3d06abd`.
Native Oncoanalyser source:3.0.0, `7c74c87a43749952b38c9a18915947570f0595a0`.
This is an investigation, not a change to running jobs or intake acceptance.

## Normal BAM: measured full-input stages

Source:172,246,996,661-byte BAM,1,894,142,243 records. Intake job12441709
started November10 at00:25:26 with64 CPUs/192GiB. The normal input is not yet
accepted: its prepared output is still undergoing the final complete audit.

| Stage | Evidence interval | Duration | Interpretation |
|---|---|---|---|
| Startup/reference/header/probe | Receipt00:25:29 to probe completion00:25:52 | About23s | Not the material cost |
| Complete source audit plus two checksum streams | Probe completion00:25:52 to source audit completion04:22:14 | 3h56m22s | Largest completed individual stage; verification added by this adapter |
| Name sort startup/merge before Picard starts traversal | Picard startup04:22:30 to first traversal04:50:21 | About28m | Sorting and tool startup overlap; not an isolated sort timer |
| Picard mate processing | Traversal04:50:21 to completion06:52:53 | 2h02m32s | Real mate-information preparation; limited parallelism |
| Coordinate sort | Picard completion06:52:53 to prepared BAM mtime07:36:49 | About44m | Approximation from output completion; command has no start timestamp |
| Fresh index and lightweight checks | Prepared BAM completion07:36:49 to verification start07:39:18 | About2m29s | Small compared with audits and repair |
| Complete prepared audit plus two checksum streams | Since07:39:18 | Running at checkpoint | Duration not yet known; do not mark PASS |

The actual sorting/repair/indexing interval was about3h17m, versus3h56m
for the initial audit alone. Final audit time is additional. Picard's log
reports150.38minutes from startup, including its wait for name-sorted input;
that is not an additional150minutes beyond the intervals above.

Live allocation-local process inspection at about07:51 found the final audit
process using99.8% of one CPU, while the two checksum processes used about39%
and14% and the decoder about16%. These are process-lifetime averages at that
checkpoint, not a controlled CPU profile. They establish that allocating64
cores does not parallelize the audit's serial SAM parsing. Both source and
prepared audits convert BAM to SAM, split/copy fields in a custom counter,
then parse the text again in two standard samtools checksum consumers.

## Bounded performance experiments

External diagnostic workspace:
`/global/scratch/hpc6049/wigits_runs/bam_intake_performance_20261110/r1/`.
All real read data remain there. The source BAM, active intake outputs and
frozen production package are read-only to these experiments.

- Slurm12444254:8 CPUs/16GiB/20min; `run.sbatch` invokes `benchmark.py`.
  Reads the first5,000,000 normal-BAM records into a small separate BAM.
  Two sequential repeats compare the released full audit with native binary
  all-record and mapped-record checksums and flagstat. The second repeat
  reverses execution order. Requires both checksum tables to match exactly.
- Slurm12444256:8 CPUs/16GiB/10min, afterok12444254;
  `buffered.sbatch` invokes `buffered_probe.py`. A scratch-only copy of the
  audit adds `ios::sync_with_stdio(false)` and `cin.tie(nullptr)`; no checking
  rules change. Two repeats must reproduce every audit counter and both
  checksum tables from the released implementation.
- A tiny synthetic fixture also tests whether installed samtools1.22.1
  `fixmate -p -z off` supplies the required supplementary MC tags. This is
  a compatibility test, not a real-data mate-preparation speed comparison.

| Experiment,5million records | Repeat1 | Repeat2 | Result |
|---|---:|---:|---|
| Released complete audit | 48.82s | 48.76s | Baseline |
| Native binary all-record checksum | 2.92s | 2.89s | Identical all-record checksum table |
| Native binary mapped-record checksum | 2.95s | 2.91s | Identical mapped-record checksum table |
| Native flagstat | 0.51s | 0.51s | Counts/flags only; not full compatibility audit |
| Same complete audit, C++ input buffering enabled | 8.86s | 8.89s | All counters and both checksum tables identical |

The smallest candidate changes only two C++ initialization statements:
`ios::sync_with_stdio(false); cin.tie(nullptr);`. The released program uses
default synchronized C++ input while serializing the same SAM records through
multiple consumers. Removing that input overhead sped up this complete audit
by **5.50times** (mean48.79s versus8.87s), without dropping any checks.
That is a measured subset result, not a whole-input or end-to-end speedup.
The source subset was254,554,866bytes; it samples the beginning of the
coordinate-ordered normal BAM and is not a representative random genomic
sample. Warm-cache and node differences also limit full-input extrapolation.

Slurm12444254 completed exit0 in2m05s;12444256 completed exit0 in20s.
An additional8CPU/32GiB regression job12444257 completed exit0 in1m28s:
**all24 intake tests passed** (86.150s test-suite runtime). It used the same
frozen test/source tree with only this scratch C++ buffering change. The tests
retain deliberate read-loss/sequence/quality mutation rejection and failed
intake launch gating. No rule was removed to improve timing.

Native checksum timing alone is not an intake replacement: it does not
implement every record-level compatibility check. Keeping the identical
complete audit with improved buffering is the smaller first change.

Installed samtools1.22.1 `fixmate` retained all12 synthetic records but left
the one eligible supplementary record without MC. It therefore fails the
current native boundary and is not a drop-in Picard replacement. Coordinate
errors in that test are expected because its output is name-sorted; the
missing supplementary MC, not order, is the relevant incompatibility. Do not
generalize this result to later samtools versions without separate testing.

## Decision boundaries and candidate improvements

1. Release/test the buffering optimization on a full-size input before replacing
   the documented Picard repair. Preserve
   all compatibility and read/quality conservation checks; avoid repeated SAM
   serialization/parsing where possible. A native binary checksum route is
   preferable to inventing a new checksum algorithm, but needs equivalent
   compatibility counters and testing before release.
2. Collect source checks during an already required preparation read, and
   output checks during writing where feasible. That can remove extra complete
   passes, but is not implemented or runtime-validated here.
3. The current `prepare()` unconditionally repairs and re-sorts accepted source
   BAMs. Introduce a separately tested no-repair route for inputs that already
   satisfy mate tags, order and index requirements; never assume MC is correct
   merely because the tag exists. This would align implementation with the
   contract's compatible-input branch. KQ0205 still needs repair because its
   eligible reads lack MC.
4. Future producer BAMs should emit consumer-compatible mate information once.
   That can avoid the name-sort/repair/coordinate-sort cycle for future imports.
   It is a proposed producer-boundary change, not a claim about current exports.
5. Retain BAM reuse as an option pending measured optimized intake and native
   REDUX acceptance. There is no full native FASTQ alignment timing for the
   same complete case here, so neither route is proven faster end to end.

No running job was cancelled, changed or restarted by this assessment. Do not
edit the active v0.2.0 package or bypass the existing PASS gate.

## Evidence and sources

- Full execution evidence: [EXECUTION.md](EXECUTION.md).
- Production receipt, `commands.log`, source audit and prepared BAM/index
  mtimes under the normal sample's directory in
  `KQ0205_full_bam_custom_20261109_r1/inputs/prepared_bams_v0.2.0/`.
- Adapter: `workflow/bam_intake.py`, `scan()` and `prepare()`;
  `workflow/bam_intake_audit.cpp`, the SAM stream reader.
- Native requirements: pinned `docs/usage.md`, "Sequencing technologies",
  recommends Picard FixMateInformation for missing mate-CIGAR tags.
- [samtools1.22 checksum](https://www.htslib.org/doc/1.22/samtools-checksum.html):
  supports BAM input, flags/tags/position checks and decoding threads.
- [samtools1.22.1 mate source](https://github.com/samtools/samtools/blob/1.22.1/bam_mate.c):
  `bam_mating_core`, primary versus supplementary processing.
