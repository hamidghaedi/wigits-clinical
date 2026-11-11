# DeepSomatic whole-genome performance review — KQ0205 on Frontenac (run D vs route A)

Date: 2026-11-11 (snapshot of job 12447879 taken at 15 h 35 m elapsed, ~09:53 local). Read-only investigation: no job was submitted, modified or cancelled; no repository file was changed. All numbers below come from `sacct`/`sstat`/`scontrol`, the stage logs, BAM headers/indexes and three sampled 10 kb BAM windows (see Appendix A). Items marked **[provisional]** rest on partial data (run D is still in its first stage); items marked **[DeepVariant design — not profiled here]** are stated from DeepVariant/DeepSomatic documentation or source knowledge, not from a measurement on this cluster.

## 1. Executive summary

- **Run D is not slow because of its Slurm configuration; it is slow because its input BAM is ~3× more expensive for `make_examples` and yields ~2× more candidates.** The REDUX tumour BAM has **84–86 % of non-duplicate reads soft-clipped (mean clip ≈ 60 bp)** in three sampled windows, vs **8–14 %** in route A's Parabricks BAM. Cause: ~80 % of this FFPE tumour's fragments are shorter than the 151 bp read; route A aligned **fastp-trimmed** FASTQs (Parabricks `@PG`), the oncoanalyser native alignment aligned **raw** FASTQs (bwa-mem2 `@PG`; oncoanalyser 3.0.0 only runs fastp for UMI/split). bwa-mem2 soft-clips the read-through adapter, and soft clips drive DeepVariant's realignment-window selection. Measured consequence: **15.3 candidates per busy-core-second in D vs 52.8 in A (3.46× more CPU per candidate)**, and 53.97 M candidates after 15.6 h vs A's 33.7 M final.
- **Stage 1 of D is CPU-bound and healthy**: `sstat` AveCPU = 40 d 22 h over 15.58 h = **63.0 of 64 cores busy**, MaxRSS 109 GB of 256 GB requested, frnt128 = 64 physical Intel cores (ThreadsPerCore=1), CPULoad 69.7. It is not waiting on GPFS. **[provisional]** By cumulative read volume (21.34 TB read; each shard reads ≈ the whole 434 GB BAM pair, as route A's trace confirms: 12.5 TB ≈ 32 × 415 GB) the stage is 77–85 % done, projecting **≈ 18–20 h** for `make_examples` alone, then **≈ 10–11 h** of inference on an A30 (A's measured 1,522 examples/s) or **≈ 19–21 h** on an RTX6000 (700–780 ex/s measured locally), plus ~25 min postprocess: **≈ 29–32 h (A30) to 38–42 h (RTX6000) end-to-end, versus route A's 11 h 36 m**.
- **Run D's pending inference job (12447880) requests an untyped `--gres=gpu:1` with a 72 h limit**, which routes it to `gpubase_14days`, a partition that contains RTX6000 nodes as well as A30/A100 nodes. Landing on an RTX6000 would double the inference time. Requesting a typed GPU and ≤ 24 h is the cheapest single fix.
- **Doubling shards from 32 to 64 doubled the BAM bytes read from GPFS (27.8 TB vs 13.3 TB projected) for a sub-linear gain**: per-shard throughput fell from 20.5 (16 shards, pilot) to 15.5 candidates/s (64 shards). Shards = allocated CPUs is still right, but the better way to get more cores is **several region-split `make_examples` jobs on the CPU partition feeding one GPU job**, which also makes `$SLURM_TMPDIR` usable (D's example set is ≈ 800 GB–1 TB, larger than the 462 GB node-local disk, so a single-job run cannot use it).
- **Memory and postprocess are over-requested**: 1.7 GB/shard observed (A: 1.6 GB/shard), so 64 shards need ~160 GB, not 256 GB; route A's postprocess took 10.6 min at 32 CPUs/96 partitions (D requests 64 CPUs/256 GB/12 h). Route A's single-task design idled its A30 for 57 % of 11.6 h; D's three-stage split already fixes that and should be kept.

## 2. Configuration comparison

| Item | Route A (dnaseq_nf `DEEPSOMATIC_FFPE`, job 12414711) | Run D (wigits-clinical `deepsomatic_full.py` 0.2.0, jobs 12447879/80/81) | dnaseq_nf module as released (`modules/deepsomatic_ffpe.nf`, `conf/gpu.config`) |
|---|---|---|---|
| DeepSomatic | 1.10.0 GPU image (same `.sif`) | same | same |
| Model | FFPE_WGS, `--call_small_model_examples`, same `flags_for_calling` | same | same |
| Reference | GRCh38.d1.vd1 (GDC; 2,779 contigs, decoys + viral) | HMF `GRCh38_masked_exclusions_alts_hlas.fasta` (3,367 contigs; 57 tumour reads on alt contigs, 0 on HLA) | GRCh38.d1.vd1 |
| Tumour BAM | `KQ0205-DNA.analysis_ready.bam`, 242.6 GB, Parabricks fq2bam 4.6.0-1 from **fastp-trimmed** FASTQ, MarkDuplicates, BQSR table generated but qualities not recalibrated (binned 37/25/11) | `KQ0205-DNA.redux.bam`, 316.7 GB, bwa-mem2 `-Y -K 100000000` from **raw** FASTQ + REDUX; same flowcell HMNFNDSXF; binned qualities | n/a (module consumes Parabricks BAM) |
| Normal BAM | 172.2 GB, same provenance | 117.0 GB REDUX; flowcell HNMH2DSXF | n/a |
| Mapped / unmapped reads (idxstats) | T 4.587 G / 7.7 M; N 1.893 G / 1.4 M | T 4.731 G / **305.9 M** (REDUX unmapping); N 2.033 G / 48.7 M | — |
| Depth (BamMetrics / Parabricks metrics) | T 63.6× / N 74.8× | T 64.1× / N 77.5× (OverlappingReadPercent 37 % tumour) | — |
| Soft-clipped non-dup reads (3 windows) | T 8–14 % (mean clip 44–55 bp); N 3–6 % | **T 84–86 % (mean clip 59–60 bp)**; N 7–12 % | — |
| Duplicate-flag rate (chr2 window) | T 10.1 %, N 8.2 % | T 18.8 %, N 14.6 % | — |
| `--num_shards` | 32 (= task.cpus) | 64 | `${task.cpus}` = 32 |
| `--regions` | none (whole genome) | none; runner refuses any `--regions` | optional `deepsomatic_ffpe_intervals` BED |
| `--intermediate_results_dir` | `ds_tmp` in the Nextflow work dir (GPFS `/global/scratch`); 433 GB still on disk | `results/intermediate` (GPFS); 793 GB examples + 0.67 GB CVO so far | `ds_tmp` (GPFS) |
| `--sample_name_tumor/normal` | KQ0205-DNA / KQ0205-gDNA | same | same |
| `call_variants --batch_size` | 256 | 256 (added post-submission to `commands.json`, 2026-11-11 05:41 UTC; before that the native default 1,024 per replica) | 256 (`deepsomatic_call_batch_size`, canaried 128 vs 256 on 8,192 examples) |
| GPU | 1 × **A30** (frnt144, AMD 48-core node), untyped request | pending, untyped `--gres=gpu:1`, 16 CPUs, 128 GB, 72 h → `gpubase_14days` | untyped `--gres=gpu:1`, exclude list frnt110,143,153–156,190,191,201 |
| `postprocess --cpus/--num_partitions` | 32 / 96 | 64 / 192 | `${task.cpus}` / 3×cpus |
| CPU / mem / time requested | 32 CPU, 160 GB, 18 h, one task | examples 64 CPU/256 GB/60 h (CPU partition); inference 16 CPU/128 GB/72 h; postprocess 64 CPU/256 GB/12 h | 32 CPU, 120 GB (+40 GB/retry), 18 h × attempt |
| Candidates / examples | 33,705,090 / 27,149,056 (final) | 53,969,760 / ≈ 46 M at 15.6 h, still running | KQ0315 reference: 20.85 M examples |
| Stage wall times | examples 6 h 27 m; inference 4 h 57 m; postprocess 10 m 35 s; total 11 h 35 m 46 s | examples ≥ 15 h 35 m (running) | KQ0293 6.76 h / 7.85 h; KQ0294 (127×) 16.74 h examples on RTX6000 node |
| CPU efficiency | 1,732 % over the whole task; 27.5 cores busy during examples | 63.0 cores busy | — |
| Peak RSS | 60.7 GB (trace), 50.5 GB (sacct MaxRSS) | 109.0 GB (sstat, running) | — |

20 Mb pilot (job 12446769, chr17:1-10 M + chr19:1-10 M, 16 shards, RTX6000 node frnt153): 636,429 candidates in 38 m 13 s (user 564 min = 14.8 cores busy; 20.5 candidates/s/shard), inference 513,024 examples in 12 m 13 s at default batch 1,024 (700 ex/s), postprocess 1 m 53 s, MaxRSS 23.2 GB.

## 3. Where the time goes

### Route A (job 12414711, 32 CPUs + A30)
| Phase | Wall | CPU (user) | Notes |
|---|---:|---:|---|
| make_examples | 6 h 27 m (387.0 min) | 10,634.6 min = 177 core-h | 27.5/32 cores busy; 1,452 candidates/s; 42.8 s per 2,048-candidate chunk per shard (16,884 chunks) |
| call_variants | 4 h 57 m (297.2 min) | 908.9 min | 1,522 examples/s steady state (0.066 s/100 at end); GPU-bound; 32 CPUs idle |
| postprocess | 10 m 35 s | 267.5 min | 32 CPUs / 96 partitions; CVO sorting 4.0 min single-threaded |
| total | 11 h 35 m 46 s | — | GPU idle for 57 % of the allocation |

### Run D so far (job 12447879, 64 CPUs, CPU partition)
| Metric | Value | Source |
|---|---|---|
| Elapsed at snapshot | 15 h 34 m 49 s | sacct |
| AveCPU | 40 d 22 h 10 m = 982 core-h → 63.0 cores busy | sstat |
| Candidates / examples written | 53,969,760 / ≈ 46 M (per-task last log lines; 808 k–890 k per task, well balanced) | results/examples.log |
| Throughput | 962 candidates/s total; 15.0–15.5 per shard; 132.6 s mean per 2,048-candidate chunk (27,274 chunks, max 273 s) | examples.log |
| Task 0 time course | 78 k candidates at 1 h 20 m → 798 k at 14 h 49 m, linear (no slowdown, no stall) | examples.log |
| Reads / writes | 21.34 TB read, 0.88 TB written (793 GB of `make_examples_somatic.tfrecord-*`) | sstat, `ls` |
| MaxRSS | 109.0 GB (1.70 GB per shard) | sstat |
| Node | frnt128: 2 × 32 physical cores, Intel AVX-512, 498 GB, CPULoad 69.7 | scontrol show node |

**CPU-bound, not I/O-bound.** 63 of 64 cores are busy continuously; the GPFS read rate is only ≈ 380 MB/s aggregate (6 MB/s per shard) and the write rate ≈ 16 MB/s. There is no sign of waiting: the per-chunk time is stable across 15 h and across tasks.

**Why each shard reads the whole BAM [DeepVariant design — not profiled here, but consistent with both jobs' I/O counters].** `make_examples` assigns 1 kb partitions round-robin to tasks, so every task performs index seeks across the entire genome and decompresses nearly every BGZF block. Route A's trace `rchar` = 12.5 TB ≈ 32 × 415 GB; D's 21.3 TB so far against 64 × 434 GB = 27.8 TB. This is why read volume scales with `--num_shards`, and why **[provisional]** D's stage-1 progress can be estimated at 77 % (if all reads are BAM) to 85 % (if, as in A, ~10 % of `rchar` belongs to the later stages) → projected stage-1 wall ≈ 18–20 h, finishing around 12:30–14:40 local on 2026-11-11. Genome position of the running shards could not be read directly (login-node `ssh frnt128` returns nothing; no overlap step was attached to the job, by instruction).

**Projected remainder [provisional].** Final candidates ≈ 63–70 M (1.9–2.1× A); final examples ≈ 54–60 M. Inference at A's A30 rate (1,522 ex/s) ≈ 10–11 h; at the RTX6000 rate measured locally (700 ex/s pilot; 770 ex/s in `resource_sizing_20261106/DEEPSOMATIC.md`) ≈ 19–21 h. Postprocess ≈ 20–30 min (A: 10.6 min for 27 M examples). End-to-end ≈ 29–32 h on an A30, 38–42 h on an RTX6000, excluding queue waits between the three dependent jobs.

## 4. Root-cause findings

### F1 — Input BAM composition (adapter read-through soft clips) is the dominant cost driver. Evidence strength: high for the association; the mechanism is DeepVariant design knowledge.
- Per-shard, per-candidate CPU is 3.1× (chunk times) to 3.5× (AveCPU vs user time) higher in D than in A, on physical cores in both cases; shard count alone cannot explain it (16-shard pilot shows the same slow rate, 20.5 cand/s/shard, on an unloaded 20-core node).
- The same 10 kb windows contain nearly identical non-duplicate read counts (e.g. chr2: 12,649 vs 13,442) and identical depth (64× vs 63.6×), but **84.2 / 85.7 / 84.4 % soft-clipped reads in the REDUX tumour vs 8.2 / 7.7 / 13.9 % in the Parabricks tumour**, with 78–83 % of tumour inserts < 151 bp in both. The normal (inserts < 151 bp in 4–5 %) shows a much smaller gap (7–12 % vs 3–6 %).
- Provenance: Parabricks `@PG CL` reads `--in-fq .../fastp/trimmed/KQ0205-DNA.trimmed_1.fast…`; bwa-mem2 `@PG CL` reads `KQ0205-DNA_S2_L00{1..4}_R1_001.fastq.gz`. oncoanalyser 3.0.0 (`7c74c87a43`) has `modules/local/fastp/{umi,split}` only.
- Mechanism **[DeepVariant design — not profiled here]**: the realigner's window selector counts soft-clipped bases, mismatches and indels as evidence; with > 80 % of reads clipped, nearly every 1 kb partition triggers local de Bruijn assembly + re-alignment, the most expensive step of `make_examples`. Partially aligned adapter bases near read ends also add low-VAF mismatches, which at the FFPE_WGS threshold `vsc_min_fraction_snps = 0.029` (≈ 2 alt reads at 64×) become candidates — consistent with ≈ 2× the candidate count. SAGE (the HMF caller) is designed for untrimmed, soft-clipped input; DeepVariant's documented workflows assume adapter-trimmed reads (not verified against DeepSomatic documentation in this review).
- Not the cause: hyperthreading (ThreadsPerCore=1), memory (109/256 GB), I/O wait (63/64 cores busy), reference contig count (57 alt-contig reads total), base-quality binning (identical 37/25/11 bins in both), BQSR (not applied in either).

### F2 — `--num_shards 64` doubled the GPFS read volume and gave sub-linear scaling. Evidence: `sstat` MaxDiskRead 21.3 TB; A trace `rchar` 12.5 TB; per-shard 20.5 → 15.5 cand/s from 16 → 64 shards (different nodes, so partly confounded by CPU model: pilot node frnt153 Intel 20-core vs frnt128 Intel 64-core).

### F3 — Inference GPU type is left to the scheduler. The 72 h request puts job 12447880 in `gpubase_14days`, whose nodes (sinfo, 2026-11-11) are 1 × dgx `a100:8`, 1 × `rtx6000:4`, 2 × `rtx6000:2`, 6 × `a30:2`. Measured local rates: A30 ≈ 1,500 ex/s (two completed runs), RTX6000 ≈ 700–770 ex/s. A100/L4/L40S/V100 are present on Frontenac but **no local DeepSomatic throughput measurement exists for them** (the RTX6000_DEEPSOMATIC_CHECK doc explicitly forbids extrapolating).

### F4 — Over-sized memory/postprocess requests delay scheduling without speeding anything up. 256 GB for examples (observed 109 GB) and 64 CPU/256 GB/12 h for a ~25 min postprocess.

### F5 — The runner's resource model was calibrated on the wrong basis. `deepsomatic_full.py::resources()` scales the 20 Mb pilot (chr17/chr19 p-arms, 31.8 k candidates/Mb — gene-dense, ~1.5–3× the genome average) by genome length and assumes 4× speed-up from 16 → 64 CPUs; it produced 60 h / 72 h / 12 h requests that are generous in time but say nothing about where the time goes.

### F6 — The three-stage design is correct and should be kept. Route A idled an A30 for 6.6 of 11.6 h; D runs stage 1 on the CPU partition. The cost is two dependency hops (each adds a queue wait).

## 5. Tuning recommendations, ranked by expected impact

| # | Change | Expected effect | Risk / caveat | Basis |
|---|---|---|---|---|
| 1 | **Decide, in the data-wiring contract, how adapter read-through is handled before DeepSomatic** — either trim adapters before the oncoanalyser alignment (fastp as in route A; check whether oncoanalyser's `fastp/split` path or an upstream pre-trim step is acceptable to the HMF toolchain) or accept the untrimmed REDUX BAM and budget for it. | Trimmed-like input: `make_examples` CPU ≈ 177 core-h (A) instead of ≈ 1,200 core-h **[provisional]**; candidates/examples halve, so inference halves too. | Changes the alignment that SAGE/PURPLE/ESVEE also consume; HMF's tools expect their own alignment conventions, so this must be validated on the whole-pipeline level, not just for DeepSomatic. Do **not** work around it with `--realign_reads=false` or by altering `vsc_*` thresholds without a calibrated comparison. | Evidenced (F1) |
| 2 | **Type the inference GPU and shorten the limit**: `--gres=gpu:a30:1 --time=16:00:00` (fits `gpubase_24hrs`, which also has 3 idle a30:2 nodes at snapshot). Keep `--batch_size 256` (canaried). | Avoids a 2× slower RTX6000 placement; opens more nodes. ≈ 10–11 h for 54–60 M examples **[provisional]**. | A100 (dgx) would likely be faster but is unmeasured and mostly allocated; run a bounded canary (the 4,096-example RTX6000 check pattern) before relying on it. L4 (24 GB) is also unmeasured. | Evidenced (A30/RTX6000 rates) |
| 3 | **Region-split `make_examples` into N parallel CPU jobs** (e.g. 4–8 jobs × 32 CPUs, each `--regions` = a contig group balanced by BAM bytes, each with its own `--intermediate_results_dir`), then one `call_variants` over all example files, then one `postprocess_variants`. | Wall ≈ total core-h / (N × 32) + imbalance: untrimmed input ≈ 1,200 core-h → ≈ 9–10 h with 4 × 32 cores, ≈ 5 h with 8 × 32 **[provisional]**. Each job reads only its contigs, so total BAM reads drop back to ≈ 32 × BAM. Each job's examples (≈ 100–200 GB) fit `$SLURM_TMPDIR` and can be copied to GPFS at the end. | (a) `call_variants --examples` and `postprocess_variants --infile/--small_model_cvo_records` must accept comma-separated lists of sharded specs — **verify with `--help` inside the 1.10.0 container before building on it**. (b) The runner currently rejects any `--regions`; the contract's "no region filter" intent is preserved only if the union of groups is exactly the full contig list, including unplaced/decoy contigs — assert this in code. (c) Per-group CVO/sidecar naming must stay unique. (d) More jobs = more queue waits. | Design; mechanism evidenced by F2 |
| 4 | **Split inference across 2 GPU jobs on disjoint example-file lists** (e.g. shards 0–31 and 32–63), postprocess consumes both outputs. | Halves inference wall (≈ 5–6 h on two A30s). | Same flag-syntax verification as #3; the dnaseq_nf lesson on `MirroredStrategy` means each job must still see exactly one GPU (`CUDA_VISIBLE_DEVICES` is already enforced by the runner). | Design |
| 5 | **Right-size memory**: examples 2.5 GB × shards (→ 160 GB for 64, 96 GB for 32); inference 64–128 GB (dnaseq_nf uses 120 GB after a 48 GB host-memory failure on Narval); postprocess 32 CPUs / 96 GB / 2 h with `num_partitions = 3 × cpus`. | Faster scheduling; no runtime change. | Keep headroom for higher-depth cases (KQ0294 127×). | Evidenced (RSS, A's postprocess) |
| 6 | **Keep `--num_shards = allocated CPUs`, but prefer 32 shards per job** (and more jobs) over 64 shards in one job. | Removes the 2× BAM re-read and the observed 25 % per-shard slowdown. | None beyond #3. | Evidenced (F2) |
| 7 | **Intermediate results on `$SLURM_TMPDIR`** only in the region-split layout (#3); in the single-job layout D's ≈ 0.8–1 TB of examples exceed the 462 GB node-local disk. | Lower GPFS small-write load; little wall-time gain because the stage is CPU-bound. | Copy-back step adds minutes; inference still reads from GPFS. | Evidenced (sizes) |
| 8 | **Thread hygiene for the small model inside `make_examples`** (e.g. `OMP_NUM_THREADS=1`, `TF_NUM_INTRAOP_THREADS=1`, `TF_NUM_INTEROP_THREADS=1` per shard). | Possibly removes the mild oversubscription (CPULoad 69.7 on 64 cores). | **Hypothesis — not evidenced**; test on the 20 Mb pilot first. | Hypothesis |
| 9 | **Replace the runner's pilot-scaled time model** with an input-based one: ≈ 2.8 core-h per BAM-pair GB for untrimmed REDUX input (≈ 1,200 core-h / 434 GB) vs ≈ 0.43 core-h/GB for trimmed input (177 / 415) **[provisional, single case each]**; or extend the `resource_sizing_20261106` formula (`T_examples = 0.03 × G h` at 32 CPUs) with a ≈ 3× multiplier for untrimmed input. Better still, derive the request from the first 30 min of shard-chunk times. | Realistic walltime requests; fewer 60 h/72 h jobs in the 14-day partitions. | Needs ≥ 2 more cases to calibrate. | Evidenced (this review + sizing doc) |
| 10 | Housekeeping: the 433 GB `ds_tmp` of route A (work/22/70ba71…) and the ≈ 800 GB D intermediates should have an explicit retention rule; `examples.%j.log` is empty by design (the runner writes `results/<stage>.log`). | Scratch hygiene; monitoring clarity. | None. | Observed |

Not recommended: running `make_examples` on a GPU node (A's design: 57 % GPU idle); untyped GPU requests for multi-hour inference; deriving batch size from GPU memory (dnaseq_nf policy, keep 256).

## 6. Proposed per-stage resource table (30–60× FFPE tumour / 30–80× normal, Frontenac)

Two input scenarios, because scenario is the biggest factor. "Untrimmed REDUX" = today's oncoanalyser native alignment; "trimmed-like" = adapter-trimmed reads (route A behaviour). All projections are **[provisional]**, scaled from KQ0205 (64×/78×, 434 GB BAM pair); scale CPU-hours roughly with BAM pair size.

| Stage | Layout | CPU | Mem | Walltime request | GPU | Expected wall, untrimmed REDUX | Expected wall, trimmed-like |
|---|---|---:|---:|---:|---|---:|---:|
| make_examples | single job (current) | 64 | 160 GB | 30 h | none (CPU partition) | 18–20 h | 4–5 h |
| make_examples | **region-split 4 jobs** | 4 × 32 | 96 GB each | 14 h each | none | ≈ 9–10 h | ≈ 2 h |
| make_examples | region-split 8 jobs | 8 × 32 | 96 GB each | 8 h each | none | ≈ 5 h | ≈ 1.5 h |
| call_variants | one job | 8–16 | 96–128 GB | 16 h (≤ 24 h) | **`gpu:a30:1`**, batch 256 | 10–11 h (54–60 M ex.) | 5 h (27 M ex.) |
| call_variants | two jobs, disjoint shards | 2 × 8 | 96 GB each | 8 h each | 2 × `gpu:a30:1` | ≈ 5–6 h | ≈ 2.5 h |
| postprocess_variants | one job | 32 | 96 GB | 2 h | none | ≈ 25 min | ≈ 11 min |
| **End-to-end (no queue)** | current D | | | | | **≈ 29–32 h** (A30) / 38–42 h (RTX6000) | ≈ 10 h |
| | 4-way split + 1 A30 | | | | | ≈ 20–21 h | ≈ 7.5 h |
| | 8-way split + 2 A30 | | | | | ≈ 11 h | ≈ 4.5 h |

Scratch: budget ≈ 2 GB of GPFS per million candidates for examples (A: 433 GB for 27 M examples; D: 793 GB for ≈ 46 M) plus the call-outputs; delete intermediates once the VCF is accepted.

## 7. What could not be verified

- Genome position / exact completion fraction of job 12447879 (no node access without attaching to the job); the 77–85 % figure is inferred from I/O counters.
- The causal mechanism soft clips → realignment windows → CPU is from DeepVariant design knowledge; no profiler was run. The association (same depth, 10× the soft-clip rate, 3.1–3.5× CPU per candidate, ~2× candidates) is measured.
- Contribution of the reference (HMF vs GDC decoys) and of REDUX consensus/unmapping to the candidate count — not isolated; the 306 M unmapped reads are not read by region queries and the alt-contig read count is 57.
- A100 / L4 / L40S / V100 inference throughput on Frontenac — no local measurement; the dnaseq_nf Narval figure (813 ex/s, GPU type not stated in the comment) is not transferable.
- Whether `call_variants`/`postprocess_variants` 1.10.0 accept comma-separated multi-pattern inputs (needed for recommendations 3–4) — check `--help` in the container.
- Final candidate/example counts and the inference/postprocess times of run D (still running).
- BAM comparisons are from three 10 kb windows (chr2:30.00–30.01 M, chr7:55.00–55.01 M, chr17:7.60–7.61 M) and the idxstats totals, not whole-BAM statistics.

## Appendix A — files, logs and jobs inspected (all read-only)

Slurm: `sacct -j 12447879,12447880,12447881,12414711,12446769` (Elapsed, AllocCPUS, ReqMem, MaxRSS, AveCPU, TotalCPU, NodeList, ReqTRES, AllocTRES, State, Submit/Start/End, Timelimit); `sstat -j 12447879 --allsteps`; `scontrol show job 12447879 / 12447880`; `scontrol show node frnt128 frnt144 frnt153 frnt118`; `sinfo -p gpubase_6hrs,gpubase_24hrs,gpubase_14days,gpu-L4`; `squeue -u hpc6049`.

Run D root `/global/scratch/hpc6049/wigits_runs/KQ0205_full_deepsomatic_native_20261110/`: `case.json`, `resources.json`, `submission.json`, `batch256_update.json`, `batch256_native_dry_run.log`, `input_header_audit_20261111.json`, `examples.12447879.log` (0 bytes), `package_v0.2.0/{deepsomatic_full.py,deepsomatic_full.sbatch,release.json}`, `results/{commands.json,commands.before_batch256.json,examples.json,native_dry_run.log,examples.log (29,503 lines)}`, `results/intermediate/` (128 files, 733 GiB `du`, 792.5 GB examples).

Pilot `/global/scratch/hpc6049/wigits_runs/deepsomatic_seam_20261110/`: `case.json`, `caller.12446769.log`, `caller_test_v0.1.0/{execution.json,runner.log,logs/make_examples_somatic.log,logs/call_variants.log,logs/postprocess_variants.log}`.

Route A `/global/scratch/hpc6049/mohcc_runs/gi_frontenac_20261103/KQ0205/`: `work/22/70ba719541a0567e75b335e070cb85/{.command.sh,KQ0205.deepsomatic_ffpe.run.log (21,867 lines),ds_tmp/}`, `launch/nextflow_trace.txt` (DEEPSOMATIC_FFPE row), `inputs/*.bam` sizes.

BAMs (header, `idxstats`, three 10 kb windows via `samtools 1.20`): `KQ0205_full_native_20261109_onco3.0.0/results/KQ0205_full_native/alignments/{KQ0205-DNA,KQ0205-gDNA}/*.redux.bam`, `gi_frontenac_20261103/KQ0205/inputs/{KQ0205-DNA,KQ0205-gDNA}.analysis_ready.bam`; BamMetrics `results/KQ0205_full_native/bamtools/*/*.bam_metric.summary.tsv`; REDUX `*.redux.{bqr,duplicate_freq,jitter_params}.tsv` (headers only).

dnaseq_nf (`/global/project/hpcg6049/somatic_pipeline/dnaseq_nf`, checkout `run/kq0043-reporting-validation-20260831` @ 2b893ae): `modules/deepsomatic_ffpe.nf`, `modules/deepsomatic_call_examples.nf`, `conf/gpu.config` (envWhitelist, DEEPSOMATIC_FFPE / DEEPSOMATIC_CALL_EXAMPLES_* / DEEPSOMATIC_POSTPROCESS_EXAMPLES blocks, process_gpu clusterOptions), `nextflow.config` lines 169–237, `docs/LESSONS.md` (DeepSomatic entries), `docs/qc/reviews/35_v1.4.24_fastq2deepsomatic_full_wgs.md`, `21_deepsomatic_native_postprocess_recovery.md`.

somatic_pipeline_sm: `docs/qc/resource_sizing_20261106/DEEPSOMATIC.md`, `docs/rac2027/RTX6000_DEEPSOMATIC_CHECK.md`.

wigits-clinical (`hartwig_adoptation`): `git log -- workflow/deepsomatic_full.py` (d233320 "Set DeepSomatic inference batch size to 256", b7510cd "Run whole-genome DeepSomatic in dependent CPU and GPU stages"); flag grep of `workflow/deepsomatic_full.py`. oncoanalyser 3.0.0 clone `~/.nextflow/assets/.repos/nf-core/oncoanalyser/clones/7c74c87a43…/modules/local/fastp/`.

## Appendix B — key raw numbers

- A make_examples: `real 387m0.204s user 10634m37.678s sys 154m5.561s`; call_variants `real 297m12.754s user 908m52.820s`; postprocess `real 10m35.397s user 267m32.286s`; last inference line `Predicted 27149056 examples in 106051 batches [0.066 sec per 100]`.
- D sstat (15 h 34 m 49 s): `AveCPU 40-22:10:29, MaxRSS 108989868K, AveRSS 104869880K, MaxDiskRead 21344422831484, MaxDiskWrite 882626976378`.
- D per-task candidates at snapshot: min 808,126 (task 54), max 890,032 (task 36); sum 53,969,760.
- Window stats (non-duplicate reads; soft-clipped %, mean clip bp, insert < 151 bp %): chr2 T-REDUX 12,649 / 84.2 / 59.7 / 79.8; T-Parabricks 13,442 / 8.2 / 49.6 / 82.8; chr7 12,977 / 85.7 / 60.2 / 82.5 vs 13,269 / 7.7 / 54.5 / 84.9; chr17 14,605 / 84.4 / 58.9 / 77.6 vs 14,782 / 13.9 / 43.6 / 80.1. Normal: chr2 7.3 % vs 2.9 %, chr7 8.1 % vs 3.9 %, chr17 11.6 % vs 5.6 %.
