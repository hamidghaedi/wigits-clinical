# Full-genome DeepSomatic runtime investigation

Read-only inspection on 11 November2026, around00:20 EST. Running jobs and
frozen packages were not changed. Full-run component0.2.0 is `b7510cd`.

## Current execution

Job12447879 is generating examples after6h20m onfrnt128, using64 CPUs and
256GiB. Jobs12447880 (GPU inference) and12447881 (postprocessing) remain
dependency-pending. These are the three native caller phases; there is no
new alignment, BAM conversion or full-file validation in this execution.

Slurm accounting indicates approximately63 busy CPU equivalents and103GiB
maximum RSS. A lightweight process inspection found64 substantial
make_examples processes. One inspected process had289 threads, of which
one was running and288 sleeping. The thread count alone is therefore not
evidence of active oversubscription. The workers continue producing examples.

## Historical full-genome comparisons

Times below are from the successful attempt's native caller log, not requested
walltime or failed earlier attempts. All historical runs used32 shards.

| Case | CPU examples | GPU inference | Postprocessing | Whole task |
|---|---:|---:|---:|---:|
| KQ0205 | 6h27m | 4h57m | 10m35s | 11h35m45s |
| KQ0293 | 6h46m | 7h51m | 14m37s | 14h52m37s |
| KQ0294 | 10h02m | 2h46m | 9m15s | 12h58m04s |
| KQ0299 | 9h22m | 8h29m | 17m26s | 18h10m04s |

At6h16m the new KQ0205 log had approximately21million candidates and18million
examples; the old finished run's last periodic counters totaled approximately
33.7million candidates and27.1million examples. These counters are neither final
variant counts nor a fraction of genome completed. They cannot provide a reliable
completion ETA, particularly across different alignments and references.

The current median logged time per2,000-candidate increment was135.4s per
shard versus44.42s in the old run. Shards, loci, reference and hardware differ;
this is observational evidence of slower example production, not a controlled
threefold performance comparison.

## Inputs and settings

- Correct full FASTQ-run REDUX tumour and matched-normal BAMs, not the regional
  pilot or custom-reference route. No regions restriction or downsampling.
- Official Google DeepSomatic1.10.0 GPU image and FFPE_WGS model match the prior
  caller. Both runs enable the FFPE small model, SNP53/indel36 thresholds,
  alt-aligned pileups and automatic read trimming.
- Both current BAM headers have all3,367 reference contigs with matching names
  and lengths, four read groups and their intended sample labels. Neither has
  comparable contig MD5 tags; header inspection does not establish sequence
  identity. Reference provenance is recorded by the native source run.
- Combined BAM size is403.91GiB versus386.40GiB previously, a4.53% increase.
  New bamtools coverage is64.082x tumour and77.507x normal, versus prior
  metrics63.597568x and74.787936x. Different metric tools limit exact comparison;
  these values do not suggest accidental doubling of the input coverage.
- Previous BAMs used GRCh38.d1.vd1; current native REDUX BAMs use the HMF masked
  reference. Old KQ0205 ran on AMDfrnt144, current CPU phase on Intelfrnt128.
  These changes have not been separated experimentally. Busy CPUs and slower
  production do not establish the specific computational bottleneck.
- Current limits60h/72h/12h are conservative allocation ceilings derived from
  the regional pilot, not expected runtimes. They are not calibrated predictions
  from the historical full-genome cases.

## Concrete setting discrepancy: inference batch size

The current generated call_variants command omits `--batch_size` and therefore
uses the native1024 default. Historical full runs explicitly used256. The DNA
repository's nextflow.config documents an A30 out-of-memory failure at1024
(job12070537), followed by a128-versus256 canary with matching genotype argmax
(jobs12071126/12071127, comparison12071133). This is not a256-versus1024
equivalence test. The recent20Mb pilot did succeed with1024 on RTX6000, so a
future failure is not certain.

This omission cannot explain the running CPU phase. It is a separate GPU-memory
risk to address before inference using the established256 setting while
preserving all completed examples. No queued job or package was changed in this
investigation. No claim is made that the CPU slowdown is solved.

## Evidence locations

Current root:
`/global/scratch/hpc6049/wigits_runs/KQ0205_full_deepsomatic_native_20261110/`.
Read `case.json`, `resources.json`, `results/commands.json`,
`results/examples.log` and `input_header_audit_20261111.json`.

Historical root:
`/global/scratch/hpc6049/mohcc_runs/gi_frontenac_20261103/`.
Successful work directories beneath each case:

- KQ0205: `work/22/70ba719541a0567e75b335e070cb85/`
- KQ0293: `work/b0/397fd99b04622a09056d2de00156c4/`
- KQ0294: `work/66/7fea4e7e447d26e82deeceed017e9b/`
- KQ0299: `work/98/460f10671ee7f956090840932a7818/`

Use each `.command.sh`, `*.deepsomatic_ffpe.run.log` and
`launch/nextflow_trace.txt`. Preserve successful-attempt boundaries.

Checks: read-only squeue/sstat/sacct; bounded log inspection; remote process
inspection; samtools1.22.1 header-only comparison against the native FASTA
dictionary. No heavy scan or new compute experiment was submitted.
