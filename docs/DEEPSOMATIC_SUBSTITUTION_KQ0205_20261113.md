# DeepSomatic → oncoanalyser `sage_dir` substitution — KQ0205 (2026-11-13)

Status: executed end-to-end on Frontenac; oncoanalyser pass 2 produced PAVE 1.9, PURPLE 4.4, LINX, CHORD, SIGS, CUPPA,
PEACH, LILAC, QSEE and ORANGE 5.0.1 outputs from the **complete DeepSomatic 1.10.0 PASS call set** (14,436 records) via the
sage-compat adapter (`workflow/sage_compat/`, commit 9a47651). Reviewer: Codex. Contract: `docs/DEEPSOMATIC_WIRING_CONTRACT.md`.

Run directory: `/global/scratch/hpc6049/wigits_runs/KQ0205_full_deepsomatic_substituted_20261113/` (`adapter/`, `results/KQ0205_full_deepsomatic_substituted/`, `work/`, `logs/`, `launch.sh`, `samplesheet.csv`).
Comparator: route C (native SAGE), `/global/scratch/hpc6049/wigits_runs/KQ0205_full_native_20261109_onco3.0.0/results/KQ0205_full_native/`.

## 1. Scope of this run

**DNA-only, two-pass Option A** (contract §5). Route C's REDUX, AMBER, COBALT, ESVEE and VirusInterpreter outputs were supplied as
existing inputs; SAGE_GERMLINE, PAVE (somatic + germline), PURPLE, LINX, CHORD, SIGS, CUPPA, PEACH, LILAC, CIDER, BAMTOOLS,
TEAL, QSEE and ORANGE ran natively. The only replaced output is the SAGE_SOMATIC discovery VCF.

Deliberately excluded or disabled, and why:

| Item | Setting | Reason |
|---|---|---|
| RNA (`KQ0205-RNA`) | no RNA rows in the samplesheet → ISOFOX, STAR, RNA CIDER, SAGE_APPEND not run | SAGE_APPEND re-genotypes the PURPLE somatic VCF on the RNA BAM and rebuilds read contexts from `INFO/RC_INFO` (`sage/vcf/CandidateSerialisation.java`); the adapter cannot produce `RC_INFO` without fabrication, so RNA depth per variant is impossible for this VCF. ORANGE therefore has no RNA section. |
| `SAGE_APPEND_*` | `ext.when = false` in `workflow/conf/deepsomatic_substituted.config` | Guard against the same failure if RNA rows are ever added to this samplesheet. |
| NEO | `--processes_exclude neo` | NEO parses `FORMAT/RC_CNT` as a 6-element list and throws when absent (`neo/…/SomaticVariant.java:69-80`); not emitted by the adapter. |
| SAGE_VISUALISER | `--processes_exclude sage_visualiser` | Re-derives read contexts from the somatic VCF and BAMs; not needed by ORANGE (optional input) and untested on a VCF without `RC_*`. |
| VIRUSBREAKEND | existing `virusinterpreter_dir` row (route C output) | Independent of SAGE; 13 h 44 m in route C. Reused, not re-run. |
| TINC | n/a | Adapter emits no `TINC_RECOVERED`; `fit_variant.tsv` has no consumer outside SAGE. |

Consumer coverage for the somatic small-variant VCF in this run:

| Consumer | Needs from the somatic VCF | This run | Note |
|---|---|---|---|
| PAVE 1.9 | PASS filter, `TIER` (PON thresholds per tier), `MH`/`REP_*` (realignment), `LPS` (phasing) | **ran** (degraded: no `LPS` → no phased MNV/inframe reinterpretation) | Default `FilterType.PASS`: non-PASS records are dropped regardless of tier (see §4). |
| PURPLE 4.4 | `AD` both samples, `TIER`, `REP_C`, `TNC` (kataegis), `LPS` (phased reporting) | **ran** (degraded: no `LPS`) | 4,870 fitting variants (TRACE log). |
| LINX | PURPLE outputs only | **ran** | SV path untouched (route C ESVEE). |
| CHORD | PURPLE somatic VCF + SV VCF; indel context from the reference | **ran** | HR_proficient in both routes. |
| SIGS | PASS SNVs + `TNC` | **ran** | Allocation differs strongly (§6). |
| CUPPA | PASS SNVs/`TNC`, indels/`REP_*`, drivers, SVs | **ran** (DNA classifier only) | Top class changes (§6). |
| ORANGE 5.0.1 | PURPLE/LINX/… dirs; `-sage_dir` declared but unread | **ran** (JSON + PDF, DNA-only) | |
| NEO | `FORMAT/RC_CNT` | **excluded** | Would throw. |
| SAGE_APPEND (RNA depth) | `INFO/RC_INFO` | **excluded** (no RNA rows; also `ext.when=false`) | Needs genuine read contexts. |
| SAGE_VISUALISER | read contexts | **excluded** | Optional for ORANGE. |
| VirusInterpreter | n/a (BAM-based) | reused from route C | |

Path to closing the `RC_INFO`/`RC_CNT` gap: Codex's native-context work on this branch (commit f00524a, `docs/DEEPSOMATIC_EVIDENCE_ADAPTER.md`,
`workflow/merge_sage_evidence.py`, `workflow/build_sage_evidence.sbatch`, `workflow/patches/sage-5.0.2-evidence-export.patch`) measures
genuine SAGE read contexts and `RC_CNT` for DeepSomatic alleles with an isolated SAGE evidence-export overlay (180/183 alleles matched on the
regional pilot; `test_evidence_consumers.py` exercised SAGE_APPEND, NEO's count reader and PAVE on that subset). Unresolved before those
fields can be carried into this `sage_dir`: whole-genome execution of the evidence-export stage (only the ≤1,000-allele probe scope exists),
alleles whose candidates SAGE never discovers (the chr17/chr19 indels needed exact-CIGAR seeds), count semantics (SAGE fragment support vs
DeepSomatic read counts — e.g. 0 vs 91 ALT at the chr17 deletion), phasing (`LPS`) for DeepSomatic-only calls, and the owner's decision on
whether a merged file is still "the DeepSomatic population". Until then the substituted route is DNA-only without NEO.

## 2. Commands, jobs, wall times

| Stage | Slurm job | Wall | Outcome |
|---|---|---|---|
| Adapter attempt 1 | 12466683 | 5 s | failed: Java self-test assertion (test expectation wrong; code unchanged) |
| Adapter attempt 2 | 12466703 | 2 m 06 s | failed: `bcftools norm` rejects DeepSomatic `FORMAT/NAF` (declared Number=R, written A); fixed by `bcftools annotate -x FORMAT/NAF` |
| Adapter attempt 3 (called-GT rule) | 12466732 | 6 m 22 s | succeeded; 14,436 retained, every panel non-PASS record dropped because DeepSomatic GERMLINE records carry GT 0/0 — outputs archived as `adapter_run1_calledgt_rule/` |
| **Adapter (final)** | **12466791** | **5 m 22 s** | 15,193 records written; validation clean |
| samtools verification of normal support | 12466856 | 11 s | confirmed pileup at 6 sites (§5) |
| Pass-2 head, first launch | 12466868 | 1 h 18 m (cancelled) | PAVE_SOMATIC (12466882, 1 m 43 s) completed; five 32-CPU tasks stayed pending >1 h on a saturated cluster; head and its pending tasks cancelled by me (mine only) |
| **Pass-2 head, relaunch `-resume`** | **12467530** | started 15:04:55 local; ORANGE written 16:30 | PAVE_SOMATIC cached; see task table |

Pass-2 tasks (relaunch; Slurm id, elapsed, CPUs/mem): SAGE_GERMLINE 12467656 1 m 37 s 8/32G · CIDER 12467657 6 m 17 s 8/32G ·
PAVE_GERMLINE 12467729 40 s 6/36G · PURPLE 12467838 53 s 8/48G (TRACE) · PEACH 12468399 13 s · LINX_GERMLINE 12468401 30 s ·
CHORD 12468402 1 m 07 s · SIGS 12468403 13 s · LINX_SOMATIC 12468405 14 s · CUPPA 12468413 36 s · BAMTOOLS 12467655 22 m 01 s 16/48G,
12467658 7 m 45 s 16/48G · LINX_VISUALISER 12468412 3 m 41 s · LINXREPORT 12468546 13 s · LILAC 12468400 33 s 16/48G ·
QSEE 12468563 23 s · ORANGE 12468583 8 s · TEAL_PREP 12467659 (running at report time, 16/48G; TEAL is not an ORANGE input).
Wall time is dominated by queue waits, not compute.

Commands:
```
sbatch … workflow/sage_compat/run_sage_compat_adapter.sbatch workflow/sage_compat/configs/KQ0205_deepsomatic_substituted_20261113.json
export WIGITS_WORKFLOW_DIR=$WIGITS_REPO/workflow
bash workflow/run_oncoanalyser_frontenac.sh KQ0205_full_deepsomatic_substituted_20261113 \
     samplesheets/KQ0205_full_deepsomatic_substituted_20261113.csv \
     -c workflow/conf/deepsomatic_substituted.config --processes_exclude neo,sage_visualiser
```
`launch.sh` in the run directory records the exact Nextflow command (nf-core/oncoanalyser -r 3.0.0, `-profile apptainer`,
frontenac.config + frontenac_limits.config + reference_data.config + deepsomatic_substituted.config, `--mode wgts --genome GRCh38_hmf`).
Samplesheet: redux_dir (tumour, normal → route C `alignments/<sample>/`), sage_dir (tumour → `adapter/sage_dir/`), amber_dir, cobalt_dir
(tumor), esvee_dir (tumor_normal), virusinterpreter_dir (tumor); sha256 `949afc91ff5087faf165b5af76919d2ea0cbe5d89c0be7d3291ac19640fd42ff`.

Verification of the seam: `.nextflow.log` creates `SAGE_CALLING:SAGE_SOMATIC` with `maxForks=0` and never submits a task for it (no
`Submitted process … SAGE_SOMATIC` in either head log); the PAVE_SOMATIC work dir (`work/66/939313…`) stages
`KQ0205-DNA.sage.somatic.vcf.gz` as a symlink resolving to `/global/scratch/hpc6049/wigits_runs/KQ0205_full_deepsomatic_substituted_20261113/adapter/sage_dir/KQ0205-DNA.sage.somatic.vcf.gz`; PAVE logged
`loaded … 14435` records across chromosomes and PURPLE logged `loaded 14435 somatic variants from pave_somatic/…`.

Versions: DeepSomatic 1.10.0 (FFPE_WGS model), bcftools 1.22, hmftools-sage 5.0.2 image (classes used: `ReferenceData`, `TierSelector`,
`RepeatInfo`, `Microhomology`; OpenJDK 21-internal), PAVE 1.9, PURPLE 4.4 (build 2026-06-08), ORANGE 5.0.1, oncoanalyser 3.0.0
(`7c74c87a43`), Nextflow 26.04.4, Apptainer 1.3.5. Resources `hmf_pipeline_resources.38_v3.0.0--8` (SHA256s in the VCF header).

## 3. Adapter VCF header (excerpt; `##contig` lines omitted)

```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=GERMLINE,Description="DeepSomatic: Non somatic variants">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (tumour: DeepSomatic AD; normal: pileup of the normal REDUX BAM)">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allelic frequency (tumour: DeepSomatic VAF; normal: pileup ALT/DP)">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth (tumour: DeepSomatic DP; normal: usable reads spanning POS in the normal REDUX BAM)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype (tumour: DeepSomatic call; normal: not called, ./.)">
##INFO=<ID=DS_FILTER,Number=1,Type=String,Description="Original DeepSomatic FILTER for retained non-PASS records (HOTSPOT/PANEL tier only)">
##INFO=<ID=DS_NAD,Number=R,Type=Integer,Description="DeepSomatic FORMAT/NAD of the original record projected to this allele (normal counts as computed by DeepSomatic; informational)">
##INFO=<ID=DS_NDP,Number=1,Type=Integer,Description="DeepSomatic FORMAT/NDP of the original record (informational)">
##INFO=<ID=DS_ORIGINAL_ALLELE,Number=1,Type=String,Description="Original DeepSomatic record when the representation changed (bcftools norm old-rec-tag): CHR|POS|REF|ALT|USED_ALT_IDX">
##INFO=<ID=MH,Number=1,Type=String,Description="Microhomology (SAGE Microhomology.findHomology on reference/ALT haplotype)">
##INFO=<ID=REP_C,Number=1,Type=Integer,Description="Repeat sequence count (reference-only RepeatInfo.findMaxRepeat as SAGE setMaxRefRepeat)">
##INFO=<ID=REP_S,Number=1,Type=String,Description="Repeat sequence">
##INFO=<ID=TIER,Number=1,Type=String,Description="Tier: [HOTSPOT, PANEL, HIGH_CONFIDENCE, LOW_CONFIDENCE] (SAGE 5.0.2 TierSelector with native resources)">
##INFO=<ID=TNC,Number=1,Type=String,Description="Tri-nucleotide context (reference, centred on POS)">
##sage_compat_adapter_annotator_version=0.1.0
##sage_compat_adapter_date=2026-11-13T18:23:46Z
##sage_compat_adapter_filter_policy=FILTER=PASS for DeepSomatic PASS; non-PASS records retained only when TIER is HOTSPOT or PANEL (DS_FILTER carries the original filter); RefCall/NoCall/other records dropped
##sage_compat_adapter_normal_bam=/global/scratch/hpc6049/wigits_runs/KQ0205_full_native_20261109_onco3.0.0/results/KQ0205_full_native/alignments/KQ0205-gDNA/KQ0205-gDNA.redux.bam
##sage_compat_adapter_normal_pileup=bam=/global/scratch/hpc6049/wigits_runs/KQ0205_full_native_20261109_onco3.0.0/results/KQ0205_full_native/alignments/KQ0205-gDNA/KQ0205-gDNA.redux.bam;min_mapq=10;min_baseq_substitutions=20;excluded=unmapped,secondary,supplementary,duplicate,qcfail;DP=usable reads spanning POS;indel/complex alleles matched between aligned anchors POS-1 and POS+len(REF) without base-quality filter
##sage_compat_adapter_normalisation=bcftools norm -f ref -m -any -c e --old-rec-tag DS_ORIGINAL_ALLELE; alleles retained only when the DeepSomatic tumour GT calls them
##sage_compat_adapter_not_emitted=LPS,LPS_RC,RC_INFO,RC_CNT,RC_QUAL,RC_JIT,RC_IPC,RC_NM,RC_MH,RC_REPC,RC_REPS,UMI_CNT,AED,AFBQ,AFMQ,ARCBQ,ARMQ,ASBQ,FSB,MUC,RSB,SAC,PON_COUNT,PON_MAX,PON_MULTI,PON_READS,GND_FREQ,MQF,TQP,MED,MSG,NEARBY_INDEL,TINC_RECOVERED (not fabricated; see README)
##sage_compat_adapter_reference=/global/project/hpcg6049/somatic_pipeline/data/references/oncoanalyser/3.0.0/genomes/GRCh38_hmf/GRCh38_masked_exclusions_alts_hlas.fasta dict_sha256=3fc739f04ff92156e9fa6d7ba3ed562f08e4f6a77ae3e37b9531d15394be91ae
##sage_compat_adapter_resources=hotspots=KnownHotspots.somatic.38.vcf.gz sha256=4381ed54ac42877d62fa4afe080463c476fc393a6516fe6ca4ac9aa3ea945aca;driver_gene_panel=DriverGenePanel.38.tsv sha256=1a5d01c17f8d954d323af71b30159efada46771e96512c93d6a5a6926cf7b8f3;high_confidence_bed=HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed.gz sha256=bd38220c2c8cadfa4f9b8ba2cfb50dc1b38fec758fb3b79cb44c16fdce59635b;ensembl_data_dir=/global/project/hpcg6049/somatic_pipeline/data/references/oncoanalyser/3.0.0/hmf_pipeline_resources.38_v3.0.0--8/common/ensembl_data
##sage_compat_adapter_slurm_job=12466791
##sage_compat_adapter_source_vcf=/global/scratch/hpc6049/wigits_runs/KQ0205_full_deepsomatic_native_20261110/results/deepsomatic.vcf.gz sha256=bf24cdea92b03cb2c332fd73031d1c3855dfaaabb1bac03d405b6516aad20132 deepvariant_version=1.10.0
##sage_compat_adapter_tools=sage_image=quay.io-biocontainers-hmftools-sage-5.0.2--hdfd78af_0.img;bcftools=1.22;java=openjdk version 21-internal 2023-09-19
##sage_compat_adapter_tumor_bam=/global/scratch/hpc6049/wigits_runs/KQ0205_full_native_20261109_onco3.0.0/results/KQ0205_full_native/alignments/KQ0205-DNA/KQ0205-DNA.redux.bam
##sage_compat_adapter_version=0.1.0 (wigits-clinical workflow/sage_compat; annotator uses hmftools sage 5.0.2 classes TierSelector,RepeatInfo,Microhomology)
##source=deepsomatic-1.10.0+sage-compat-adapter-0.1.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	KQ0205-gDNA	KQ0205-DNA
```

## 4. Allele population (ledger)

Source: `/global/scratch/hpc6049/wigits_runs/KQ0205_full_deepsomatic_substituted_20261113/adapter/ALLELE_LEDGER.tsv` (51,318 rows = every PASS record after splitting + every non-PASS record inside the
SAGE panel/hotspot regions), `ledger_summary.json`.

| Step | Count |
|---|---|
| DeepSomatic records (whole VCF) | 74,162,883 = PASS 14,436 · GERMLINE 4,443,224 · NoCall 1,330,600 · RefCall 68,374,623 (no LowQual) |
| PASS after `bcftools norm` (split / left-aligned) | 14,436 input → 19 multiallelic records split, 1 realigned; 19 uncalled second ALTs dropped (`dropped_alt_not_in_called_gt`) → **14,436 PASS alleles** (26 carry `multiallelic_split`, 1 `realigned_or_trimmed`, 2 both) |
| Non-PASS inside panel+hotspot regions (7,976 intervals) | 36,862 after split: RefCall 35,825 + NoCall 281 → `dropped_nonpass_no_variant_call` 36,106; GERMLINE 757 → `retained_nonpass_PANEL` 757 (DeepSomatic GT 0/0, ALT-supported) |
| Adapter VCF | **15,193** records: TIER HIGH_CONFIDENCE 9,202 · LOW_CONFIDENCE 5,233 · PANEL 758 (1 PASS + 757 GERMLINE) · HOTSPOT 0; PASS = 10,667 SNV + 3,769 indel; 31 PASS alleles with tumour ALT < 3 (kept) |
| TIER parity vs native SAGE (exact alleles) | 6,459 / 6,459 agree (100 %) |
| PAVE 1.9 output | 14,435 records (chrM PASS record dropped: non-HumanChromosome; **757 GERMLINE records dropped**: oncoanalyser runs PAVE without `-filter_type`, so `PaveConfig.Filter = PASS` and `ChromosomeTask.processVariant` returns on any non-PASS record — the contract's "HOTSPOT/PANEL non-PASS retained" clause (§3, §4.2.5) holds only for `-filter_type PANEL`) |
| PAVE filters on the 14,435 | **PASS 6,980 · PON 2,374 · PONGnomad 841 · PON;PONGnomad 4,240** (PON on 6,614; PONGnomad on 5,081) |
| PURPLE 4.4 somatic VCF | 14,435 records; **PASS 6,980** (6,798 SNV, 182 indel; TIER HIGH_CONF 5,904 / LOW_CONF 1,076); **REPORTED 0**; HOTSPOT 0; somatic fitting variants **4,870** (route C: 7,347 loaded, 3,596 fitting) |
| Overlap with route C PURPLE PASS (6,134) | shared 1,162 · substituted-only 5,818 · route-C-only 4,972 |

Route C for comparison: SAGE PASS 7,347 (HC 5,783 / LC 1,549 / PANEL 15 / HOTSPOT 0; 6,063 SNV/MNV + 1,284 indel) →
PAVE/PURPLE PASS 6,134 (SNV 4,984, MNV 124, indel 1,026), PON 304, PONGnomad 523, PON;PONGnomad 386. SAGE PASS alleles present in the
DeepSomatic PASS set: 1,305 / 7,347.

## 5. Finding: germline leakage in the DeepSomatic PASS set, invisible in DeepSomatic's own NAD

The contract required normal `AD` from a pileup of the normal REDUX BAM rather than DeepSomatic's `NAD`. That pileup (MAPQ ≥ 10,
BQ ≥ 20 for substitutions) shows ALT support in **5,826 of 14,436 PASS alleles**; 4,874 have normal ALT ≥ 3 reads and normal VAF ≥ 0.10
(germline-like), and **3,813 of those carry DeepSomatic `NAD` ALT = 0** (DS normal depth also ~half of the pileup depth). Verified with
`samtools mpileup -q 10 -Q 20` (job 12466856) at six sites, e.g. chr1:1955062 C>A normal 93 reads with ≈44 `A` (DS NAD 54,0; tumour AD
52,35), chr1:4171906 ATTGTT>A normal 96/96 reads carrying the deletion (DS NAD 109,0), chr1:1239125 C>CCTACATGCAT normal ≈26/86 with the
insertion (DS NAD 98,0). Tumour VAF of PASS alleles: 518 < 0.05, 3,263 < 0.10, 4,241 < 0.20, 2,351 < 0.35, 4,063 ≥ 0.35.

Consequences in this run: (i) PURPLE's fit excludes variants with normal ALT reads (`SomaticPurityFitter.java:212-216`), so the pileup column
did the germline screen the contract intended; (ii) PAVE's PON/gnomAD filters removed 7,455 of the 14,435 records, and only **42** PURPLE PASS
variants still have normal ALT ≥ 3 / VAF ≥ 0.10 (35 of them with DS NAD = 0), 6,713 of 6,980 have zero normal ALT reads; (iii) DeepSomatic's
3,769 PASS indels collapse to 182 PURPLE PASS indels (3,900 PASS alleles sit in reference repeats with `REP_C` ≥ 4). The adapter kept all
PASS records as PASS, as the contract requires — the filtering happened natively downstream. Why DeepSomatic 1.10 reports `NAD` ALT = 0 at
sites where the normal BAM plainly carries the allele (REDUX consensus reads, soft-clipping, normal-channel realignment?) is an open question
for the caller review, not something this adapter can or should patch.

## 6. PURPLE / downstream comparison (substituted vs route C)

| Metric | Substituted (DeepSomatic PASS) | Route C (SAGE) |
|---|---|---|
| Purity / ploidy | **0.25 / 2.92** (range 0.21–0.29, ploidy 2.70–4.85) | 0.24 / 2.90 (0.20–0.28, 2.70–3.05) |
| Fit status / method | NORMAL / NORMAL, score 0.5087, somaticPenalty 0.0443, WGD true | NORMAL / NORMAL, score 0.4858, somaticPenalty 0.0308, WGD true |
| QC status | WARN_DELETED_GENES, WARN_HIGH_COPY_NUMBER_NOISE | WARN_DELETED_GENES, WARN_HIGH_COPY_NUMBER_NOISE |
| CN segments / unsupported | 1,307 / 911 | 1,384 / 986 |
| Deleted genes, contamination, AMBER depth, TINC | 649, 0.008, 122, 0 | 649, 0.008, 122, 0 |
| TMB / status | **2.2595 /Mb LOW** | 2.0549 /Mb LOW |
| TML / status | **23 LOW** | 14 LOW |
| MSI indels /Mb, status | **0.0332 MSS** | 0.0560 MSS |
| SV TMB | 90 | 90 |
| Somatic fitting variants | 4,870 | 3,596 |
| Driver catalog (somatic) | ATRX LOH (chrX q21.1, ENST00000373344) | ATRX LOH — identical, no diff |
| REPORTED variants (PURPLE PASS) | 0 | 0 |
| ORANGE | JSON + PDF produced (DNA-only); fit 0.25/2.92, 0 reportable somatic variants | JSON + PDF (DNA+RNA) |
| CHORD | HR_proficient (p_hrd 0.00) | HR_proficient (p_hrd 0.02) |
| SIGS top allocations | Sig12 34 %, Sig1 19 %, misalloc 19 %, Sig3 13 % | Sig3 58 %, Sig9 25 %, misalloc 17 % |
| CUPPA DNA top class | CNS: Glioma 0.77, Mesothelium 0.12 | Prostate 0.64 (DNA); combined with RNA: NET Colorectum/Small intestine 0.45 |

The copy-number/purity layer is essentially unchanged (same AMBER/COBALT/ESVEE inputs); the SNV population differs (1,162 shared PASS
alleles), which moves SIGS and CUPPA-DNA materially and raises TML (23 vs 14) while the driver catalog and reportable set are unchanged
(none). No acceptance threshold is invented here; these are the measured differences.

## 7. Output paths (substituted run)

- PURPLE somatic VCF: `/global/scratch/hpc6049/wigits_runs/KQ0205_full_deepsomatic_substituted_20261113/results/KQ0205_full_deepsomatic_substituted/purple/KQ0205-DNA.purple.somatic.vcf.gz` (+ `.tbi`)
- Driver catalog: `/global/scratch/hpc6049/wigits_runs/KQ0205_full_deepsomatic_substituted_20261113/results/KQ0205_full_deepsomatic_substituted/purple/KQ0205-DNA.purple.driver.catalog.somatic.tsv` (germline: `….driver.catalog.germline.tsv`)
- Purity / QC: `/global/scratch/hpc6049/wigits_runs/KQ0205_full_deepsomatic_substituted_20261113/results/KQ0205_full_deepsomatic_substituted/purple/KQ0205-DNA.purple.purity.tsv`, `/global/scratch/hpc6049/wigits_runs/KQ0205_full_deepsomatic_substituted_20261113/results/KQ0205_full_deepsomatic_substituted/purple/KQ0205-DNA.purple.qc`, `….purple.purity.range.tsv`
- Copy number: `/global/scratch/hpc6049/wigits_runs/KQ0205_full_deepsomatic_substituted_20261113/results/KQ0205_full_deepsomatic_substituted/purple/KQ0205-DNA.purple.cnv.gene.tsv`, `….purple.cnv.somatic.tsv`, `….purple.segment.tsv`
- `/global/scratch/hpc6049/wigits_runs/KQ0205_full_deepsomatic_substituted_20261113/results/KQ0205_full_deepsomatic_substituted/purple/purple.version` (4.4), PURPLE TRACE log: `/global/scratch/hpc6049/wigits_runs/KQ0205_full_deepsomatic_substituted_20261113/work/09/4d2f30*/.command.log`
- ORANGE: `/global/scratch/hpc6049/wigits_runs/KQ0205_full_deepsomatic_substituted_20261113/results/KQ0205_full_deepsomatic_substituted/orange/KQ0205-DNA.orange.json`, `/global/scratch/hpc6049/wigits_runs/KQ0205_full_deepsomatic_substituted_20261113/results/KQ0205_full_deepsomatic_substituted/orange/KQ0205-DNA.orange.pdf`
- Adapter: `/global/scratch/hpc6049/wigits_runs/KQ0205_full_deepsomatic_substituted_20261113/adapter/sage_dir/KQ0205-DNA.sage.somatic.vcf.gz` (+ `.tbi`, `header.txt`), `/global/scratch/hpc6049/wigits_runs/KQ0205_full_deepsomatic_substituted_20261113/adapter/ALLELE_LEDGER.tsv`,
  `/global/scratch/hpc6049/wigits_runs/KQ0205_full_deepsomatic_substituted_20261113/adapter/ledger_summary.json`, intermediates in `/global/scratch/hpc6049/wigits_runs/KQ0205_full_deepsomatic_substituted_20261113/adapter/work/` (pass.norm.vcf.gz, panel.bed, candidates.tsv,
  annotations.tsv, filter_counts.txt, logs)
- Logs: `/global/scratch/hpc6049/wigits_runs/KQ0205_full_deepsomatic_substituted_20261113/logs/adapter.<job>.log`, `/global/scratch/hpc6049/wigits_runs/KQ0205_full_deepsomatic_substituted_20261113/logs/head.<job>.log`, `/global/scratch/hpc6049/wigits_runs/KQ0205_full_deepsomatic_substituted_20261113/.nextflow.log`
- Repo record: `tests/results/deepsomatic-substitution/KQ0205_20261113/` (ledger_summary.json, adapter VCF header, PURPLE purity/qc/drivers, launch.sh)

## 8. Deviations from the task and open issues

1. RNA omitted; NEO and SAGE_VISUALISER excluded; SAGE_APPEND disabled by config (reasons in §1). VirusBreakend reused from route C
   (not in the task's list of existing inputs; 13.7 h saved; SAGE-independent).
2. The ledger lists every PASS record and every non-PASS record inside the panel regions row by row; the 73.7 M RefCall/NoCall/off-panel
   GERMLINE records are represented by FILTER counts only (a 74 M-row TSV was judged useless).
3. Non-PASS HOTSPOT/PANEL retention is implemented per the contract but has no downstream effect under oncoanalyser's PAVE invocation
   (default PASS-only filter) — contract §3/§4.2.5 should be corrected. One PASS record on chrM is dropped by PAVE (non-human contig).
4. `REP_C/REP_S` use an unbounded reference search starting at POS (indel POS+1) — SAGE bounds it by the read-context core; `MH` is
   computed on the reference-derived ALT haplotype. Both are reference-only approximations (README).
5. PURPLE TRACE did not emit `logFilteredFittingCandidate` lines in 4.4; the fitting count (4,870) comes from the DEBUG line
   `somatic fitting variants(…)`.
6. The pass-2 head had to be relaunched once (first launch's 32-CPU tasks never scheduled); `-resume` reused PAVE_SOMATIC. TEAL_PREP was
   still running when this report was written; TEAL is not an ORANGE input. Final Nextflow exit and the execution trace will be appended to
   the run directory by the launcher (`logs/nextflow_head_*.log`, `results/pipeline_info/`).
7. The normal `.bai` is older than the published BAM (publishDir copy order); htsjdk and samtools warn, reads are correct.
8. Open: cause of DeepSomatic `NAD` = 0 at germline sites (§5); whole-genome `RC_INFO`/`RC_CNT` via Codex's evidence export (§1);
   owner decision on whether the PON/gnomAD-filtered population (6,980 PASS) or the raw DeepSomatic PASS set (14,436) is "the DeepSomatic
   population" for SomaReport.
