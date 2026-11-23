# sage-compat adapter (DeepSomatic → oncoanalyser `sage_dir`)

Implements §4 of `docs/DEEPSOMATIC_WIRING_CONTRACT.md`: turns a paired DeepSomatic VCF into a
SAGE-shaped somatic VCF that oncoanalyser 3.0.0 accepts through the samplesheet row
`filetype=sage_dir, sample_type=tumor, sequence_type=dna`, so that SAGE_SOMATIC is skipped and
PAVE 1.9 / PURPLE 4.4 / LINX / ORANGE run natively on the DeepSomatic PASS population.
Only the SAGE_SOMATIC discovery output is replaced; REDUX, AMBER, COBALT, ESVEE, SAGE_GERMLINE and
everything downstream stay native (Option A, two-pass, in the contract).

## Files

| File | Role |
|---|---|
| `SageCompatAnnotator.java` | Single-file Java program run with the SAGE 5.0.2 jar on the class path (`java --class-path sage.jar SageCompatAnnotator.java …`). Reuses SAGE's own classes: `ReferenceData` + `TierSelector` (INFO/TIER), `RepeatInfo.findMaxRepeat` (REP_C/REP_S), `Microhomology.findHomology` (MH); computes TNC from the reference; pileups the normal REDUX BAM with htsjdk for the normal `AD/DP/AF`; writes the bgzipped, tabix-indexed VCF with htsjdk. `--self-test` checks repeat, homology and pileup classification on synthetic data. Derived from Codex's `NativeSageRegions.java` / `ReferenceRepeatAnnotator.java` (dnaseq_nf-purple-adapter-20261109/research/purple_adapter). |
| `sage_compat_adapter.py` | `build` (allele selection + ledger), `merge-bed`, `finalize` (ledger merge, VCF validation, TIER parity vs native SAGE), `purple-summary` (PURPLE outcome counts after pass 2). Standard library only. |
| `run_sage_compat_adapter.sbatch` | Slurm job orchestrating bcftools → panel BED → build → annotate → finalize. `sbatch run_sage_compat_adapter.sbatch configs/<case>.json`. |
| `configs/<run_id>.json` | Inputs for one run (paths only). |
| `tests/test_sage_compat_adapter.py` | Unit test on a synthetic normalised VCF (`python3 tests/test_sage_compat_adapter.py`). |

## Inputs

- DeepSomatic paired-mode VCF (`.vcf.gz` + `.tbi`), single tumour genotype column, called on the REDUX tumour/normal
  BAMs with the HMF `GRCh38_masked_exclusions_alts_hlas.fasta`.
- The REDUX normal BAM (`<normal_id>.redux.bam` + `.bai`).
- HMF resources `hmf_pipeline_resources.38_v3.0.0--8`: `KnownHotspots.somatic.38.vcf.gz`, `DriverGenePanel.38.tsv`,
  `common/ensembl_data`, and the SAGE high-confidence BED oncoanalyser passes to SAGE
  (`HG001_GRCh38_GIAB_highconf_…_noCENorHET7.bed.gz`, `conf/hmf_data.config:119`).
- Containers from the shared cache: `hmftools-sage-5.0.2` (sage.jar + JDK), `bcftools-1.22`.

## Allele population and filter policy

1. `bcftools view -f PASS | bcftools norm -f <ref> -m -any -c e --old-rec-tag DS_ORIGINAL_ALLELE`
   splits multiallelics and left-aligns/trims; every representation change is recorded in the ledger
   (`representation` column) and in `INFO/DS_ORIGINAL_ALLELE`.
2. A split allele is kept only when the DeepSomatic tumour GT calls it (`1` present after bcftools' GT
   projection); uncalled ALTs of multiallelic sites, symbolic ALTs and records without AD are dropped with a reason.
3. Non-PASS records are extracted only inside the SAGE panel+hotspot regions (dumped from `ReferenceData.PanelWithHotspots`
   with `panel-bed`) and kept only when `TierSelector` returns HOTSPOT or PANEL **and** the tumour GT calls the allele
   (in the KQ0205 VCF these are `GERMLINE` calls; `RefCall`/`NoCall` have no called ALT). Their DeepSomatic FILTER
   string is kept as the VCF FILTER (declared in the header) and copied to `INFO/DS_FILTER`.
   Note (verified on KQ0205, 2026-11-13): oncoanalyser invokes PAVE without `-filter_type`, so `PaveConfig.Filter = PASS` and
   `ChromosomeTask.processVariant` drops every non-PASS record; the HOTSPOT/PANEL clause (`ChromosomeTask.java:185-188`) applies
   only under `-filter_type PANEL`. Retaining them keeps the adapter contract-conformant but has no downstream effect.
4. PASS records are never removed for low ALT support; `low_alt_support` (ALT < 3) is a ledger flag only.

## Output (`<out_dir>/sage_dir/<tumor_id>.sage.somatic.vcf.gz` + `.tbi`)

- `##contig` lines and record order follow the reference `.dict`.
- Sample columns: `<normal_id>` then `<tumor_id>`; FORMAT `GT:AD:DP:AF` on both.
  Tumour: DeepSomatic GT/AD/DP/VAF. Normal: GT `./.` (as SAGE writes), AD/DP/AF from the pileup
  (MAPQ ≥ 10 = SAGE `DEFAULT_MIN_MAP_QUALITY`; base quality ≥ 20 = SAGE `MATCHING_BASE_QUALITY` for substitutions;
  unmapped/secondary/supplementary/duplicate/QC-fail reads excluded; indel/complex alleles matched between the aligned
  anchors POS-1 and POS+len(REF)). Thresholds are written to `##sage_compat_adapter_normal_pileup`.
- INFO: `TIER`, `REP_C`/`REP_S` (omitted when no qualifying repeat, as SAGE), `TNC`, `MH` (indels with homology),
  `DS_FILTER` (retained non-PASS), `DS_ORIGINAL_ALLELE`, `DS_NAD`/`DS_NDP` (DeepSomatic's own normal counts, informational).
- QUAL = DeepSomatic QUAL. Header: `##source=deepsomatic-<ver>+sage-compat-adapter-<ver>` and a
  `##sage_compat_adapter_*` block (versions, resource SHA256s, source VCF SHA256, pileup thresholds, Slurm job id).
- **Not emitted, by design (contract §4.2.7):** `LPS`, `RC_INFO`, `RC_CNT`, `RC_QUAL`, `RC_JIT`, `UMI_CNT`, quality FORMAT
  fields, `PON_*`, `GND_FREQ`, `MQF`, `TQP`, `TINC_RECOVERED`. Consequences: no phased MNV/inframe reinterpretation,
  SAGE_APPEND (RNA) and NEO cannot consume the file (omit RNA rows / `--processes_exclude neo`).

### Known approximation
`REP_C/REP_S` are searched on the reference only, starting at POS (indel: POS+1), rightwards and unbounded;
SAGE searches the same way but inside the read-context core (`VariantReadContextBuilder.setMaxRefRepeat`), so a few
counts can differ at repeat boundaries (Codex observed 3/15 k). `MH` is computed on the ALT haplotype built from the
reference, which equals SAGE's result for a read that carries the allele without errors.

## Ledger

`ALLELE_LEDGER.tsv` has one row per allele that was considered: every DeepSomatic PASS record (after splitting) and
every non-PASS record inside the panel regions. Columns: `record_id` (original CHR|POS|REF|ALTs|used_alt_idx), allele,
`source_class`, DeepSomatic FILTER/QUAL/GT/AD/DP/VAF/NAD/NDP, `representation`
(`unchanged`, `multiallelic_split`, `realigned_or_trimmed`), `low_alt_support`, `disposition`
(`candidate_pass`, `candidate_nonpass`, `dropped_alt_not_in_called_gt`, `dropped_symbolic_or_missing_alt`,
`dropped_missing_tumour_ad`, `dropped_duplicate_allele`), the adapter annotations (TIER, REP_C, REP_S, TNC, MH, normal
pileup) and `final_disposition` (`retained_pass`, `retained_nonpass_<TIER>`, `dropped_nonpass_tier_<TIER>`, or the
build disposition). The bulk classes (RefCall, NoCall, GERMLINE outside the panel) are counted in
`ledger_summary.json["deepsomatic_records_by_filter"]` rather than listed row by row.

`ledger_summary.json`: FILTER counts of the whole DeepSomatic VCF, dispositions, retained tiers, VCF validation result,
TIER parity against the native SAGE VCF on the exact-allele intersection, and (after `purple-summary`) PURPLE PASS /
REPORTED / PON-filtered counts and the driver catalog.

## Observed on KQ0205 (2026-11-13)

See `docs/DEEPSOMATIC_SUBSTITUTION_KQ0205_20261113.md`: 14,436 PASS → 15,193 adapter records → PAVE 6,980 PASS / 7,455 PON-or-gnomAD →
PURPLE purity 0.25 / ploidy 2.92, 0 REPORTED. The normal pileup exposed ALT support in 5,826 PASS alleles (3,813 with DeepSomatic NAD = 0).

## Running pass 2

See `docs/DEEPSOMATIC_SUBSTITUTION_KQ0205_<date>.md` for the samplesheet pattern: `redux_dir` rows (tumour, normal),
`sage_dir` row (tumour → `<out_dir>/sage_dir`), existing `amber_dir`, `cobalt_dir`, `esvee_dir`, `virusinterpreter_dir`
rows, no RNA rows, `--processes_exclude neo`.
