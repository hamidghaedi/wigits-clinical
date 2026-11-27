# Oncoanalyser → pVACtools interface specification

**Status: canonical for the neoantigen layer (gap G28).** Every field the layer consumes is listed here with
its source file, its meaning and the conversion rule. Nothing reads an oncoanalyser output that is not in this
document. Changing a mapping means changing this file and `tests/test_neo_interface.py` in the same commit.

Verified against **oncoanalyser 3.0.0** outputs of `KQ0004_full_native_20261114_onco3.0.0`
(`results/KQ0004/`) and the nf-core test run `subject_a`, read on 2026-11-23. Tool versions as recorded in
each run's `pipeline_info/nf_core_oncoanalyser_software_mqc_versions.yml`: SAGE 5.0.2, PAVE 1.9, PURPLE 4.4,
LINX 2.3.1, LILAC 2.0, ISOFOX 2.0.1, NEO 1.3, ORANGE 5.0.1.

Paths below are relative to `RESULTS=<run_root>/results/<case>`; `<T>` is the tumour DNA sample id
(e.g. `KQ0004-2023-AAAD-2-DNA`), `<N>` the normal (`…-gDNA`), `<R>` the tumour RNA (`…-2-RNA`).

## 0. Sample and transcript conventions

| Item | Value observed | Rule |
|---|---|---|
| PURPLE somatic VCF sample columns | `<N>`, `<T>` (2 columns, normal first) | index by sample name from the `#CHROM` line, never by position |
| SAGE_APPEND directories | `sage_append/<T>/`, `sage_append/<N>/` | the RNA genotype lives in the **`<T>` directory's** VCF, whose `#CHROM` line carries `<N>`, `<T>`, `<R>` |
| Transcript annotation | HMF Ensembl data under `data/references/oncoanalyser/3.0.0/hmf_pipeline_resources.38_v3.0.0--8/common/ensembl_data/` (`ensembl_gene_data.csv`, `ensembl_trans_exon_data.csv`, …) | VEP must run with a cache of the **same Ensembl release** as this bundle; transcript ids are versionless `ENST…` in PAVE/NEO/ISOFOX, so compare on the unversioned id |
| Reference genome | `GRCh38_masked_exclusions_alts_hlas.fasta` (see `docs/directory_layout.md` §3) | pVACseq, VEP and IGV all use this FASTA; never substitute GRCh38.d1.vd1 |

## 1. Somatic short variants — `purple/<T>.purple.somatic.vcf.gz`

PURPLE republishes the SAGE call set with PAVE annotation and purity-adjusted fields. This is the pVACseq
input, not `sage/`.

| Field | Source | Consumed as | Note |
|---|---|---|---|
| `FILTER` | VCF | eligibility: `PASS` only | |
| `INFO/TIER` | SAGE | eligibility / evidence (`HOTSPOT`, `PANEL`, `HIGH_CONFIDENCE`, `LOW_CONFIDENCE`) | |
| `INFO/IMPACT` | PAVE, 10 fields `[Gene, Transcript, CanonicalEffect, CanonicalCodingEffect, SpliceRegion, HgvsCoding, HgvsProtein, OtherReportableEffects, WorstCodingEffect, GenesAffected]` | pre-filter to coding consequence classes | canonical transcript only |
| `INFO/PAVE_TI` | PAVE, `Gene\|GeneName\|Transcript\|Effects\|SpliceRegion\|HGVS.c\|HGVS.p\|RefSeqId\|Exon\|Codon`, comma-separated per transcript | cross-check VEP's transcript set | the non-canonical transcripts PAVE kept |
| `INFO/PURPLE_VCN`, `PURPLE_CN`, `PURPLE_MACN`, `PURPLE_AF` | PURPLE | clonality features; CCF = `PURPLE_VCN / PURPLE_CN` reported alongside `SUBCL` | |
| `INFO/SUBCL` | PURPLE | subclonal likelihood (0 = clonal) | absent when zero |
| `INFO/BIALLELIC`, `BIALLELIC_PROB` | PURPLE | evidence | |
| `INFO/LPS`, `LPS_RC` | SAGE local phase set + read count | proximal-variant grouping (§5) | read-backed, within fragment length |
| `INFO/HOTSPOT`, `NEAR_HOTSPOT`, `REPORTED`, `GND_FREQ`, `MAPPABILITY`, `TNC`, `REP_C/REP_S`, `NEARBY_INDEL` | SAGE/PAVE | artefact and context flags (FFPE screen uses `TNC` + `MUC`) | |
| `FORMAT/AD`, `DP`, `AF` (sample `<T>`) | SAGE | pVACseq `--tdna-*`: `AD` → tumour DNA alt depth, `DP` → depth, `AF` → VAF | pVACseq takes `AD/DP/AF`; use the vcf-readcount-annotator form of the names |
| `FORMAT/AD`, `DP` (sample `<N>`) | SAGE | normal support eligibility (alt ≥ 1 in normal → exclude unless hotspot) | |
| `FORMAT/MUC` | SAGE min unique fragment coordinates in alt reads | FFPE/duplication screen | low `MUC` + C>T/G>A (`TNC`) → Hold |
| `FORMAT/FSB`, `RSB` | SAGE strand bias (ref,alt) | artefact screen | |
| `FORMAT/RC_CNT` | `[Full, PartialCore, Core, Realigned, Reference, Total]` | evidence detail; the field NEO requires (see G23) | |

## 2. RNA support — `sage_append/<T>/<T>.sage.append.vcf.gz`

Same records as §1, re-genotyped on the tumour RNA BAM; the `<R>` column is appended.

| Field (sample `<R>`) | Consumed as (pVACseq) | Rule |
|---|---|---|
| `AD[1]` | `RAD` tumour RNA alt depth | |
| `DP` | `RDP` tumour RNA depth | |
| `AF` | `RAF` tumour RNA VAF | |
| `AD`, `DP` combined | `rna_coverage_status` ∈ `SUPPORTED` (DP ≥ 10 and alt ≥ 2), `COVERED_NO_ALT` (DP ≥ 10, alt < 2), `NOT_COVERED` (DP < 10) | **an absent mutant allele at low depth is not evidence of non-expression**; the three states must stay distinguishable downstream |
| `FORMAT/LPSA` | RNA-side local phase set | corroborates §5 phasing |

Observed scale (whole genome, all records incl. non-coding): 112–472 records at RNA `DP ≥ 10` per KQ case,
1,227 in KQ0205. Coding-candidate coverage is reported per case in the test record.

## 3. HLA typing — `lilac/<T>.lilac.tsv`, `lilac/<T>.lilac.qc.tsv`

| Field | Consumed as | Rule |
|---|---|---|
| `Allele` | pVACseq `--allele` | `A*02:01` → `HLA-A*02:01`; **de-duplicate** (homozygous alleles appear twice); class I = `A/B/C`, class II = `DPA1/DPB1/DQA1/DQB1` |
| `TumorCopyNumber` | allele retention | < 0.5 → treat the allele's candidates as presentation-compromised |
| `SomaticMissense`, `SomaticNonsenseOrFrameshift`, `SomaticSplice`, `SomaticInframeIndel` | allele disruption | any > 0 → flag |
| `RefTotal/RefUnique`, `TumorTotal`, `RnaTotal` | typing confidence / allele expression | |
| qc `Status` | gate | `PASS`, `WARN_*` (proceed, record), `FAIL` (no class I predictions for that case — subject_a is `FAIL`) |
| qc `ScoreMargin`, `NextSolutionAlleles`, `HlaYAllele`, `SomaticVariantsMatched` | typing confidence | |

LILAC 2.0 does **not** report `DRB1/3/4/5`; DR typing needs an external typer (phase 2).

## 4. Expression — `isofox/<T>.isf.gene_data.tsv`, `isofox/<T>.isf.transcript_data.tsv`

| Field | Consumed as | Rule |
|---|---|---|
| gene `GeneId`, `GeneName`, `AdjTPM` | pVACseq `GX` (gene expression) | join on `GeneId` (ENSG); `AdjTPM` is the fitted value, `RawTPM` is unfitted |
| gene `MedianTpmCancer`, `PercentileCancer`, `MedianTpmCohort`, `PercentileCohort` | cohort context for review | `-1` means "no cohort value" |
| transcript `TransName`, `AdjTPM` | pVACseq `TX` (transcript expression) | join on `TransName` (ENST, versionless) |
| transcript `SJSupported`, `UniqueSJFragments` | transcript-structure confidence | |

## 5. Phasing — SAGE phase sets (+ germline)

The maintained adapter is `workflow/neo_phasing.py` (`sage-lps-to-pvac-hp-v1`).
It copies all selected somatic anchors into the single-tumour proximal VCF,
retaining original LPS/LPS_RC and genotypes. Neighbourhood is ±100bp. The main
selected VCF and its measured DNA/RNA/expression values are unchanged.

Heterozygous phase comes only from exact-allele matches in the native **somatic**
SAGE tumour namespace (`sage/somatic/*.sage.somatic.vcf.gz`), falling back explicitly
to PURPLE somatic when the raw file is absent. One LPS with positive LPS_RC and
positive measured tumour alternate AD and PURPLE tumour heterozygosity becomes an opaque, chromosome-scoped FORMAT/HP
token. Multiple/overlapping LPS, duplicate records and missing support stay
unresolved. Unknown HP is omitted from that record's FORMAT, not encoded as
HP='.' (pVACtools can compare two missing values as equal). Source HP is retained
as SW_ORIGINAL_HP, not trusted across caller namespaces. Raw SAGE may carry GT=./.; its positive AD/LPS_RC supports phase while the
PURPLE genotype is retained. A contradictory called source genotype is unresolved.
No GT phasing, clean
normal counts, long-range phase or RNA phase is inferred.

Germline LPS identifiers belong to a separate call and **cannot** be joined by
integer equality. A nearby PASS germline allele gets HP only if that exact allele
has supported phase in the tumour somatic source. Actual tumour homozygous-alt
handling remains pVACtools' native genotype rule. A missing tumour column never
uses a normal genotype as its substitute. PASS filters are retained. Nonselected
somatic context is recorded in the evidence ledger but not exported: pVACtools
classifies context absent from the main VCF as germline and would incorrectly
alter WT peptides. Distinct germline alleles at an anchor position are also
recorded as unsupported context: pVACtools tests position rather than exact
allele when assigning context origin. No selected candidate is removed.

The adapter writes `<case>.phasing_evidence.json` with per-allele source/support,
derived HP/status, omitted context and limitations; it is referenced in the
interface manifest. After prediction, `phasing_usage.json` records actual pVAC
converter-table rows separately from available phase evidence. Neither technical
completion nor an HP tag proves incorporation into every peptide.

Pinned [SAGE 5.0.2 phasing documentation](https://github.com/hartwigmedical/hmftools/blob/sage-v5.0.2/sage/README.md#7-phasing)
describes local read evidence. [pVACtools 5.4 proximal parser](https://github.com/griffithlab/pVACtools/blob/v5.4.0/pvactools/lib/proximal_variant.py)
requires anchors and applies same-transcript **missense** proximal consequences;
non-missense context remains an explicit library limitation. The
[input converter](https://github.com/griffithlab/pVACtools/blob/v5.4.0/pvactools/lib/input_file_converter.py)
determines somatic/germline context by presence in the main VCF. Synthetic
container tests exercise these actual methods and the resulting protein sequence.

## 6. Native predictions to compare against — `neo/`

| File | Key fields | Use |
|---|---|---|
| `neo/finder/<T>.neo.neo_data.tsv` | `NeId`, `VariantType` (`MISSENSE`, `FRAMESHIFT`, `INFRAME_DELETION`, `INFRAME_FUSION`, `OUT_OF_FRAME_FUSION`), `VariantInfo` (`chr:pos:ref:alt` for point variants, `chrUp:posUp:orientUp;chrDown:…` for fusions), `VariantCopyNumber`, `CopyNumber`, `SubclonalLikelihood`, `UpstreamAA`/`NovelAA`/`DownstreamAA`, `WildtypeAA`, `NmdMin/Max`, `UpTranscripts`/`DownTranscripts` | candidate-level join key; `VariantInfo` joins to §1 on `chr:pos:ref:alt` |
| `neo/scorer/<T>.neo.peptide_scores.tsv` | `NeIds`, `Allele` (`C1203` form), `Peptide`, `FlankUp/FlankDown`, `Score`, `Rank`, `Likelihood`, `LikelihoodRank`, `ExpLikelihood`, `ExpLikelihoodRank`, `EffectiveTpm`, `RecogSim`, `OtherAlleleRecogSim`, `AlleleCN`, `AlleleDisrupted` | peptide × allele comparison with pVACseq |
| `neo/scorer/<T>.neo.neoepitope.tsv` | per-neoepitope `RnaFrags`, `RnaDepth`, `TpmUp/Down`, `ExpectedTpm`, `RawEffectiveTpm`, `EffectiveTpm`, `TpmCancerUp/Down`, `PeptideCount`, `TpmSource` | RNA/expression evidence already computed by NEO |
| `neo/annotated_fusions/<T>.isf.neoepitope.tsv` | ISOFOX fusion support for NEO fusion candidates | fusion RNA evidence |

Allele-name normalisation (the join key across NEO, LILAC and pVACseq):
`C1203` ⇄ `C*12:03` ⇄ `HLA-C*12:03` — insert `*` after the gene, `:` before the last two digits.

## 7. Fusion neoepitopes (DNA) — `linx/somatic_annotations/<T>.linx.neoepitope.tsv`

`geneNameUp/Down`, `chromosomeUp/Down`, `positionUp/Down`, `orientationUp/Down`, `junctionCopyNumber`,
`copyNumber`, `insertSeq`, `chainLength`, `transcriptsUp/Down`. Phase 1 keeps these as **NEO-scored** candidates
(no pVACfuse); they join to §6 on the composite `VariantInfo`.

## 8. Case-level context

| Source | Fields | Use |
|---|---|---|
| `purple/<T>.purple.purity.tsv` | `purity`, `ploidy`, `status`, `wholeGenomeDuplication`, `msStatus`, `msIndelsPerMb`, `tml`, `tmlStatus`, `tmbPerMb`, `tmbStatus`, `svTumorMutationalBurden`, `runMode` | CCF maths, case header, candidate-count context |
| `purple/<T>.purple.qc` | QC status | gate (the purity TSV has **no** `qcStatus` column — read the `.purple.qc` file) |
| `orange/<T>.orange.json` → `immuneEscape` | `hasHlaEscape`, `hasAntigenPresentationPathwayEscape`, `hasIFNGammaPathwayEscape`, `hasCD58InactivationEscape`, `hasPDL1OverexpressionEscape`, `hasEpigeneticSETDB1Escape` | case-level flags on every candidate (gap G10 renders the same table in the reporter) |
| `orange/<T>.orange.json` → `isofox.novelSpliceJunctions`, `isofox.fusions` | RNA-only event inventory | phase 2 input |
| `pipeline_info/nf_core_oncoanalyser_software_mqc_versions.yml` | tool versions | provenance in every candidate table |

## 9. Outputs of the adapter

`workflow/neo_export.py` writes, per case, into `RUN_ROOT/<run_id>/neo_ext/`:

| File | Content |
|---|---|
| `<case>.pvacseq_input.vcf.gz` | VEP-annotated somatic PASS variants with `AD/DP/AF`, `RAD/RDP/RAF`, `GX`, `TX` |
| `<case>.proximal.vcf` | §5 anchors/context before VEP and indexing |
| `<case>.phasing_evidence.json` | Original phase/read/genotype evidence, derived HP, unresolved/omitted context |
| `phasing_usage.json` (prediction directory) | Actual pVAC converter uptake, distinct from available source evidence |
| `<case>.hla_alleles.txt` | de-duplicated `HLA-…` class I (and class II) allele list |
| `<case>.interface_manifest.json` | source file paths + sha256, tool versions, record counts, rule version of this document |

The manifest is what makes a candidate table auditable: every number in it must be traceable to a file listed
in the manifest.
