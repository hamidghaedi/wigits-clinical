# Neoantigen discovery on top of oncoanalyser — proposal (draft for discussion)

Date: 2026-11-23. Status: adopted as the working plan on branch `neo/discovery-layer` (gap G28); phases below are being executed in order. Written against the
ChatGPT review thread of 2026-11-05/11 (INTerpath-001 → bioinformatics workflow → published pipelines →
oncoanalyser integration), the oncoanalyser 3.0.0 smoke-test outputs on Frontenac
(`test_20261014_onco3.0.0`, subject_a), and the current wigits-clinical repository (commit 7a8f207).

Goal as stated in the thread: a reproducible, expert-reviewed neoantigen shortlist per tumour/normal WGS +
tumour RNA case (WGTS → HLA typing → annotation → RNA support + phasing → prediction and ranking →
pVACtools/pVACview review), as a research output feeding downstream validation — not a manufacturing-ready
vaccine design.

## 1. What changes when you look at our actual outputs

The thread's oncoanalyser output table was written from the docs. The smoke-test results show the native
layer already covers more of the "evidence matrix" than it assumed:

| Evidence the thread wants | Where it already is (files from `results/subject_a/`) | Gap remaining |
|---|---|---|
| Candidate peptides, per HLA allele | `neo/scorer/*.neo.peptide_scores.tsv`: `Peptide`, `FlankUp/Down`, `Allele`, `Score`, `Rank`, `Likelihood`, `LikelihoodRank`, `ExpLikelihood` (expression-weighted), `RecogSim` (TCR-recognition similarity), `OtherAlleleRecogSim` | Class I only, 8–12-mers; Hartwig's own presentation model, not NetMHCpan/MHCflurry |
| Clonality / copy number per candidate | `neo/finder/*.neo_data.tsv`: `VariantCopyNumber`, `CopyNumber`, `SubclonalLikelihood`; `NmdMin/Max` (nonsense-mediated decay), `UpstreamAA/NovelAA/DownstreamAA`, `WildtypeAA`, transcript lists | CCF not stated explicitly; derive from PURPLE purity + VCN |
| HLA allele retained in tumour? | `peptide_scores.tsv`: `AlleleCN`, `AlleleDisrupted`; `lilac/*.lilac.tsv`: per-allele `TumorCopyNumber`, `SomaticMissense/Nonsense/Splice/...`, `RnaTotal` | Class II alleles absent (LILAC is class I) |
| RNA support per somatic variant | `sage_append/<tumor>/*.sage.append.vcf.gz` FORMAT `AD`, `DP`, `AF`, `RC_CNT`, `MUC` (unique fragments), strand bias `FSB/RSB` — genotyped on the RNA BAM | Mapping to pVACseq `RAD/RDP/RAF` is a format conversion only |
| Phasing of proximal variants | SAGE `LPS`/`LPSA` local phase sets (read-backed, within read length); NEO consumes them | Germline–somatic phasing and variants farther apart than a read need WhatsHap-style phasing |
| Fusion neoepitopes (DNA) | `linx/somatic_annotations/*.linx.neoepitope.tsv` (junction CN, transcripts) → NEO scores them with frame (`OUT_OF_FRAME_FUSION` seen in the test) | — |
| Fusion/splice (RNA-only) | `isofox/*.isf.pass_fusions.tsv`, `*.alt_splice_junc.tsv`, `*.retained_intron.csv` as evidence; `neo/annotated_fusions/*.isf.neoepitope.tsv` annotates NEO candidates with RNA fusion support | Peptide reconstruction for RNA-only events (pVACfuse/pVACsplice) — second phase |
| Expression | `isofox/*.isf.gene_data.tsv`, `*.transcript_data.tsv` (TPM); NEO `EffectiveTpm` | Map to pVACseq `GX`/`TX` |
| Immune escape summary | ORANGE JSON immune-escape table (G10 in the gap register: computed, not rendered) | Render it; feed as a case-level gate |

Consequence: Phase 1 is less "build an adapter to recover evidence" and more "export what NEO/LILAC/
SAGE_APPEND already computed, add an independent predictor, and put a review process on top".

## 2. Two project-level decisions that come before any code

**A. NEO conflicts with the DeepSomatic-substitution workstream.** Our current substituted route runs
`--processes_exclude neo` because NEO parses SAGE's six-element `RC_CNT` and throws when it is absent
(G23; `DEEPSOMATIC_WIRING_CONTRACT.md`). A neoantigen workflow needs NEO (or a replacement) on every
case. Options: (i) vaccine-track cases run the native SAGE path end to end; (ii) finish the authentic
`RC_CNT`/context recovery G23 asks for so NEO runs on substituted calls; (iii) treat pVACseq as the
primary predictor and NEO as optional corroboration. (i) is cheapest and keeps the two tracks honest
about what they are; (iii) loses the Hartwig likelihood/recognition features. Recommend (i) now, (ii) later.

**B. Reference genome.** oncoanalyser runs on `GRCh38_masked_exclusions_alts_hlas.fasta`; our archive BAMs
are GRCh38.d1.vd1 (G17, `NATIVE_BAM_ENTRY.md`). LILAC in particular depends on the HLA-region handling
of the HMF genome. For the vaccine track, either realign from FASTQ on the HMF genome (what the smoke
test did) or validate LILAC/SAGE behaviour on d1.vd1 BAMs before trusting HLA calls. IGV review and
pVACseq must use the same FASTA as the calls.

## 3. Ideas, ordered by how much each buys for how little

1. **Interface specification first (the thread's own recommendation), as a repo document.** One table per
   oncoanalyser source file: columns consumed, unit, conversion rule to the pVACseq/pVACview field, test
   that proves the mapping on subject_a. The §1 table above is its skeleton. This is also the only way to
   keep the extension decoupled from oncoanalyser upgrades (3.0.0 → 3.x changes NEO columns).

2. **An `oncoanalyser → pVACseq` adapter as a small Nextflow subworkflow under `workflow/`**, consuming the
   published `results/<case>/` directory (not the work dir), emitting: VEP-annotated somatic VCF with
   tumour/normal `AD/DP/AF` and RNA `RAD/RDP/RAF`, `GX/TX` expression TSVs from ISOFOX, HLA allele string
   from LILAC (`A*03:01` → `HLA-A*03:01`), a proximal-variant VCF from SAGE phase sets + germline calls.
   We already hold the annotation stack: SomaSieve's VEP 115.2 container, offline cache, FASTA, REVEL and
   AlphaMissense plugins (all on `cancer_dbs/vep`, `insilico`). Transcript set must be pinned to what
   NEO/PAVE used (HMF Ensembl data in `hmf_pipeline_resources.38_v3.0.0--8/common/ensembl_data`) so the
   two predictors talk about the same isoforms.

3. **Run pVACseq class I as the second, independent presentation predictor; compare, don't replace.**
   NEO's `Likelihood`/`RecogSim` versus NetMHCpan-EL/MHCflurry percentile ranks on the same peptides gives a
   concordance table per case; disagreements are the review queue. pVACtools must run locally (IEDB
   tools + MHCflurry in the griffithlab container; NetMHCpan is academic-licence software — record the
   licence in `VERSIONS.md`). No patient data to the IEDB web API.

4. **Unified candidate table, tiered not scored.** Join NEO + pVACseq + LILAC + SAGE_APPEND + PURPLE into
   one row per (variant, transcript, peptide, allele) with the thread's evidence categories, plus tier
   (Prioritise / Review / Hold / Exclude) and explicit reasons. Eligibility filters (PASS, normal support,
   allele retained, RNA covered) are separate from ranking. Keep the full table; publish the shortlist.
   Schema versioned in the repo (`docs/evidence/` already exists for the clinical layer — same pattern).

5. **Review stage = pVACview + IGV + a written decision per candidate**, adopting the ImmunoNX manual-review
   README as the SOP. Review inputs: the pVACseq aggregated TSV + metrics JSON (native pVACview format),
   plus our unified table as a side panel. Decisions go to `wigits_releases/<run_id>/neo_review/` with
   reviewer, timestamp, rationale — the same release discipline as the clinical reports.

6. **Validation path we can start immediately: HCC1395.** ImmunoNX publishes an HCC1395 worked example
   (322 → 78 candidates). We already have SEQC2 HCC1395 FASTQs on Frontenac (`data/input/SRR7890960`,
   `SRR7890964`, `seqc2_wgs_ffpe_manifest.tsv`) and a chr22 SEQC2 test fixture in dnaseq_nf. Running
   oncoanalyser WGTS on HCC1395 (RNA from SEQC2 is public) and comparing our shortlist with ImmunoNX's
   is a concrete, external benchmark for Phase 1. It also exercises the FFPE arm of the SEQC2 set.

7. **FFPE is our differentiator and our risk.** INTerpath-001 excluded ~17% of screened patients for
   tissue inadequacy; the MOHCC bladder cohort is largely FFPE. Feed the existing FFPE evidence into
   eligibility: DeepSomatic FFPE model calls as a second opinion on candidate variants, SAGE strand-bias
   and `MUC` (unique fragments) thresholds, the trim-divergence damage signal from dnaseq_nf. A candidate
   whose variant is C>T/G>A with low unique-fragment support is "Hold" regardless of binding score.

8. **Case-level gates from the immune-escape table (G10).** B2M/HLA LOH, JAK1/2, antigen-presentation
   pathway defects: render ORANGE's table in our reporter and use it to annotate the whole shortlist
   ("class I presentation compromised") rather than only per-allele `AlleleDisrupted`.

9. **Phase 2 items, in order of value:** HLA class II typing (HLA-HD or arcasHLA on normal WGS/RNA) +
   NetMHCIIpan; germline–somatic proximal phasing beyond SAGE's read-length window (WhatsHap on tumour
   BAM, pVAC proximal VCF); RNA-only fusions (ISOFOX pass_fusions → pVACfuse) and splice (alt_splice_junc →
   pVACsplice); CCF from PURPLE purity/VCN; diversity-aware selection (distinct mutations × alleles × class).

10. **Bladder-specific framing.** INTerpath-005 (MIBC) is the relevant trial; TMB in MIBC is intermediate,
    APOBEC-driven. Ideas: report APOBEC-signature share (SIGS output already exists) alongside candidate
    counts; prioritise clonal truncal candidates for a neoadjuvant setting; keep the class-I-retained
    fraction of candidates as a case-level number. These are research hypotheses, not selection rules.

## 4. What this is not

Not a clinical deliverable. The clinical layer (SERVE/PROTECT/ROSE/reporter, G01–G14) is about
actionability of drivers; this track is about passenger-mutation immunogenicity. Keep them as separate
subworkflows with separate release folders and separate gap-register rows, sharing only the oncoanalyser
results directory, the reference contract and the review discipline.

## 5. Proposed first steps (if adopted)

1. Gap-register row G28 "Neoantigen discovery and review extension" pointing to this document.
2. `docs/neoantigen_interface.md` — the §1 table completed against subject_a files, with a pytest-style
   check per mapping.
3. Decision on §2A (native SAGE path for vaccine-track cases) recorded in AGENTS/CONTRIBUTING.
4. HCC1395 WGTS oncoanalyser run (one Slurm head job; references now durable in
   `data/references/oncoanalyser/3.0.0/`).
5. pVACtools container into the shared SIF cache; IEDB/MHCflurry data under `data/knowledgebases/`.
6. Adapter subworkflow + unified table on subject_a and HCC1395; compare with ImmunoNX's 78.

## References named in the thread (to verify before citing in repo docs)

ImmunoNX (Griffith lab, WUSTL; protocol repo `griffithlab/ImmunoNX_protocol`, ManualReview_README.md);
pVACtools docs (pvactools.readthedocs.io; pVACview *Genome Medicine* 2024); nextNEOpi (`icbi-lab/nextNEOpi`,
*Bioinformatics* 2022; DSL1); nf-neoant (`bioinfo-pf-curie/nf-neoant`); nf-core/epitopeprediction;
NeoDisc (Huber et al., *Nat Biotechnol* 2024/25; restricted licence); NeoQual-nf (MSKCC); hmftools
`neo/README.md`, `lilac/README.md`; INTerpath-001 ESMO 2026 LBA1.
