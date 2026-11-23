# What the native Hartwig metadata actually controls

Reviewed 9 November 2026 for the KQ0205 chr17/chr19 execution test. This
supersedes the earlier implication that samplesheet parser acceptance validated
the downstream meaning of `cancer_type`. A valid ICD-O label is not automatically
a valid expression-cohort key.

## Pinned implementations

| Component | Version / source |
|---|---|
| Oncoanalyser | 3.0.0, `7c74c87a43749952b38c9a18915947570f0595a0` |
| ORANGE | 5.0.1, tag commit `5fb9c98d488eb8f2daabe550b3a1b51f0f7d8cd1` |
| NEO | 1.3, tag `99290c0cadb723d1786ec8026bc7c7c74b446d75` |
| HMF resources | Staged `hmf_pipeline_resources.38_v3.0.0--8` |

The samplesheet passes **one string unchanged** to two different consumers.
No automatic OncoTree, ICD-O, DOID or synonym conversion occurs here.

| Consumer | Meaning of `cancer_type` | Actual accepted vocabulary / behavior |
|---|---|---|
| Oncoanalyser parser | Case-level string in `info` | Recognizes the key and checks syntax; does not validate disease semantics |
| ORANGE | Passed as `-primary_tumor_location` | Display text, stored as a tumour label with an empty DOID; `Carcinoid` is displayable but establishes no formal ontology mapping |
| NEO | Passed as `-cancer_type` | Exact, case-sensitive column name in the supplied transcript-median reference; not an arbitrary histology or OncoTree label |

For NEO, when sample RNA expression is supplied, `NeoScorerTask` chooses
`TpmSource.SAMPLE`. Without sample RNA expression, an exact cohort key selects
`CANCER_TYPE`; an absent/unmatched key selects `COHORT` (pan-cancer). An
unmatched key is not a negative biological result. Cancer-specific reference
TPM remains unavailable for that key even when sample RNA is present.

## The staged NEO reference is the vocabulary authority

Reference:
```text
/global/project/hpcg6049/somatic_pipeline/data/references/oncoanalyser/3.0.0/hmf_pipeline_resources.38_v3.0.0--8/misc/neo/tpm_cohort/isofox.hmf_38_151_2600.transcript_medians.csv.gz
```

Its expression columns, after GeneId/GeneName/TransName, are:

```text
All; Stomach; Mesothelium; Lung; Liver; Nervous system; Thyroid gland;
Bone/Soft tissue; Lymphoid tissue; Small intestine; Bile duct; Kidney;
Esophagus; Ovary; Uterus; Colorectum; Breast; Head and neck;
Gallbladder; Skin; Pancreas; Prostate; Urothelial tract; Other
```

`Carcinoid` is **absent**. `Small intestine` is present. `Small_intestine` is
absent; silently substituting underscores would not fix a cohort mismatch.
`All` is the pan-cancer column, not a specific tumour diagnosis.

There is also a native interface limitation: the samplesheet schema forbids
whitespace in `info`, and the NEO/ORANGE modules interpolate this argument
without shell quoting. Thus a multiword reference term is not reliably carried
by the current native samplesheet route. Do not introduce shell-escaping tricks
or claim an underscored alias works. A later interface change would need its
own explicit mapping and consumer test; none is introduced in this test package.

## KQ0205 handling

- The current sheet retains `cancer_type:Carcinoid` **as ORANGE display text**,
  supported by recorded morphology 8240/3. It is not represented as a NEO
  cohort assignment, a resolved primary site or a confirmed OncoTree/DOID code.
- No specific NEO cohort is assigned. Native sample-expression behavior applies
  when RNA analysis supplies it; native pan-cancer fallback applies otherwise.
- The source has C17.9 alongside the written primary site “Other and ill-defined
  digestive organs”. Preserve both source values; do not resolve their conflict
  by assigning `Small intestine` in the test sheet.
- Disease Ontology nodes exist in the resource file, including a small-intestine
  carcinoid node, but ontology availability alone does not select that diagnosis.

## Other known information and where it belongs

All recorded clinical/laboratory rows and source paths are preserved outside Git:
```text
/global/scratch/hpc6049/hartwig_test_fastqs_20261109/KQ0205/case_metadata.json
```

| Information | Native placement / disposition |
|---|---|
| Case, sample identities and tumour/normal DNA/RNA roles | Dedicated samplesheet columns, fully populated |
| R1/R2 paths | `filepath`, fully populated |
| Combined test-library and lane labels | `info.library_id` / `info.lane`; explicitly technical merged labels |
| Confirmed DNA flowcells | `info.flowcell`; missing RNA flowcell stays unspecified |
| Tumour FFPE, blood normal, frozen RNA, library kits and run dates | Preserved in source metadata; no invented native `info` keys |
| Histology, grade, stage and pathological cellularity | Preserved as clinical observations; not turned into PURPLE constraints or forced results |
| Collection dates | Source precision is month-only; preserve it. ORANGE's direct `-sampling_date` requires `yyMMdd`, and the native CSV exposes no sampling-date key |
| Sex, age, OncoTree and DOID assignment | Not supplied by this bound clinical snapshot; do not invent them |
| DOID ontology configuration | ORANGE has a separate CLI concept; it is not a samplesheet `info` key |

Two ORANGE 5.0.1 limitations matter before treating a native PDF as a case report:

1. `formConfiguredPrimaryTumorDoid` has ontology loading commented out and
   iterates an empty node set. Passing a DOID does not establish successful
   resolution in this version. The native Oncoanalyser module supplies the
   ontology file but does not pass `primary_tumor_doids` from the samplesheet.
2. When sampling date is omitted or malformed, ORANGE defaults to today's
   date. The recorded month cannot safely be expanded to a made-up day.
   A native smoke-test PDF must not be treated as a factually complete patient
   report on the strength of successful rendering.

These are documented interface/tool limitations, not reasons to change the
source clinical records or disable the analysis. This package changes no
interpretation rules or production reporting template.

## Verification evidence

- Native JSON schema and unmodified `Utils.parseInput`: passed, including
  explicit checks of the cancer label and both recorded DNA flowcells.
- The actual staged cohort header was inspected; no new reference download.
- Source inspection used the component tags above, not only current `master`.
- `CheckCancerType.java` exercises the installed NEO JAR on one transcript
  copied from that public reference: exact cohort match, missing `Carcinoid`,
  rejected underscore alias and pan-cancer values retained for an unmatched key.
- `check_metadata_consumers.sbatch` runs this small binary check under Slurm
  and captures the installed ORANGE bytecode and both JAR hashes. Execution
  results are recorded in `EXECUTION.md` after completion.

Completed binary probe **12441043**, exit **0:0**: `Carcinoid` and
`Small_intestine` did not match, `Small intestine` did match, and the unmatched
key retained pan-cancer TPM while returning unavailable cancer-specific TPM.
Installed ORANGE bytecode corroborated the source-level display/DOID finding.
Exact executed JAR hashes and the small harness bind correction are recorded
in `EXECUTION.md`. This does not claim that a full clinical report was validated.

Private evidence:
```text
/global/scratch/hpc6049/hartwig_test_fastqs_20261109/KQ0205/terminology_review/
```

## Primary sources

- [Oncoanalyser input documentation](https://nf-co.re/oncoanalyser/3.0.0/docs/usage/#info-field-keys)
- [Native NEO module](https://github.com/nf-core/oncoanalyser/blob/7c74c87a43749952b38c9a18915947570f0595a0/modules/local/neo/scorer/main.nf)
- [Native ORANGE module](https://github.com/nf-core/oncoanalyser/blob/7c74c87a43749952b38c9a18915947570f0595a0/modules/local/orange/main.nf)
- [NEO exact TPM lookup](https://github.com/hartwigmedical/hmftools/blob/neo-v1.3/neo/src/main/java/com/hartwig/hmftools/neo/score/TpmMediansCache.java)
- [NEO sample/cohort source selection](https://github.com/hartwigmedical/hmftools/blob/neo-v1.3/neo/src/main/java/com/hartwig/hmftools/neo/score/NeoScorerTask.java)
- [ORANGE disease-label handling](https://github.com/hartwigmedical/hmftools/blob/orange-v5.0.1/orange/src/main/java/com/hartwig/hmftools/orange/algo/OrangeAlgo.java)
- [ORANGE configuration and date parsing](https://github.com/hartwigmedical/hmftools/blob/orange-v5.0.1/orange/src/main/java/com/hartwig/hmftools/orange/OrangeConfig.java)
