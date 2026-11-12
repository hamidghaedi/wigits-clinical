# DeepSomatic allele-preserving native evidence adapter

Development on `hartwig_adoptation`, 11 November 2026 (execution timestamps
extend into 12 November UTC). This is the context increment for
[G23](gap_register.md), following the measured-count/reference annotation
adapter. Production/default callers and executing jobs remain unchanged.

## Why the original 115 alleles were missing

The regional DeepSomatic pilot has 183 PASS alleles. Native SAGE 5.0.2 output
matched 68. Its `specific_var_file` restricts discovery; it is not a force-call
interface. `ReadContextCounters.filterCandidates` removes candidates before
normal evidence collection, and `RegionTask.finaliseResults` filters again.
Disabling only the final hard filter does not bypass the earlier exclusion.

An isolated evidence-export overlay now preserves discovered requested
candidates across those two gates while recording their real native filter
outcomes. It does not change thresholds, reads, qualities, native counting,
duplicate handling, or DeepSomatic's decisions. Native SAGE remains available
as the control. The overlay is never installed over the shared image/JAR.

The measured result is 180/183 allele matches: 176 unique contexts and four
with two contexts. This recovers 112 of the originally missing 115. The four
multiple-context records use the first native CandidateWriter context: native
`AltContext.selectCandidates` and `RefContextCache.altCandidates` order by
full matches, core matches and event count. Other contexts remain in the raw
evidence file; their overlapping counts are not added.

Subsequent observed exact-CIGAR contexts supply the remaining three indels.
Thus all183 PASS alleles have authentic contexts, including all115 initially
unmatched alleles. Context availability is distinct from fragment support,
caller agreement and complete workflow equivalence; see the count discrepancies
below. The accepted master remains the complete15281-row candidate derivative.

## Inputs, outputs and execution

1. `workflow/build_sage_evidence.sbatch` builds the two-class overlay using
   the pinned sources in `workflow/patches/sage-5.0.2-evidence-export.sources.json`
   and the adjacent patch. Retain the upstream GPL licence with the source
   package. It records original and overlay JAR SHA256s.
2. `workflow/probe_sage_context.py --evidence-jar <jar>` measures requested
   pilot alleles with native tumour/normal REDUX BAMs, actual BQR/jitter/MS
   tables, reference and resource bundle. The mode requires explicit alleles
   and forbids restricted-scope TINC. Its current bounded scope is at most
   1000 PASS alleles; it is not the whole-genome workflow entry.
3. `workflow/merge_sage_evidence.py --source <annotated-deepsomatic.vcf.gz>
   --evidence <native-evidence.vcf.gz> --candidates <native-candidates.tsv.gz>
   --out <complete.vcf>` preserves the complete candidate derivative and adds
   genuine evidence. No caller intersection is used.
4. Discovery-missing indels may use `InspectMissingContexts.java`,
   `prepare_context_seeds.py`, `EnrichSeedContext.java` and native SAGE_APPEND
   to measure an actual observed exact-CIGAR context. The optional merge
   `--supplement` accepts these separately measured records only where native
   discovery supplied no record. It refuses replacement and a CIGAR that does
   not encode the exact allele.
   For multiple observed contexts, `--primary` consumes the actual native
   `AltContext.addReadContext/selectCandidates` result from the same reads.
   Every chosen context must still have an observed exact-CIGAR record.
   Native quality, full/core matching and event-count ordering are retained.
5. `workflow/test_evidence_consumers.py` exercises native SAGE context
   reconstruction, RNA append, NEO's count reader and PAVE. The explicit
   measured-PASS test view is separate from the complete preserved VCF. Its
   execution record lists every unresolved PASS allele. A successful subset
   test is not full downstream readiness.

Caller QUAL, FILTER, INFO and existing FORMAT values stay unchanged. Colliding
native fields are retained as `SAGE_*`; `SAGE_FILTER` does not replace FILTER.
Native `RC_CNT` is measured, never constructed from AD. Missing context is
explicit, never a reference genotype or evidence of absent biological support.

NEO sums the first four native `RC_CNT` categories and reads DP. Its input
projection therefore uses the matching measured `SAGE_DP`, retaining the
original caller depth as `DS_DP`. RNA keeps its own native append depth. The
master and PAVE/PURPLE input retain the original caller depths. This prevents
mixing a fragment-support numerator with a caller read-depth denominator.

Consumer projection in merge component0.1.2 also uses the actual native context
TIER/TNC/repeat/MH annotations, retaining interim originals under
`DS_ORIGINAL_INFO_*` and in the unchanged master. Three PASS reference-repeat
counts differ from the earlier reference-only approximation, and two reference
repeats are absent from the native read-context annotation. `SAGE_INFO_FIELDS`
records actual native presence, so absence is preserved rather than carrying
forward an inconsistent approximation. Caller QUAL/FILTER/genotypes stay intact.

## Repeat and realignment safeguards

A reference-flank/ALT substring can match a reference read inside a tandem
repeat. Native context-builder acceptance and a successful roundtrip alone do
not prove that the context encodes the requested allele. The initial
haplotype-only seed experiment is retained as rejected diagnostic evidence.
Only an observed exact CIGAR event is accepted by the seed exporter; the
serialized CIGAR is checked again at merge. No artificial read or CIGAR is made.

Google 1.10 native local-realignment diagnostics are being inspected for the
remaining indels. The diagnostic logger overwrites per-region BAMs across
samples, and its serialization can omit pairing flags. Diagnostic BAMs are
for inspection only, never substitute tumour/normal counting inputs. Native
SAGE_APPEND measurements use the original REDUX BAMs and actual sidecars.

The corrected diagnostic preserves the original 1000-base partition grid.
It recovers the original 91 supporting reads for the chr17 deletion and22
for the chr19 insertion. The native primary contexts' read identities match
the candidate protobuf's tumour sample labels. This supports the two contexts;
it does not change the source calls or imply that SAGE and Google count alike.

`export_deepsomatic_read_context.py` is an optional bounded diagnostic observer
after native allele-counter normalization. It unpacks the image's existing
runtime (including its original symlinks and protobuf path) into scratch,
calls the native entry without changing scientific options, and writes the
already computed reads separately by sample role. It does not install software,
alter reads or replace counting inputs. Job12454496 completed successfully,
exporting14 sample/region BAMs from the7 diagnostic partitions. The underlying
caller reported286 candidates,239 examples and83 small-model examples, matching
the unobserved native diagnostic's totals. This is diagnostic infrastructure,
not a required extra whole-genome example-generation pass.

## Counts must not be forced to agree

Native append on the three observed contexts completed (job12454494). For the
chr17 deletion it measures zero ALT fragments at tumour depth30, while Google's
local assembly supplies91 ALT reads at caller depth141. Both remain visible
with their actual source/counting units. Do not overwrite caller AD with zero,
label this absence of biological support, or copy91 into a fabricated SAGE
fragment-support vector. For the chr19 insertion, native append measures one
tumour ALT fragment at depth17; the Google diagnostic has22 supporting reads.
The first chr17 insertion measures7 tumour ALT fragments at depth89.

These discrepancies remain scientific comparison findings. Successful context
serialization and native consumer execution do not validate count equivalence
or neoantigen sensitivity. The complete downstream comparison must examine
their effect, including NEO, before calling the replacement feature-equivalent.

## Verified scope and remaining work

Native consumer job 12454446 completed successfully: all 180 supplied alleles
reconstructed and survived RNA append; NEO's installed count reader accepted
all 180; native PAVE completed and retained 180 records. The master retained
all 15281 normalized candidate rows, including all 183 PASS alleles, unchanged
in allele, QUAL, FILTER and original genotype payload.

Job12454495 subsequently passed all four native checks for all183 PASS alleles
with no unresolved contexts. Job12454502 repeats that test using the corrected
native-context annotation projection (component0.1.2); its final receipt and
allele-set audit are recorded in the linked execution record below.

This does **not** yet establish:

- downstream equivalence for the recovered indel contexts and count discrepancies;
- whole-genome execution of the evidence-export stage;
- complete phasing/MNV behavior for DeepSomatic-only calls (native informative
  phase-set retention still depends on native passing variants);
- equivalent TINC recovery (the control remains an independent diagnostic;
  its recovered labels must not be copied onto DeepSomatic calls);
- full NEO prediction or PURPLE/LINX/ORANGE results from the substituted arm.

The zero-native-fragment-support audit contains three PASS alleles (one deletion
and two SNVs), not just the deletion discussed above. Their caller AD/DP and
native filters/support remain separate. None is removed or reclassified here.

Keep RNA append, NEO and all downstream branches enabled in the design. Do not
publish this experimental output as a feature-complete `sage_dir` while these
items remain open. Full DeepSomatic inference/postprocessing is independent
and was not changed by this work.

Exact job commands, outputs and updated counts are recorded in
[the execution record](../tests/results/deepsomatic-distributed/EXECUTION.md).

## Pinned scientific sources

SAGE tag `sage-v5.0.2`: `evidence/ReadContextCounters.java`,
`pipeline/RegionTask.java`, `candidate/AltContext.java`,
`candidate/RefContextCache.java`, `candidate/ReadContextCandidate.java`,
`common/VariantReadContextBuilder.java`, `vcf/CandidateSerialisation.java`.
NEO tag `neo-v1.3`: `SomaticVariant.fromContext` and native
`VariantReadSupport`. These are the installed oncoanalyser 3.0 tool versions.
Google r1.10: `make_examples_core.py` (per-sample realignment order),
`make_examples_somatic.py` (normal then tumour) and `realigner/realigner.py`
(native diagnostic output). Source snapshots and hashes are retained alongside
the private execution package; the overlay's two source hashes are versioned.
