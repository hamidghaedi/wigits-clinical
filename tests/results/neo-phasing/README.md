# SAGE phase → pVACseq verification

27 November 2026. Adapter source `58fa710`; generic dependent-entry source
`6f65f11`; final origin-ambiguity guard `c3aec02`, on `codex/reusable-native-intake`. Native Oncoanalyser3.0.0,
SAGE5.0.2, VEP115.2, pVACtools5.4.0 and the existing GRCh38 HMF resources.
Original KQ0061/Fir packages and completed native/CH results are untouched.

| Check | Evidence / result |
|---|---|
| Offline primary suite | `python -m unittest discover -s tests`:222 tests,3 existing skips, exit0 |
| Existing interface/export/CH feature tests | `python -m pytest tests/test_neo_interface.py tests/test_neo_export.py workflow/tests/test_chip_curation.py -q`:54 passed |
| SAGE adapter suite | `python -m unittest discover -s workflow/sage_compat/tests`:3 passed |
| Pinned pVAC positive/negative controls | Slurm12516227, exit0,4s. Real parser accepts cis plus actual tumour homozygote, rejects unknown/different phase and different transcript. Real VcfConverter and protein generator produce the expected mutant and WT sequences. |
| Reference export through common launcher | Head12516219,1m11; child12516220,10s, exit0. Only NEO_PREPARE; export status PREPARED, prediction NOT_RUN. |
| Selected input preservation |72/72 record bodies identical, including RNA/DNA/expression. Source headers identify the new adapter. |
| Tumour phase source | Raw native SAGE GT is `./.`; measured positive AD/LPS_RC plus retained PURPLE tumour GT supports10 anchors.56 unresolved and6 native-genotype-only; no genotype inference. |
| Actual converter uptake |28 somatic variant–transcript rows for one two-allele phase group; formerly0. Native source/GT/LPS/support preserved in ledger. |
| Actual sequence change |2426 FASTA entries retained;28 MT sequences change with context, WT sequences unchanged. Every changed real sequence also matches the source proximal amino-acid change at its recorded protein-position offset; synthetic controls independently check expected sequences. |
| Full prediction/review | Head12516221 COMPLETED0:0,26m11; child12516222 COMPLETED0:0,25m12.69 candidate identities retained,24 measured/source-value columns identical;2 best peptides changed with context. |
| Dependent review-only entry | Head12516225 COMPLETED0:0,1m42, afterok:12516221 against its future output. Child12516233 COMPLETED0:0,33s;69/69 rows and all66 columns match the prediction-stage review exactly, without repeating models. |
| Final-source export equivalence | Slurm12516232 exits0 in6s. Final guarded exporter produces selected/proximal VCFs byte-identical to those consumed by the full prediction test. |
| Existing source product preservation |23 original export/CH products retain both SHA256 and mtime; no native/CH task runs. |

Commands use the maintained launcher with `--entry neo-prepare`, then
`--entry neo-predict --existing-export <new export>`, then `--entry neo-review`
with explicit new export/prediction paths. Dependency/time/resume controls retain
their meanings. No native alignment/caller/CH task is in these recovery graphs.

The earlier synthetic probe12516211 failed because its coordinate/flanking
window extended before the synthetic contig start; the test coordinates were
corrected in714e45d. Probe12516214 passed. The first real export exposed native
SAGE's uncalled GT; its source-support handling was corrected in58fa710 and
tested offline and in the pinned tool. Only the newly submitted superseded
reference prediction test12516215/12516216 was cancelled. Historical test evidence
is retained; KQ0061 and other active analyses were never cancelled or modified.

Runtime evidence stays outside Git at
`/global/scratch/hpc6049/neo_phasing_adapter_20261127/`; small receipts are retained
in the coordinator's `data/storage_records/neo_phasing_adapter_20261127/`.
Full outputs and per-allele evidence do not enter Git.

This validates the supported interface, not clinical/vaccine performance or
long-range phasing. Ambiguous/multiple LPS remain unresolved. Germline LPS
numbers cannot establish somatic phase. Nonselected somatic context is recorded
rather than mislabelled as germline by pVACtools. Distinct germline alleles at
selected positions are also recorded as unsupported context rather than being
mislabelled somatic by the position-based native converter. Non-missense context is still
unsupported by pVACtools5.4. No calls are removed from the selected candidate VCF,
no filters/thresholds change, and native NEO remains independent.
