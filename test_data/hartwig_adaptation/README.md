# KQ0205 chr17 + chr19 test FASTQs

**Current source location:** `wigits-clinical`, branch `hartwig_adoptation`.
This directory was migrated intact from the DNA adaptation branch; see the
[migration record](../../docs/migrations/DNA_HARTWIG_ADOPTATION_20261110.md).
The preparation and baseline evidence below describe the historical chromosome
test, not the newer full-genome FASTQ or BAM executions.

Initial test data for branch `hartwig_adoptation`. Extract **all** primary reads
overlapping either chromosome and retrieve their mates, including mates outside
those chromosomes or unmapped. No downsampling, duplicate removal, QC-fail or
mapping-quality filter. Secondary/supplementary alignment records are excluded
so an individual sequencing read is emitted once. Both original source BAMs and
existing analyses remain unchanged.

`inputs.csv` binds FFPE tumour WGS, blood-normal WGS and frozen tumour RNA.
All three BAMs and indexes are local on Frontenac; no transfers are needed.
`prepare_fastqs.py` uses samtools 1.22.1 indexed pair retrieval, collation and
FASTQ conversion. It prefers OQ qualities when present and otherwise uses QUAL.
These BAM-derived reads can retain earlier trimming or hard clipping and are
not a reconstruction of original sequencer FASTQs.

## Native samplesheet

[`samplesheet.KQ0205.chr17_chr19.csv`](samplesheet.KQ0205.chr17_chr19.csv)
is the completed Oncoanalyser FASTQ samplesheet: one analysis group and three
rows for tumour DNA, normal DNA and tumour RNA. Each row contains its R1/R2
paths separated by `;`, as required by native Oncoanalyser 3.0.0. Use this file
as `--input` for the native test; no analysis is launched by creating it.

`library_id:<sample>_chr17_chr19;lane:merged` identifies each combined BAM-derived
test library. These are technical test read-group labels, not claims about the
original physical library or sequencing lane. No flowcell, clinical metadata or
original lane numbers have been invented. Confirmed DNA flowcells are included
from the laboratory manifest; the missing RNA flowcell is left unspecified.
The tumour DNA row supplies the case-level `cancer_type:Carcinoid`, based on the
recorded morphology 8240/3, not a new OncoTree assignment or primary-site choice.
**This is an ORANGE display label, not a validated NEO cohort assignment.**
The staged NEO reference has no `Carcinoid` column. See the consumer-specific
terminology review below before interpreting this field.
`info` accepts documented keys, not arbitrary clinical fields. Specimen
preservation remains in `inputs.csv`; it is not an extra native column.

The remaining recorded case information is available to the adaptation agent at:
```text
/global/scratch/hpc6049/hartwig_test_fastqs_20261109/KQ0205/case_metadata.json
```
This private reference includes the bound clinical and laboratory source rows,
histology/site codes, grade, pathology cellularity, dates, specimen preparation
and library kits, with source paths. Source discrepancies are preserved. Native
Oncoanalyser does not automatically consume this JSON. Clinical payload stays
outside Git; no new metadata-ingestion code is introduced here.

Supported `info` keys are `library_id`, `lane`, `flowcell`, `cancer_type`,
`read_group_overrides`, `longitudinal_sample`, and `generate_redux_tsvs_only`.
The last three are optional operational settings, not additional clinical fields.
References: [native input documentation](https://nf-co.re/oncoanalyser/3.0.0/docs/usage/#info-field-keys)
and [morphology terminology](https://seer.cancer.gov/manuals/2026/SPCSM_2026_Appendix_E.pdf).

Consumer meanings, exact reference vocabulary, source limitations and test
evidence are recorded in [METADATA_CONSUMERS.md](METADATA_CONSUMERS.md).
Schema/parser acceptance alone does not establish downstream semantic validity.

Format authority: native `assets/schema_input.json`, `docs/usage.md` (FASTQ
inputs), and `lib/Utils.groovy` at Oncoanalyser 3.0.0 commit
`7c74c87a43749952b38c9a18915947570f0595a0`. Samplesheet metadata can be complete
before extraction finishes: the datasets are ready only when all three
`completion.json` files exist and their Slurm jobs have completed successfully.
Validation passed against the official JSON schema and the unmodified native
`Utils.parseInput` parser: one group, tumour DNA, normal DNA and tumour RNA.
This parser-only check launched no analysis processes.

Outputs are external to Git:
```text
/global/scratch/hpc6049/hartwig_test_fastqs_20261109/KQ0205/
  KQ0205-DNA/KQ0205-DNA.chr17_chr19_R{1,2}.fastq.gz
  KQ0205-gDNA/KQ0205-gDNA.chr17_chr19_R{1,2}.fastq.gz
  KQ0205-RNA/KQ0205-RNA.chr17_chr19_R{1,2}.fastq.gz
```

Each sample also has `commands.log`, `completion.json`, and separate singleton
and other-read FASTQs. Any such reads are retained, not silently discarded.
Completion requires matching paired names/counts and accounting for every
selected primary BAM record. The generated subset BAM is removed only after
that check succeeds. A pre-existing output directory is never overwritten.

The owner authorizes this standalone Slurm preparation for the adaptation test;
it does not change production Nextflow entry points or scientific results.

## Exact preparation launch

Run from this directory after committing the scripts. Logs go on scratch.
```bash
runtime=/global/scratch/hpc6049/hartwig_test_fastqs_20261109/KQ0205
mkdir -p "$runtime/logs" "$runtime/package"
# Stage committed bytes so later branch edits do not change a queued job.
for filename in prepare_fastqs.py prepare_fastqs.sbatch inputs.csv; do
    git show "HEAD:test_data/hartwig_adaptation/$filename" > "$runtime/package/$filename"
done
git rev-parse HEAD > "$runtime/package/source_commit.txt"
package_dir="$runtime/package"
python3 - "$package_dir" "$runtime" <<'PY'
import csv, subprocess, sys
from pathlib import Path
package, runtime = map(Path, sys.argv[1:])
for row in csv.DictReader((package / 'inputs.csv').open()):
    sample = row['sample_id']
    # BAI mapped-record counts across chr17/19: 208.6M / 91.1M / 88.4M.
    # Conservative, uncalibrated time allowance: 2h + 2h per 50M records,
    # rounded up. Mates can add records; limits are 12h/6h/6h with headroom.
    walltime = '12:00:00' if sample == 'KQ0205-DNA' else '06:00:00'
    command = ['sbatch', '--parsable', '--account=def-hpcg6049_cpu',
               '--time=' + walltime, '--output=' + str(runtime / 'logs' / (sample + '.%j.log')),
               str(package / 'prepare_fastqs.sbatch'), str(package / 'prepare_fastqs.py'),
               row['bam'], sample, str(runtime / sample)]
    print(sample, subprocess.check_output(command, text=True).strip(), flush=True)
PY
```

Allocation is eight CPUs and 24 GiB RAM per library. Indexed BAI counts include
non-primary alignments and are an input-size proxy, not expected FASTQ counts.
Time allowances are conservative estimates, not measured performance claims.
The request fits the site's general CPU nodes; no partition is hardcoded.

`test_prepare_fastqs.py` is a small synthetic end-to-end check, including mates
on other chromosomes, an unmapped mate, low-MAPQ duplicate reads, OQ recovery,
secondary exclusion and separate singleton retention. Run it with the same
samtools module loaded before preparing the real data.
