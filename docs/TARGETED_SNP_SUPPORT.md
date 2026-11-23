# Targeted SNP support

`snp-support/v0.1.0` adds a read-only indexed count probe for an explicit native
AMBER SNPCheck VCF. It uses the existing bcftools 1.21 container; no installation
or caller run is required. Inputs, reference, container and output are arguments.

The DNA pipeline's `modules/tn_genotype_concordance.nf` already uses indexed
bcftools mpileup for common-site evidence. Its genotype and Picard/Broad
haplotype-map workflow is a different scope from the requested native 100-site
probe. This adapter reuses the count approach without importing that pipeline,
running bcftools call or claiming its calibrated tumour-aware verdict.

Run `workflow/probe_snp_support.sbatch <python> <released-script>` followed by
`--targets <native-snpcheck.vcf.gz> --reference <exact-native.fa>`
`--tumor-bam <indexed.bam> --normal-bam <indexed.bam>`
`--tumor-id <SM> --normal-id <SM> --container <bcftools.sif> --out <scratch-dir>`.
The wrapper requests 1 CPU, 4GiB and 20 minutes; it never changes a BAM/report.

Only supplied one-base BED regions are fetched. MQ>=50/BQ>=30, disabled BAQ,
depth cap10000 and explicit excluded alignment flags are recorded in the exact
command. Default bcftools mate-overlap handling remains. These counts are not
claimed identical to AMBER/SAGE fragments or unique molecules. REF/ALT, other
allele counts and DP remain distinct. Requested ALT absence in a measured
pileup is zero support; missing whole-locus observations remain missing.

Outputs are the raw mpileup VCF, exact regions, stderr/version/command metadata,
one row for every requested target, and exploratory coarse AF concordance.
No genotype caller, identity reassignment, mixture fraction or new formal QC
gate is introduced. The exploratory bands and depth condition are documented
in the source, and are not an assay validation. Molecular identity or severe
discordance still requires source/laboratory reconciliation.

Sources: [bcftools mpileup](https://samtools.github.io/bcftools/bcftools.html#mpileup),
[AMBER4.3 SNPCheck and contamination](https://github.com/hartwigmedical/hmftools/blob/amber-v4.3/amber/README.md).
