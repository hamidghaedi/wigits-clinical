#!/bin/bash
# Call, annotate and curate clonal haematopoiesis variants for one case.
#
#   chip_call.sh <case> <blood_sample> <tumour_sample> <oncoanalyser_results_dir> <out_dir> [mean_blood_depth]
#
# Runs downstream of nf-core/oncoanalyser and reads its REDUX alignments. SAGE is run over the blood
# alignment with the blood in the sample-of-interest slot and the tumour as reference, as oncoanalyser
# does for germline calling, and without the -germline flag, which penalises allele fractions away from
# inherited values. Every SAGE record is kept: the FILTER verdict is moved into INFO so that PAVE
# annotates all of them. Curation then follows Vlasschaert et al., Blood 2023;141(18):2214-2223.
#
# Paths default to the Frontenac layout and can be overridden from the environment.
set -euo pipefail

CASE=$1; BLOOD=$2; TUMOUR=$3; RESULTS=$4; OUT=$5; MEAN_DEPTH=${6:-}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PANELS=$HERE/panels

CACHE=${CACHE:-/global/project/hpcg6049/somatic_pipeline/containers/NXF_SINGULARITY_CACHEDIR}
REFROOT=${REFROOT:-/global/project/hpcg6049/somatic_pipeline/data/references/oncoanalyser/3.0.0}
RES=${RES:-$REFROOT/hmf_pipeline_resources.38_v3.0.0--8}
REF=${REF:-$REFROOT/genomes/GRCh38_hmf/GRCh38_masked_exclusions_alts_hlas.fasta}
SAGE_SIF=${SAGE_SIF:-$CACHE/quay.io-biocontainers-hmftools-sage-5.0.2--hdfd78af_0.img}
PAVE_SIF=${PAVE_SIF:-$(ls "$CACHE"/*pave* | head -1)}
SAMTOOLS_SIF=${SAMTOOLS_SIF:-$CACHE/samtools-1.22.1--h96c455f_0.sif}
THREADS=${THREADS:-8}
BIND=${BIND:-/global}
TMP=${SLURM_TMPDIR:-${TMPDIR:-/tmp}}
export APPTAINER_TMPDIR=$TMP TMPDIR=$TMP

HIGH_CONF=$RES/dna/variants/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed.gz
ENSEMBL=$RES/common/ensembl_data
STAGE=$OUT/stage
mkdir -p "$OUT" "$STAGE" "$OUT/pave"

echo "=== resources"
python3 "$HERE/build_chip_panel.py" "$PANELS/chip_genes.tsv" "$RES/common/DriverGenePanel.38.tsv" "$OUT/ChipGenePanel.38.tsv"
python3 "$HERE/build_chip_hotspots.py" "$RES/dna/variants/KnownHotspots.somatic.38.vcf.gz" "$PANELS" "$OUT/KnownHotspots.chip.38.vcf"
apptainer exec -B "$BIND" "$SAMTOOLS_SIF" sh -c "cd $OUT && bgzip -f KnownHotspots.chip.38.vcf && tabix -fp vcf KnownHotspots.chip.38.vcf.gz"

echo "=== stage both alignments into one directory"
# SAGE looks for each sample's base-quality recalibration table beside the reference alignment.
for s in "$BLOOD" "$TUMOUR"; do
  for f in "$RESULTS"/alignments/"$s"/"$s".redux.*; do ln -sf "$f" "$STAGE/$(basename "$f")"; done
done

echo "=== sage"
apptainer exec -B "$BIND" -B "$TMP" "$SAGE_SIF" sage -Xmx40G \
  -tumor "$BLOOD" -tumor_bam "$STAGE/$BLOOD.redux.bam" \
  -reference "$TUMOUR" -reference_bam "$STAGE/$TUMOUR.redux.bam" \
  -ref_sample_count 0 -ref_genome "$REF" -ref_genome_version 38 \
  -hotspots "$OUT/KnownHotspots.chip.38.vcf.gz" -driver_gene_panel "$OUT/ChipGenePanel.38.tsv" \
  -high_confidence_bed "$HIGH_CONF" -ensembl_data_dir "$ENSEMBL" -sequencing_type ILLUMINA \
  -panel_only -threads "$THREADS" -output_vcf "$OUT/$CASE.chip.raw.vcf.gz" 2>&1 | tail -4
[ -s "$OUT/$CASE.chip.raw.vcf.gz" ] || { echo "SAGE wrote nothing"; exit 1; }

echo "=== keep every record: move the SAGE verdict into INFO"
python3 - "$OUT/$CASE.chip.raw.vcf.gz" "$OUT/$CASE.chip.allpass.vcf" <<'PY'
import gzip, sys
n = 0
with gzip.open(sys.argv[1], "rt") as fh, open(sys.argv[2], "w") as out:
    for line in fh:
        if line.startswith("##"):
            out.write(line)
        elif line.startswith("#CHROM"):
            out.write('##INFO=<ID=SAGE_FILTER,Number=1,Type=String,'
                      'Description="FILTER value assigned by SAGE before it was moved to INFO">\n')
            out.write(line)
        else:
            f = line.rstrip("\n").split("\t")
            f[7] += ";SAGE_FILTER=" + f[6].replace(";", ",")
            f[6] = "PASS"
            out.write("\t".join(f) + "\n")
            n += 1
print(f"records {n}")
PY
apptainer exec -B "$BIND" "$SAMTOOLS_SIF" sh -c "cd $OUT && bgzip -f $CASE.chip.allpass.vcf && tabix -fp vcf $CASE.chip.allpass.vcf.gz"

echo "=== pave"
apptainer exec -B "$BIND" -B "$TMP" "$PAVE_SIF" pave -Xmx40G \
  -sample "$BLOOD" -input_vcf "$OUT/$CASE.chip.allpass.vcf.gz" \
  -ref_genome "$REF" -ref_genome_version 38 \
  -gnomad_freq_dir "$RES/dna/variants/gnomad" -gnomad_no_filter \
  -clinvar_vcf "$RES/dna/variants/clinvar.38.vcf.gz" -driver_gene_panel "$OUT/ChipGenePanel.38.tsv" \
  -mappability_bed "$RES/dna/variants/mappability_150.38.bed.gz" -ensembl_data_dir "$ENSEMBL" \
  -sequencing_type ILLUMINA -threads "$THREADS" -output_dir "$OUT/pave" \
  -output_vcf "$OUT/pave/$CASE.chip.annotated.vcf.gz" 2>&1 | tail -2
[ -s "$OUT/pave/$CASE.chip.annotated.vcf.gz" ] || { echo "PAVE wrote nothing"; exit 1; }

echo "=== curate"
python3 "$HERE/chip_candidates.py" --annotated-vcf "$OUT/pave/$CASE.chip.annotated.vcf.gz" \
  --outdir "$OUT" --case "$CASE" --blood-sample "$BLOOD" --tumour-sample "$TUMOUR" \
  --panel-file "$PANELS/chip_genes.tsv" --panel-dir "$PANELS" ${MEAN_DEPTH:+--mean-depth "$MEAN_DEPTH"}
python3 "$HERE/chip_review_table.py" --annotated-vcf "$OUT/pave/$CASE.chip.annotated.vcf.gz" \
  --variants-tsv "$OUT/$CASE.chip_variants.tsv" --outdir "$OUT" --case "$CASE" \
  --blood-sample "$BLOOD" --tumour-sample "$TUMOUR" --panel-dir "$PANELS"
echo "CHIP_CALL_DONE $CASE"
