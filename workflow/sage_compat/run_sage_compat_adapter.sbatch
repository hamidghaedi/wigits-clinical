#!/usr/bin/env bash
#SBATCH --job-name=sage-compat-adapter
#SBATCH --account=def-hpcg6049_cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
# Usage: sbatch workflow/sage_compat/run_sage_compat_adapter.sbatch <config.json>
# Produces <out_dir>/sage_dir/<tumor_id>.sage.somatic.vcf.gz(.tbi), ALLELE_LEDGER.tsv, ledger_summary.json.
set -eo pipefail
config=$1
test -s "$config"
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
module load apptainer/1.3.5
set -u
get() { python3 -c 'import json,sys; v=json.load(open(sys.argv[1])).get(sys.argv[2], sys.argv[3] if len(sys.argv)>3 else None); print("" if v is None else v)' "$config" "$@"; }

SRC=$(get source_dir)
OUT=$(get out_dir)
VCF=$(get deepsomatic_vcf)
FASTA=$(get reference_fasta)
DICT="${FASTA%.fasta}.dict"; test -s "$DICT" || DICT="${FASTA}.dict"
RES=$(get hmf_resources)
TUMOR=$(get tumor_id); NORMAL=$(get normal_id)
TBAM=$(get tumor_bam); NBAM=$(get normal_bam)
SAGE_IMG=$(get sage_image); BCF_IMG=$(get bcftools_image)
SAGE_VCF=$(get native_sage_somatic_vcf "")
MIN_MAPQ=$(get min_mapq 10); MIN_BASEQ=$(get min_baseq 20)
HOTSPOTS="$RES/dna/variants/KnownHotspots.somatic.38.vcf.gz"
PANEL="$RES/common/DriverGenePanel.38.tsv"
HIGHCONF="$RES/dna/variants/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed.gz"
ENSEMBL="$RES/common/ensembl_data"
for f in "$VCF" "$VCF.tbi" "$FASTA" "$DICT" "$HOTSPOTS" "$PANEL" "$HIGHCONF" "$NBAM" "$SAGE_IMG" "$BCF_IMG" "$SRC/SageCompatAnnotator.java" "$SRC/sage_compat_adapter.py"; do test -e "$f" || { echo "missing: $f" >&2; exit 1; }; done

mkdir -p "$OUT/work" "$OUT/sage_dir" "$OUT/logs"
W="$OUT/work"
log() { echo "$(date -u +%FT%TZ) $*"; }
BCF() { apptainer exec --cleanenv --bind /global:/global "$BCF_IMG" bcftools "$@"; }
JAR=$(apptainer exec --cleanenv "$SAGE_IMG" bash -c 'dirname "$(readlink -f "$(command -v sage)")"')/sage.jar
JAVA() { apptainer exec --cleanenv --bind /global:/global "$SAGE_IMG" java -Xmx12g --class-path "$JAR" "$SRC/SageCompatAnnotator.java" "$@"; }
SAGE_OPTS=(-tumor "$TUMOR" -tumor_bam "$TBAM" -reference "$NORMAL" -reference_bam "$NBAM" -output_vcf "$W/unused.vcf.gz"
           -ref_genome "$FASTA" -ref_genome_version 38 -hotspots "$HOTSPOTS" -high_confidence_bed "$HIGHCONF"
           -driver_gene_panel "$PANEL" -ensembl_data_dir "$ENSEMBL")

log "self-test"; JAVA --self-test | tee "$W/selftest.log"
python3 "$SRC/tests/test_sage_compat_adapter.py" > "$W/unittest.log" 2>&1; tail -1 "$W/unittest.log"

log "filter counts over the full DeepSomatic VCF (background)"
( BCF query -f '%FILTER\n' "$VCF" | sort | uniq -c > "$W/filter_counts.txt" ) &
FC_PID=$!
( sha256sum "$VCF" | cut -d' ' -f1 > "$W/deepsomatic.sha256" ) &
SHA_PID=$!

log "PASS subset + normalisation (split multiallelics, left-align, old-rec-tag)"
# DeepSomatic 1.10 declares FORMAT/NAF as Number=R but writes A values (bcftools norm refuses); NAF is unused here.
BCF view -f PASS -Ou "$VCF" | BCF annotate -x FORMAT/NAF -Ou | BCF norm -f "$FASTA" -m -any -c e --old-rec-tag DS_ORIGINAL_ALLELE -Oz -o "$W/pass.norm.vcf.gz" 2> "$W/norm_pass.log"
cat "$W/norm_pass.log"
BCF index -t "$W/pass.norm.vcf.gz"

log "SAGE panel+hotspot regions (native ReferenceData)"
JAVA panel-bed "$W/panel_raw.bed" "${SAGE_OPTS[@]}" | tee "$W/panel_bed.log"
python3 "$SRC/sage_compat_adapter.py" merge-bed --bed "$W/panel_raw.bed" --dict "$DICT" --out "$W/panel.bed"
log "non-PASS records inside panel regions"
BCF view -R "$W/panel.bed" -e 'FILTER="PASS"' -Ou "$VCF" | BCF annotate -x FORMAT/NAF -Ou | BCF norm -f "$FASTA" -m -any -c e --old-rec-tag DS_ORIGINAL_ALLELE -Oz -o "$W/panel_nonpass.norm.vcf.gz" 2> "$W/norm_panel.log"
cat "$W/norm_panel.log"
BCF query -f '%FILTER\n' "$W/panel_nonpass.norm.vcf.gz" | sort | uniq -c > "$W/panel_nonpass_filters.txt"; cat "$W/panel_nonpass_filters.txt"

log "build candidates"
python3 "$SRC/sage_compat_adapter.py" build --pass-vcf "$W/pass.norm.vcf.gz" --panel-vcf "$W/panel_nonpass.norm.vcf.gz" \
    --dict "$DICT" --tumor-id "$TUMOR" --out-candidates "$W/candidates.tsv" --out-ledger "$W/candidates_ledger.tsv" --out-json "$W/build.json"

wait $SHA_PID
DS_VER=$(BCF view -h "$VCF" | grep -m1 '^##DeepVariant_version=' | cut -d= -f2 || true)
BCF_VER=$(BCF --version | head -1 | awk '{print $2}')
JAVA_VER=$(apptainer exec --cleanenv "$SAGE_IMG" java -version 2>&1 | head -1 | tr -d '"')
{
  echo "##source=deepsomatic-${DS_VER:-unknown}+sage-compat-adapter-0.1.0"
  echo "##sage_compat_adapter_version=0.1.0 (wigits-clinical workflow/sage_compat; annotator uses hmftools sage 5.0.2 classes TierSelector,RepeatInfo,Microhomology)"
  echo "##sage_compat_adapter_date=$(date -u +%FT%TZ)"
  echo "##sage_compat_adapter_source_vcf=$VCF sha256=$(cat "$W/deepsomatic.sha256") deepvariant_version=${DS_VER:-unknown}"
  echo "##sage_compat_adapter_tumor_bam=$TBAM"
  echo "##sage_compat_adapter_normal_bam=$NBAM"
  echo "##sage_compat_adapter_reference=$FASTA dict_sha256=$(sha256sum "$DICT" | cut -d' ' -f1)"
  echo "##sage_compat_adapter_resources=hotspots=$(basename "$HOTSPOTS") sha256=$(sha256sum "$HOTSPOTS" | cut -d' ' -f1);driver_gene_panel=$(basename "$PANEL") sha256=$(sha256sum "$PANEL" | cut -d' ' -f1);high_confidence_bed=$(basename "$HIGHCONF") sha256=$(sha256sum "$HIGHCONF" | cut -d' ' -f1);ensembl_data_dir=$ENSEMBL"
  echo "##sage_compat_adapter_tools=sage_image=$(basename "$SAGE_IMG");bcftools=$BCF_VER;java=$JAVA_VER"
  echo "##sage_compat_adapter_normalisation=bcftools norm -f ref -m -any -c e --old-rec-tag DS_ORIGINAL_ALLELE; alleles retained only when the DeepSomatic tumour GT calls them"
  echo "##sage_compat_adapter_filter_policy=FILTER=PASS for DeepSomatic PASS; non-PASS records retained only when TIER is HOTSPOT or PANEL (DS_FILTER carries the original filter); RefCall/NoCall/other records dropped"
  echo "##sage_compat_adapter_slurm_job=${SLURM_JOB_ID:-none}"
} > "$W/header_extra.txt"

log "annotate (TIER, REP_C/REP_S, TNC, MH, normal pileup) and write sage_dir VCF"
FINAL="$OUT/sage_dir/$TUMOR.sage.somatic.vcf.gz"
rm -f "$FINAL" "$FINAL.tbi"
JAVA annotate "$W/candidates.tsv" "$FINAL" "$W/annotations.tsv" "$NBAM" "$NORMAL" "$TUMOR" "$W/header_extra.txt" "$MIN_MAPQ" "$MIN_BASEQ" "${SAGE_OPTS[@]}" | tee "$W/annotate.log"
test -s "$FINAL.tbi" || BCF index -t "$FINAL"

log "waiting for filter counts"; wait $FC_PID; cat "$W/filter_counts.txt"
log "finalize ledger + validation"
python3 "$SRC/sage_compat_adapter.py" finalize --candidates-ledger "$W/candidates_ledger.tsv" --annotations "$W/annotations.tsv" \
    --filter-counts "$W/filter_counts.txt" --final-vcf "$FINAL" --dict "$DICT" --normal-id "$NORMAL" --tumor-id "$TUMOR" \
    ${SAGE_VCF:+--sage-vcf "$SAGE_VCF"} --out-ledger "$OUT/ALLELE_LEDGER.tsv" --out-json "$OUT/ledger_summary.json"
BCF view -h "$FINAL" | grep -v '^##contig' > "$OUT/sage_dir/header.txt"
cp "$config" "$OUT/adapter_config.json"
log "done: $FINAL"; ls -la "$OUT/sage_dir"
