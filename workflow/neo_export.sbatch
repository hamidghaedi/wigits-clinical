#!/usr/bin/env bash
#SBATCH --account=def-hpcg6049_cpu
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --job-name=neo_export
#
# Export one oncoanalyser case to pVACseq inputs, then annotate with VEP.
#
#   sbatch workflow/neo_export.sbatch <results/<case> dir> [outdir]
#
# Defaults to the KQ0004 development case. Paths come from workflow/conf/frontenac_paths.env
# (docs/directory_layout.md). VEP notes recorded in docs/neoantigen_interface.md:
#   * the merged offline cache is BAM-edited, so --use_given_ref is mandatory, otherwise VEP
#     replaces the VCF's REF allele with the transcript reference;
#   * pVACseq requires the Wildtype and Frameshift plugins, which ship with pVACtools and are
#     installed into $VEP_PLUGINS by `pvacseq install_vep_plugin`, not with the VEP image;
#   * pVACseq also requires --tsl and --biotype, and it rejects the run unless the PROXIMAL
#     variants VCF is itself VEP-annotated, bgzipped, tabixed and carries exactly one sample.
set -eo pipefail
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
module load apptainer/1.3.5
set -u

REPO="${REPO:-/global/project/hpcg6049/somatic_pipeline/wigits-clinical}"
source "${REPO}/workflow/conf/frontenac_paths.env"

CASE_DIR="${1:-/global/scratch/hpc6049/wigits_runs/KQ0004_full_native_20261114_onco3.0.0/results/KQ0004}"
CASE="$(basename "${CASE_DIR}")"
OUTDIR="${2:-${RUN_ROOT}/neo_ext/${CASE}}"
VEP_SIF="${CONTAINER_DIR}/ensembl-vep_release_115.2.sif"
VEP_CACHE="/global/scratch/hpc6049/cancer_dbs/vep/cache"
VEP_PLUGINS="${KB_ROOT}/vep_plugins/pvactools"
FASTA="${ONCO_REF_ROOT}/3.0.0/genomes/GRCh38_hmf/GRCh38_masked_exclusions_alts_hlas.fasta"

mkdir -p "${OUTDIR}"
export TMPDIR="${SLURM_TMPDIR:-$PWD}"

echo "[neo_export] case=${CASE} outdir=${OUTDIR}"
python3 "${REPO}/workflow/neo_export.py" --results-case-dir "${CASE_DIR}" --outdir "${OUTDIR}" --case "${CASE}"

SEL="${OUTDIR}/${CASE}.somatic.selected.vcf"
echo "[neo_export] VEP annotation"
apptainer exec -B /global -B "${TMPDIR}" "${VEP_SIF}" vep \
  --input_file "${SEL}" --output_file "${OUTDIR}/${CASE}.vep.vcf" --format vcf --vcf --force_overwrite \
  --offline --cache --dir_cache "${VEP_CACHE}" --cache_version 115 --merged --assembly GRCh38 \
  --fasta "${FASTA}" --use_given_ref \
  --symbol --terms SO --tsl --biotype --transcript_version --hgvs --flag_pick \
  --pick_order canonical,tsl,biotype,rank,ccds,length \
  --plugin Wildtype --plugin Frameshift --dir_plugins "${VEP_PLUGINS}" \
  --fork "${SLURM_CPUS_PER_TASK:-4}"

grep -q "WildtypeProtein" "${OUTDIR}/${CASE}.vep.vcf" \
  || { echo "[neo_export] ERROR: Wildtype plugin produced no annotation; pVACseq would reject this VCF" >&2; exit 2; }

# The proximal VCF needs the same annotation; pVACseq checks its CSQ header, not only the input's.
apptainer exec -B /global -B "${TMPDIR}" "${VEP_SIF}" vep \
  --input_file "${OUTDIR}/${CASE}.proximal.vcf" --output_file "${OUTDIR}/${CASE}.proximal.vep.vcf" \
  --format vcf --vcf --force_overwrite --offline --cache --dir_cache "${VEP_CACHE}" --cache_version 115 \
  --merged --assembly GRCh38 --fasta "${FASTA}" --use_given_ref \
  --symbol --terms SO --tsl --biotype --transcript_version --hgvs \
  --plugin Wildtype --plugin Frameshift --dir_plugins "${VEP_PLUGINS}" --fork "${SLURM_CPUS_PER_TASK:-4}"

# bgzip/tabix come from the samtools image: htslib is not on the pVACtools or VEP image PATH.
apptainer exec -B /global "${CONTAINER_DIR}/samtools-1.22.1--h96c455f_0.sif" sh -c \
  "cd ${OUTDIR} && for v in ${CASE}.vep ${CASE}.proximal.vep; do bgzip -kf \$v.vcf; tabix -fp vcf \$v.vcf.gz; done"
ls -lh "${OUTDIR}"
