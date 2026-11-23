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
#     installed into $VEP_PLUGINS by tools/install_vep_plugins.sh, not with the VEP image.
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
  --symbol --terms SO --transcript_version --hgvs --flag_pick \
  --pick_order canonical,tsl,biotype,rank,ccds,length \
  --plugin Wildtype --plugin Frameshift --dir_plugins "${VEP_PLUGINS}" \
  --fork "${SLURM_CPUS_PER_TASK:-4}"

grep -q "WildtypeProtein" "${OUTDIR}/${CASE}.vep.vcf" \
  || { echo "[neo_export] ERROR: Wildtype plugin produced no annotation; pVACseq would reject this VCF" >&2; exit 2; }

bgzip -f "${OUTDIR}/${CASE}.vep.vcf" && tabix -f -p vcf "${OUTDIR}/${CASE}.vep.vcf.gz" 2>/dev/null || true
for f in "${OUTDIR}/${CASE}.somatic.selected.vcf" "${OUTDIR}/${CASE}.proximal.vcf"; do
  [ -f "$f" ] && { bgzip -f "$f"; tabix -f -p vcf "$f.gz" 2>/dev/null || true; }
done
ls -lh "${OUTDIR}"
