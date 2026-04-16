#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB
#SBATCH -t 1-00:00:00
#SBATCH -J inv_step08
#SBATCH -o logs/inv_step08.%A_%a.out
#SBATCH -e logs/inv_step08.%A_%a.err

set -euo pipefail

# =============================================================================
# STEP08_beagle_to_dosage_by_chr.slurm
#
# Convert one per-chromosome BEAGLE file to:
#   <chrom>.sites.tsv.gz
#   <chrom>.dosage.tsv.gz
#
# Usage:
#   sbatch --array=0-27%8 STEP08_beagle_to_dosage_by_chr.slurm beagle_chr.list
# =============================================================================

source ~/.bashrc
mamba activate assembly

BEAGLE_LIST="${1:?ERROR: provide beagle_chr.list as argv1}"

IDX="${SLURM_ARRAY_TASK_ID}"
LINE=$((IDX + 1))

BEAGLE_FILE="$(sed -n "${LINE}p" "${BEAGLE_LIST}")"
if [[ -z "${BEAGLE_FILE}" || ! -s "${BEAGLE_FILE}" ]]; then
  echo "[ERROR] Could not read BEAGLE file from ${BEAGLE_LIST} at line ${LINE}" >&2
  exit 1
fi

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
INVDIR="${BASE}/inversion_localpca_v7"
SCRIPTDIR="${BASE}/inversion_codebase_clean_v7_3/clean/MODULE_5A_Inversion_Discovery"
OUTDIR="${INVDIR}/04_dosage_by_chr"

mkdir -p "${OUTDIR}" "${SCRIPTDIR}/logs"

TAG="$(basename "${BEAGLE_FILE}" .beagle.gz)"
CHROM="${TAG#main_qcpass.}"

SITES_OUT="${OUTDIR}/${CHROM}.sites.tsv.gz"
DOSAGE_OUT="${OUTDIR}/${CHROM}.dosage.tsv.gz"

echo "[INFO] ARRAY TASK ID = ${SLURM_ARRAY_TASK_ID}"
echo "[INFO] LINE          = ${LINE}"
echo "[INFO] BEAGLE_FILE   = ${BEAGLE_FILE}"
echo "[INFO] CHROM         = ${CHROM}"

if [[ -s "${SITES_OUT}" && -s "${DOSAGE_OUT}" ]]; then
  echo "[INFO] Outputs already exist for ${CHROM}, skipping"
  exit 0
fi

python3 "${SCRIPTDIR}/core_code/STEP08_beagle_to_dosage_by_chr.py" \
  --beagle "${BEAGLE_FILE}" \
  --outdir "${OUTDIR}"

[[ -s "${SITES_OUT}" ]] || { echo "[ERROR] Missing output: ${SITES_OUT}" >&2; exit 1; }
[[ -s "${DOSAGE_OUT}" ]] || { echo "[ERROR] Missing output: ${DOSAGE_OUT}" >&2; exit 1; }

echo "[DONE] ${CHROM}"
echo "  ${SITES_OUT}"
echo "  ${DOSAGE_OUT}"
