#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB
#SBATCH -t 1-00:00:00
#SBATCH -J inv_step09
#SBATCH -o logs/inv_step09.%A_%a.out
#SBATCH -e logs/inv_step09.%A_%a.err

set -euo pipefail

source ~/.bashrc
mamba activate assembly

RSCRIPT_BIN="/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript"

CHROM_LIST="${1:?ERROR: provide chrom list as argv1}"

IDX="${SLURM_ARRAY_TASK_ID}"
LINE=$((IDX + 1))

CHROM="$(sed -n "${LINE}p" "${CHROM_LIST}")"
if [[ -z "${CHROM}" ]]; then
  echo "[ERROR] Could not read chromosome from ${CHROM_LIST} line ${LINE}" >&2
  exit 1
fi

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
INVDIR="${BASE}/inversion_localpca_v7"
SCRIPTDIR="${BASE}/inversion_codebase_clean_v7_3/clean/MODULE_5A_Inversion_Discovery"

mkdir -p "${INVDIR}/05_local_pca" "${SCRIPTDIR}/logs"

SITES="${INVDIR}/04_dosage_by_chr/${CHROM}.sites.tsv.gz"
DOSAGE="${INVDIR}/04_dosage_by_chr/${CHROM}.dosage.tsv.gz"
OUTPREFIX="${INVDIR}/05_local_pca/STEP09_${CHROM}"

[[ -x "${RSCRIPT_BIN}" ]] || { echo "[ERROR] Rscript not executable: ${RSCRIPT_BIN}" >&2; exit 1; }
[[ -s "${SITES}" ]] || { echo "[ERROR] Missing ${SITES}" >&2; exit 1; }
[[ -s "${DOSAGE}" ]] || { echo "[ERROR] Missing ${DOSAGE}" >&2; exit 1; }

if [[ -s "${OUTPREFIX}.window_meta.tsv.gz" && -s "${OUTPREFIX}.window_pca.rds" && -s "${OUTPREFIX}.window_pca.tsv.gz" ]]; then
  echo "[INFO] STEP09 outputs already exist for ${CHROM}, skipping"
  exit 0
fi

echo "[INFO] CHROM=${CHROM}"
echo "[INFO] Using Rscript: ${RSCRIPT_BIN}"
"${RSCRIPT_BIN}" -e 'cat("R=", R.version.string, "\n"); cat("Rscript=", Sys.which("Rscript"), "\n"); suppressPackageStartupMessages(library(RSpectra)); cat("RSpectra OK\n")'

"${RSCRIPT_BIN}" "${SCRIPTDIR}/core_code/STEP09_local_pca_windows_by_chr.R" \
  "${SITES}" \
  "${DOSAGE}" \
  "${OUTPREFIX}" \
  2 \
  100

echo "[DONE] ${CHROM}"
