#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128GB
#SBATCH -t 1-00:00:00
#SBATCH -J inv_step10
#SBATCH -o logs/inv_step10.%j.out
#SBATCH -e logs/inv_step10.%j.err

set -euo pipefail

source ~/.bashrc
mamba activate assembly

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
INVDIR="${BASE}/inversion_localpca_v7"
SCRIPTDIR="${BASE}/inversion_codebase_clean_v7_3/clean/MODULE_5A_Inversion_Discovery"

mkdir -p "${INVDIR}/06_mds_candidates" "${SCRIPTDIR}/logs"

"${RSCRIPT_BIN}" "${SCRIPTDIR}/core_code/STEP10_window_mds_outliers.R" \
  "${INVDIR}/05_local_pca" \
  "${INVDIR}/06_mds_candidates/inversion_localpca" \
  2 20 3 500000 3
