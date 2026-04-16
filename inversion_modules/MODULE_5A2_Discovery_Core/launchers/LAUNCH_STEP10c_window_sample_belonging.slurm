#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH -t 12:00:00
#SBATCH -J inv_10c_belonging
#SBATCH -o logs/inv_10c_belonging.%A_%a.out
#SBATCH -e logs/inv_10c_belonging.%A_%a.err

set -euo pipefail

source ~/.bashrc
mamba activate assembly

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
RSCRIPT_BIN="/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript"
DISCOVERYDIR="${BASE}/inversion_codebase_clean_v7_3/clean/MODULE_5A_Inversion_Discovery"
FOLLOWUPDIR="${BASE}/inversion_codebase_clean_v7_3/clean/MODULE_5B_Inversion_Followup"
CONFIG="${FOLLOWUPDIR}/config_inversion_followup.R"

CID="${1:-all}"

"${RSCRIPT_BIN}" "${DISCOVERYDIR}/core_code/STEP10c_window_sample_belonging.R" \
    "${CONFIG}" "${CID}"
