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

# =============================================================================
# LAUNCH_STEP10c_window_sample_belonging.slurm  (v8.5)
#
# Window sample-belonging analysis (5 layers).
#
# Usage:
#   sbatch LAUNCH_STEP10c_window_sample_belonging.slurm [candidate_id]
# =============================================================================

source ~/.bashrc
mamba activate assembly

CONFIG="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_codebase_v8.5/00_inversion_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }

set -a
source "${CONFIG}"
set +a

inv_init_dirs

A2_STEPS="${DISCOVERY2DIR}/steps"

CID="${1:-all}"

"${RSCRIPT_BIN}" "${A2_STEPS}/STEP_D02_window_sample_belonging.R" \
    "${R_CONFIG}" "${CID}"

echo "[DONE] STEP10c — $(date)"
