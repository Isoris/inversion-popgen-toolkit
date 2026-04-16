#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -J inv_D_peel
#SBATCH --array=1-28
#SBATCH -o logs/inv_D_peel.%A_%a.out
#SBATCH -e logs/inv_D_peel.%A_%a.err
#SBATCH --account=lt200308-agbsci

set -euo pipefail

# ============================================================================
# LAUNCH_peel.slurm
#
# Phase 2 / 2d — Phase 9 only (peeling diagnostic) for one chromosome per
# array task. Run after LAUNCH_run_all.slurm has completed all 8 core phases.
# Adds L1 / L1b / L2 peel columns to the existing scoring_table_<chr>.tsv.
#
# Requires per-chromosome pruning table if --pruned-list is used. Build it
# via STEP_D11a_ngsrelate_perchr.sh + STEP_D11b_perchr_pruning_table.R.
# ============================================================================

source ~/.bashrc
mamba activate assembly

CONFIG="${INV_CONFIG:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_codebase_v8.5/00_inversion_config.sh}"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }

set -a
source "${CONFIG}"
set +a

inv_init_dirs

LAUNCHER_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNNER="${LAUNCHER_DIR}/run_all.R"

D_OUTDIR="${D_OUTDIR:-${INVDIR}/06_mds_candidates/snake_regions_multiscale/candidate_detection_out}"
PRUNED_LIST="${PRUNED_LIST:-${D_OUTDIR}/per_chr_pruned.tsv}"

CHR="$(printf "C_gar_LG%02d" "${SLURM_ARRAY_TASK_ID}")"

echo "================================================================"
echo "  STEP_D Phase 9 (peeling diagnostic) — ${CHR}"
echo "  Start: $(date)"
echo "================================================================"

"${RSCRIPT_BIN}" "${RUNNER}" \
  --chr         "${CHR}" \
  --sim-mat-dir "${SIM_MATS_DIR}" \
  --precomp-dir "${PRECOMP_DIR}" \
  --outdir      "${D_OUTDIR}" \
  --pruned-list "${PRUNED_LIST}" \
  --phases      "9"

echo "[DONE] STEP_D peeling — ${CHR} — $(date)"
