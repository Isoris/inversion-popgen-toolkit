#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH -t 06:00:00
#SBATCH -J inv_D_all
#SBATCH --array=1-28
#SBATCH -o logs/inv_D_all.%A_%a.out
#SBATCH -e logs/inv_D_all.%A_%a.err
#SBATCH --account=lt200308-agbsci

set -euo pipefail

# ============================================================================
# LAUNCH_run_all.slurm
#
# Phase 2 / 2d — candidate detection (matrix-based track). One chromosome
# per SLURM array task. Runs all 9 phases of run_all.R.
#
# Memory: 48G for ~9000×9000 matrices + 6 matrix-transform variants.
# Time:   ~4–6 h worst case. Use --skip-transforms for initial testing.
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
[[ -s "${RUNNER}" ]] || { echo "[ERROR] Missing runner: ${RUNNER}" >&2; exit 1; }

D_OUTDIR="${D_OUTDIR:-${INVDIR}/06_mds_candidates/snake_regions_multiscale/candidate_detection_out}"
mkdir -p "${D_OUTDIR}" "${INVDIR}/logs" logs

CHR="$(printf "C_gar_LG%02d" "${SLURM_ARRAY_TASK_ID}")"

echo "================================================================"
echo "  STEP_D run_all — ${CHR}"
echo "  Array task: ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_TASK_MAX:-?}"
echo "  SimMats:    ${SIM_MATS_DIR}"
echo "  Out:        ${D_OUTDIR}"
echo "  Start:      $(date)"
echo "================================================================"

"${RSCRIPT_BIN}" "${RUNNER}" \
  --chr         "${CHR}" \
  --sim-mat-dir "${SIM_MATS_DIR}" \
  --precomp-dir "${PRECOMP_DIR}" \
  --sv-file     "${DELLY_INV_VCF:-}" \
  --outdir      "${D_OUTDIR}" \
  --nn-scales   "0,20,40,80" \
  --phases      "all"

echo "[DONE] STEP_D run_all — ${CHR} — $(date)"
