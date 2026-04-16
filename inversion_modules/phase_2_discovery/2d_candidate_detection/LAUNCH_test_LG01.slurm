#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH -t 00:30:00
#SBATCH -J inv_D_test
#SBATCH -o logs/inv_D_test.%j.out
#SBATCH -e logs/inv_D_test.%j.err
#SBATCH --account=lt200308-agbsci

set -euo pipefail

# ============================================================================
# LAUNCH_test_LG01.slurm
#
# Quick smoke test on LG01 only: staircase + block scoring (phases 1..3),
# no matrix transforms, no NN sweep, no peeling. ~5–10 min.
# Use this first to verify the module works end-to-end before running the
# full 28-chr array.
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
mkdir -p "${D_OUTDIR}" logs

CHR="${1:-C_gar_LG01}"

echo "================================================================"
echo "  STEP_D run_all — smoke test on ${CHR}"
echo "  Phases 1..3 only (staircase + block scoring, no transforms)"
echo "  Start: $(date)"
echo "================================================================"

"${RSCRIPT_BIN}" "${RUNNER}" \
  --chr              "${CHR}" \
  --sim-mat-dir      "${SIM_MATS_DIR}" \
  --precomp-dir      "${PRECOMP_DIR}" \
  --outdir           "${D_OUTDIR}" \
  --phases           "1:3" \
  --skip-transforms  TRUE

echo "[DONE] STEP_D smoke test — ${CHR} — $(date)"
