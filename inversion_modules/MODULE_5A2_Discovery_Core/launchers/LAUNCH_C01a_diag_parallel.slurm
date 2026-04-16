#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -t 4:00:00
#SBATCH -J inv_diag
#SBATCH -o logs/inv_diag_%j.out
#SBATCH -e logs/inv_diag_%j.err
# =============================================================================
# LAUNCH_C01a_diag_parallel.slurm
#
# Runs 5 diagnostic plot scripts in parallel:
#   1. STEP_C01a_diag_profiles.R   — profile plots (01-05, 12, 24-32)
#   2. STEP_C01a_diag_matrices.R   — sim_mat heatmaps (06, 14-17)
#   3. STEP_C01a_diag_mds.R        — MDS scatter plots (07-11, 13, 18-23, 26)
#   4. STEP_C01a_diag_summaries.R  — genome-wide summaries (S01-S06)
#   5. STEP_C01b_plot_heatmaps.R   — split heatmaps (sim_mat + cores PA)
#
# All read the same precomp RDS files independently. No conflicts.
# Total wall time = max of the 5 instead of sum.
#
# Usage:
#   sbatch LAUNCH_C01a_diag_parallel.slurm [--chrom C_gar_LG06]
# =============================================================================
set -euo pipefail
source ~/.bashrc
mamba activate assembly

# ── Config ──
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${SCRIPT_DIR}/../../00_inversion_config.sh"
if [[ -f "${CONFIG}" ]]; then
  set -a; source "${CONFIG}"; set +a
fi

RSCRIPT="${RSCRIPT_BIN:-/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript}"
SNAKES_DIR="${SCRIPT_DIR}/../snakes"
PRECOMP_DIR="${MDS_DIR:-inversion_localpca_v7/06_mds_candidates}/snake_regions_multiscale/precomp"
OUTDIR="${MDS_DIR:-inversion_localpca_v7/06_mds_candidates}/snake_regions_multiscale/plots"
CORES_DIR="${MDS_DIR:-inversion_localpca_v7/06_mds_candidates}/snake_regions_multiscale"

# ── Parse optional chrom filter ──
CHROM_ARG=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --chrom) CHROM_ARG="$2"; shift 2 ;;
    *) shift ;;
  esac
done

mkdir -p "${OUTDIR}" logs

CHROM_FLAG="${CHROM_ARG:-all}"

echo "================================================================"
echo "  Diagnostic Plots (5 Parallel Jobs)"
echo "  Precomp:   ${PRECOMP_DIR}"
echo "  Cores:     ${CORES_DIR}"
echo "  Output:    ${OUTDIR}"
echo "  Chrom:     ${CHROM_FLAG}"
echo "  CPUs:      ${SLURM_CPUS_PER_TASK:-8}"
echo "  JobID:     ${SLURM_JOB_ID:-local}"
echo "  Started:   $(date)"
echo "================================================================"

# ── Launch all 5 in background ──
echo "[LAUNCH] Starting 5 parallel diag jobs..."

${RSCRIPT} "${SNAKES_DIR}/STEP_C01a_diag_profiles.R" \
  "${PRECOMP_DIR}" "${OUTDIR}" "${CHROM_FLAG}" \
  > logs/diag_profiles_${SLURM_JOB_ID:-0}.log 2>&1 &
PID1=$!
echo "  [1/5] Profiles    PID=${PID1}"

${RSCRIPT} "${SNAKES_DIR}/STEP_C01a_diag_matrices.R" \
  "${PRECOMP_DIR}" "${OUTDIR}" "${CHROM_FLAG}" \
  > logs/diag_matrices_${SLURM_JOB_ID:-0}.log 2>&1 &
PID2=$!
echo "  [2/5] Matrices    PID=${PID2}"

${RSCRIPT} "${SNAKES_DIR}/STEP_C01a_diag_mds.R" \
  "${PRECOMP_DIR}" "${OUTDIR}" "${CHROM_FLAG}" \
  > logs/diag_mds_${SLURM_JOB_ID:-0}.log 2>&1 &
PID3=$!
echo "  [3/5] MDS         PID=${PID3}"

${RSCRIPT} "${SNAKES_DIR}/STEP_C01a_diag_summaries.R" \
  "${PRECOMP_DIR}" "${OUTDIR}" "${CHROM_FLAG}" \
  > logs/diag_summaries_${SLURM_JOB_ID:-0}.log 2>&1 &
PID4=$!
echo "  [4/5] Summaries   PID=${PID4}"

${RSCRIPT} "${SNAKES_DIR}/STEP_C01b_plot_heatmaps.R" \
  "${PRECOMP_DIR}" "${CORES_DIR}" "${CORES_DIR}" "${CHROM_FLAG}" \
  > logs/diag_heatmaps_${SLURM_JOB_ID:-0}.log 2>&1 &
PID5=$!
echo "  [5/5] Heatmaps    PID=${PID5}"

echo ""
echo "[LAUNCH] Waiting for all 5 to finish..."

# Wait and report
FAIL=0
for pid in $PID1 $PID2 $PID3 $PID4 $PID5; do
  wait $pid || { echo "[ERROR] PID $pid failed"; FAIL=$((FAIL+1)); }
done

echo ""
echo "================================================================"
if [[ $FAIL -eq 0 ]]; then
  echo "  All 5 diag jobs completed successfully"
else
  echo "  WARNING: $FAIL job(s) failed — check logs/"
fi
echo "  Finished: $(date)"
echo "================================================================"

exit $FAIL
