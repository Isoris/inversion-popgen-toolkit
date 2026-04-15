#!/usr/bin/env bash
###############################################################################
# STEP_A04_beagles.sh — BEAGLE GL per-RF + whole-genome + merge
# Combines old: S10 (byRF), S11 (wholegenome), S12 (merge)
# Called by: run_step1.sh make_beagles
###############################################################################
set -euo pipefail
source "$(dirname "$0")/../config.sh"

CHUNK_LIST="${GLOBAL_DIR}/chunk_rf.list"
N=$(wc -l < "$CHUNK_LIST")

echo "[$(timestamp)] BEAGLE generation"
echo ""

# ---- Per-RF BEAGLEs (SLURM array per thinning) ----
for W in "${THIN_FINE[@]}"; do
  SITES="${THIN_DIR}/03_sites/sites.thin_${W}.majmin.tsv"
  [[ -s "$SITES" ]] || { echo "[WARN] Missing sites for thin_${W}"; continue; }
  echo "[INFO] thin_${W}: sbatch --array=0-$((N-1))%8 slurm/SLURM_A04a_beagle_rf.sh ${CHUNK_LIST} ${W}"
done

echo ""

# ---- Whole-genome BEAGLEs (one job per broad thinning) ----
for W in "${THIN_BROAD[@]}"; do
  SITES="${THIN_DIR}/03_sites/sites.thin_${W}.majmin.tsv"
  [[ -s "$SITES" ]] || { echo "[WARN] Missing sites for thin_${W}"; continue; }
  echo "[INFO] thin_${W} WG: sbatch slurm/SLURM_A04b_beagle_wg.sh ${W}"
done

echo ""
echo "[INFO] After all BEAGLE jobs complete, merge per-RF to WG:"
echo "  bash steps/STEP_A05_merge_beagles.sh"
