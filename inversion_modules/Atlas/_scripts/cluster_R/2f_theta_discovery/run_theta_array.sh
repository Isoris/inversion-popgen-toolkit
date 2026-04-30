#!/bin/bash
# =============================================================================
# run_theta_array.sh — SLURM array launcher for the θπ phase-2 pipeline
# =============================================================================
# 28-chrom array. Each task runs:
#   STEP_TR_A   (heavy precomp from ANGSD pestPG)
#   STEP_TR_B   (4-layer JSON emit)
#   STEP_TR_C   (D17 wrapper, NEW turn 4)
#   STEP_TR_D   (in-place augment with theta_d17_envelopes layer)
#
# Prerequisites:
#   - ANGSD pestPG files at PESTPG_DIR/PESTPG_SCALE (per 00_theta_config.sh)
#   - D17 scripts at the paths below — adjust if needed
#
# Usage:
#   sbatch run_theta_array.sh
# =============================================================================

#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-28
#SBATCH --output=logs/theta_%A_%a.out
#SBATCH --error=logs/theta_%A_%a.err

set -euo pipefail

# ---- Environment ------------------------------------------------------------
source ~/.bashrc
mamba activate assembly

cd /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_codebase_v8.5/phase_2_discovery/2f_theta_discovery

source 00_theta_config.sh   # exports CHROM_LIST, RSCRIPT, JSON_OUT_DIR, OUTROOT, etc.

CHROM=${CHROM_LIST[$((SLURM_ARRAY_TASK_ID - 1))]}

# D17 scripts — EDIT TO MATCH YOUR LAYOUT
D17_L1_SCRIPT="/path/to/STEP_D17_multipass_L1_only_v7.R"
D17_L2_SCRIPT="/path/to/STEP_D17_multipass_L2_v8.R"

THETA_D17_DIR="$OUTROOT/05_d17"
mkdir -p "$THETA_D17_DIR" logs

# ---- Pipeline ---------------------------------------------------------------
echo "[$(date)] [$CHROM] === STEP_TR_A precompute ==="
$RSCRIPT STEP_TR_A_compute_theta_matrices.R --chrom "$CHROM"

echo "[$(date)] [$CHROM] === STEP_TR_B classifier ==="
$RSCRIPT STEP_TR_B_classify_theta.R --chrom "$CHROM"

echo "[$(date)] [$CHROM] === STEP_TR_C D17 wrapper ==="
$RSCRIPT STEP_TR_C_theta_d17_wrapper.R \
  --chrom          "$CHROM" \
  --theta_json     "$JSON_OUT_DIR/$CHROM/${CHROM}_phase2_theta.json" \
  --out_dir        "$THETA_D17_DIR" \
  --d17_l1_script  "$D17_L1_SCRIPT" \
  --d17_l2_script  "$D17_L2_SCRIPT"

# To speed up D17 if needed, uncomment:
# --d17_l1_args "--boundary_grow_W_pct 0.005,0.02,0.05"

echo "[$(date)] [$CHROM] === STEP_TR_D augment JSON ==="
$RSCRIPT STEP_TR_D_augment_theta_json.R \
  --theta_json   "$JSON_OUT_DIR/$CHROM/${CHROM}_phase2_theta.json" \
  --d17_l1_env   "$THETA_D17_DIR/${CHROM}_theta_d17L1_envelopes.tsv" \
  --d17_l2_env   "$THETA_D17_DIR/${CHROM}_theta_d17L2_envelopes.tsv" \
  --d17_l1_bnd   "$THETA_D17_DIR/${CHROM}_theta_d17L1_boundaries.tsv" \
  --d17_l2_bnd   "$THETA_D17_DIR/${CHROM}_theta_d17L2_boundaries.tsv"

echo "[$(date)] [$CHROM] === DONE ==="
ls -lh "$JSON_OUT_DIR/$CHROM/${CHROM}_phase2_theta.json"
