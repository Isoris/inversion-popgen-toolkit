#!/bin/bash
# =============================================================================
# run_theta_array_NO_D17.sh — slim launcher, STEP_TR_A + STEP_TR_B only
# =============================================================================
# Companion to run_theta_array.sh, but skips STEP_TR_C and STEP_TR_D
# (the D17 boundary-refinement passes). Produces the 4-layer phase-2 θπ
# JSON which is what page 12 of the atlas actually renders from. The
# 5th layer (theta_d17_envelopes) can be augmented later with STEP_TR_C+D
# once D17 script paths are known.
#
# Use this for the FIRST production run. Runs in ~10–20 minutes wall time
# across all 28 chromosomes (parallel SLURM array, queue permitting).
#
# Usage:
#   sbatch run_theta_array_NO_D17.sh
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
mkdir -p logs

CHROM=${CHROM_LIST[$((SLURM_ARRAY_TASK_ID - 1))]}

# ---- Pipeline (NO D17) ------------------------------------------------------
echo "[$(date)] [$CHROM] === STEP_TR_A precompute ==="
"$RSCRIPT" STEP_TR_A_compute_theta_matrices.R --chrom "$CHROM"

echo "[$(date)] [$CHROM] === STEP_TR_B classifier ==="
"$RSCRIPT" STEP_TR_B_classify_theta.R --chrom "$CHROM"

echo "[$(date)] [$CHROM] === DONE (no D17) ==="
ls -lh "$JSON_OUT_DIR/$CHROM/${CHROM}_phase2_theta.json"
