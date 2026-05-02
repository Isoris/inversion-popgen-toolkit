#!/bin/bash
# =============================================================================
# run_ghsl_array.sh — SLURM array launcher for the GHSL phase-2 pipeline
# =============================================================================
# 28-chrom array. Each task runs:
#   STEP_C04c (local PCA)
#   STEP_C04d (D17 wrapper)
#   export_ghsl_to_json_v3 (consolidated emit)
#
# Prerequisites:
#   - STEP_C04 v6 already run, RDS files in $GHSL_V6_DIR
#   - STEP_C04b classifier already run (the annot RDS carries PASS/WEAK/FAIL)
#   - D17 scripts at the paths below — adjust if needed
#
# Usage:
#   sbatch run_ghsl_array.sh
#
# To run a single chrom interactively for debugging:
#   SLURM_ARRAY_TASK_ID=28 bash run_ghsl_array.sh    # LG28
# =============================================================================

#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --array=1-28
#SBATCH --output=logs/ghsl_%A_%a.out
#SBATCH --error=logs/ghsl_%A_%a.err

set -euo pipefail

# ---- Environment ------------------------------------------------------------
source ~/.bashrc
mamba activate assembly

RSCRIPT=/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript
BASE=/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_codebase_v8.5/phase_2_discovery
cd "$BASE"

# ---- Chromosome resolution from array task id -------------------------------
CHROM_LIST=(C_gar_LG01 C_gar_LG02 C_gar_LG03 C_gar_LG04 C_gar_LG05
            C_gar_LG06 C_gar_LG07 C_gar_LG08 C_gar_LG09 C_gar_LG10
            C_gar_LG11 C_gar_LG12 C_gar_LG13 C_gar_LG14 C_gar_LG15
            C_gar_LG16 C_gar_LG17 C_gar_LG18 C_gar_LG19 C_gar_LG20
            C_gar_LG21 C_gar_LG22 C_gar_LG23 C_gar_LG24 C_gar_LG25
            C_gar_LG26 C_gar_LG27 C_gar_LG28)
CHROM=${CHROM_LIST[$((SLURM_ARRAY_TASK_ID - 1))]}

# ---- Paths — EDIT THESE TO MATCH YOUR LAYOUT --------------------------------
GHSL_V6_DIR="2e_ghsl_discovery/02_ghsl_v6_outputs"   # <-- where STEP_C04 v6 wrote RDSes
GHSL_LOCALPCA_DIR="2e_ghsl_discovery/03_localpca"
GHSL_D17_DIR="2e_ghsl_discovery/04_d17"
JSON_OUT="2e_ghsl_discovery/05_json_out"

# D17 scripts — EDIT TO MATCH YOUR LAYOUT
D17_L1_SCRIPT="/path/to/STEP_D17_multipass_L1_only_v7.R"
D17_L2_SCRIPT="/path/to/STEP_D17_multipass_L2_v8.R"

mkdir -p "$GHSL_LOCALPCA_DIR" "$GHSL_D17_DIR" "$JSON_OUT" logs

# ---- Pipeline ---------------------------------------------------------------
echo "[$(date)] [$CHROM] === STEP_C04c local PCA ==="
$RSCRIPT 2e_ghsl_discovery/STEP_C04c_ghsl_local_pca.R \
  --chrom         "$CHROM" \
  --ghsl_matrices "$GHSL_V6_DIR/${CHROM}.ghsl_v6_matrices.rds" \
  --out_rds       "$GHSL_LOCALPCA_DIR/${CHROM}.ghsl_v6_localpca.rds" \
  --pad           1 \
  --smoothing-scale none \
  --anchor        auto \
  --min-samples-per-window 50 \
  --z-threshold           2.5 \
  --min-l2-windows        5 \
  --merge-gap             3

echo "[$(date)] [$CHROM] === STEP_C04d D17 wrapper ==="
$RSCRIPT 2e_ghsl_discovery/STEP_C04d_ghsl_d17_wrapper.R \
  --chrom          "$CHROM" \
  --localpca       "$GHSL_LOCALPCA_DIR/${CHROM}.ghsl_v6_localpca.rds" \
  --out_dir        "$GHSL_D17_DIR" \
  --d17_l1_script  "$D17_L1_SCRIPT" \
  --d17_l2_script  "$D17_L2_SCRIPT"

# To speed up D17 if needed, uncomment:
# --d17_l1_args "--boundary_grow_W_pct 0.005,0.02,0.05"

echo "[$(date)] [$CHROM] === export_ghsl_to_json_v3 ==="
$RSCRIPT 2e_ghsl_discovery/export_ghsl_to_json_v3.R \
  --chrom          "$CHROM" \
  --annot_rds      "$GHSL_V6_DIR/${CHROM}.ghsl_v6.annot.rds" \
  --persamp_rds    "$GHSL_V6_DIR/${CHROM}.ghsl_v6.per_sample.rds" \
  --karyo_rds      "$GHSL_V6_DIR/${CHROM}.ghsl_v6.karyotypes.rds" \
  --localpca_rds   "$GHSL_LOCALPCA_DIR/${CHROM}.ghsl_v6_localpca.rds" \
  --d17_l1_env     "$GHSL_D17_DIR/${CHROM}_ghsl_d17L1_envelopes.tsv" \
  --d17_l2_env     "$GHSL_D17_DIR/${CHROM}_ghsl_d17L2_envelopes.tsv" \
  --d17_l1_bnd     "$GHSL_D17_DIR/${CHROM}_ghsl_d17L1_boundaries.tsv" \
  --d17_l2_bnd     "$GHSL_D17_DIR/${CHROM}_ghsl_d17L2_boundaries.tsv" \
  --out_dir        "$JSON_OUT" \
  --primary_scale  s50 \
  --max_k          6

echo "[$(date)] [$CHROM] === DONE ==="
ls -lh "$JSON_OUT/${CHROM}/${CHROM}_phase2_ghsl.json"
