#!/usr/bin/env bash
# =============================================================================
# LAUNCH_module2b_figures.sh — Run best-seed + merge tree + all visualization
#
# v9.0 REWIRED:
#   - Sources pipeline_bridge.sh for cross-module env vars
#   - Passes --bridge to A04 for registry auto-registration
#   - Initializes sample registry on first run
#
# Usage:
#   bash launchers/LAUNCH_module2b_figures.sh [--best-k 8]
# =============================================================================
set -euo pipefail

# Source module config
source "$(dirname "$0")/../00_module2b_config.sh"

# Source cross-module bridge (sets LOAD_BRIDGE, REGISTRY_DIR, etc.)
_BRIDGE_SH="${BASE}/inversion_codebase_v8.5/utils/pipeline_bridge.sh"
if [[ -f "$_BRIDGE_SH" ]]; then
  source "$_BRIDGE_SH"
else
  m2b_log "WARN: pipeline_bridge.sh not found at $_BRIDGE_SH"
fi

SCRIPT_DIR="$(cd "$(dirname "$0")/../scripts" && pwd)"
BEST_K=8

while [[ $# -gt 0 ]]; do
  case "$1" in
    --best-k) BEST_K="$2"; shift 2;;
    *) echo "Usage: $0 [--best-k K]"; exit 1;;
  esac
done

m2b_init_dirs

# ── Define panels to process ────────────────────────────────────────────────
PANELS=(
  "wholegenome  500   all     ${N_SAMPLES}  ${SAMPLE_LIST}"
  "wholegenome  1000  all     ${N_SAMPLES}  ${SAMPLE_LIST}"
  "wholegenome  500   pruned  ${N_PRUNED}   ${PRUNED_LIST}"
)

for PANEL_LINE in "${PANELS[@]}"; do
  read -r SCOPE THIN SSET NSAMP SFILE <<< "$PANEL_LINE"
  RUN_TAG="$(build_run_tag "$SCOPE" "$THIN" "$SSET" "$NSAMP")"
  RUN_DIR="${MODULE2B_RESULTS}/${RUN_TAG}"

  if [[ ! -d "$RUN_DIR" ]]; then
    m2b_log "SKIP ${RUN_TAG} (no results directory)"
    continue
  fi

  m2b_log "============================================"
  m2b_log "Processing: ${RUN_TAG}"
  m2b_log "============================================"

  # ── STEP_A04: Best-seed selection + registry registration ────────────────
  if [[ ! -s "${RUN_DIR}/${RUN_TAG}_best_seed_by_K.tsv" ]]; then
    m2b_log "STEP_A04: Best-seed selection"
    "${RSCRIPT_BIN}" "${SCRIPT_DIR}/STEP_A04_best_seed_by_K.R" \
      --run-dir "${RUN_DIR}" \
      --eval-dir "${MODULE2B_EVALADMIX}" \
      --sample-file "${SFILE}" \
      --file-prefix "${RUN_TAG}" \
      --k-min "${K_MIN}" --k-max "${K_MAX}" \
      --seeds "$(IFS=,; echo "${SEEDS[*]}")" \
      --palette-name "${PALETTE_NAME}" \
      --best-k "${BEST_K}" \
      --bridge "${LOAD_BRIDGE}" \
      2>&1 || m2b_log "WARN: STEP_A04 had issues for ${RUN_TAG}"
  else
    m2b_log "STEP_A04: SKIP (already done)"
  fi

  # ── STEP_A05: K-hierarchy merge tree ───────────────────────────────────────
  m2b_log "STEP_A05: K-hierarchy merge tree"
  "${RSCRIPT_BIN}" "${SCRIPT_DIR}/STEP_A05_k_hierarchy_tree.R" \
    --results-dir "${RUN_DIR}" \
    --sample-file "${SFILE}" \
    --file-prefix "${RUN_TAG}" \
    --out-dir "${RUN_DIR}" \
    --k-min "${K_MIN}" --k-max "${K_MAX}" \
    2>&1 || m2b_log "WARN: STEP_A05 had issues"

  # ── STEP_B01: Admixture faceted ────────────────────────────────────────────
  m2b_log "STEP_B01: Admixture barplots"
  "${RSCRIPT_BIN}" "${SCRIPT_DIR}/STEP_B01_admixture_faceted.R" \
    --results-dir "${RUN_DIR}" \
    --eval-dir "${MODULE2B_EVALADMIX}" \
    --relatedness-file "${RELATEDNESS_RES}" \
    --pruned-samples "${PRUNED_LIST}" \
    --theme-file "${THEME_FILE}" \
    --file-prefix "${RUN_TAG}" \
    --out-dir "${MODULE2B_FIGURES}" \
    2>&1 || m2b_log "WARN: STEP_B01 had issues"

  # ── STEP_B02: Kinship network ─────────────────────────────────────────────
  m2b_log "STEP_B02: Kinship network"
  "${RSCRIPT_BIN}" "${SCRIPT_DIR}/STEP_B02_kinship_network.R" \
    --relatedness-file "${RELATEDNESS_RES}" \
    --ancestry-file "${RUN_DIR}/${RUN_TAG}_sample_ancestry.tsv" \
    --best-k "${BEST_K}" \
    --theme-file "${THEME_FILE}" \
    --out-dir "${MODULE2B_FIGURES}" \
    --prefix "${RUN_TAG}_kinship_network" \
    2>&1 || m2b_log "WARN: STEP_B02 had issues"

  # ── STEP_B03: Kinship heatmap ─────────────────────────────────────────────
  m2b_log "STEP_B03: Kinship heatmap"
  "${RSCRIPT_BIN}" "${SCRIPT_DIR}/STEP_B03_kinship_heatmap.R" \
    --relatedness-file "${RELATEDNESS_RES}" \
    --ancestry-file "${RUN_DIR}/${RUN_TAG}_sample_ancestry.tsv" \
    --best-k "${BEST_K}" \
    --pruned-samples "${PRUNED_LIST}" \
    --theme-file "${THEME_FILE}" \
    --out-dir "${MODULE2B_FIGURES}" \
    --prefix "${RUN_TAG}_kinship_heatmap" \
    2>&1 || m2b_log "WARN: STEP_B03 had issues"

  # ── STEP_B04: PCA + kinship overlay ───────────────────────────────────────
  PCANGSD_COV="${MODULE2B_PCANGSD}/${SSET}${NSAMP}/pcangsd.cov"
  if [[ -s "$PCANGSD_COV" ]]; then
    m2b_log "STEP_B04: PCA + kinship overlay"
    "${RSCRIPT_BIN}" "${SCRIPT_DIR}/STEP_B04_pca_kinship_overlay.R" \
      --pcangsd-cov "${PCANGSD_COV}" \
      --relatedness-file "${RELATEDNESS_RES}" \
      --ancestry-file "${RUN_DIR}/${RUN_TAG}_sample_ancestry.tsv" \
      --best-k "${BEST_K}" \
      --sample-file "${SFILE}" \
      --theme-file "${THEME_FILE}" \
      --out-dir "${MODULE2B_FIGURES}" \
      --prefix "${RUN_TAG}_pca_kinship" \
      2>&1 || m2b_log "WARN: STEP_B04 had issues"
  else
    m2b_log "STEP_B04: SKIP (no PCAngsd cov at ${PCANGSD_COV})"
  fi

  # ── STEP_B05: K-hierarchy Sankey ──────────────────────────────────────────
  TREE_FILE="${RUN_DIR}/${RUN_TAG}_merge_tree.tsv"
  TRACK_FILE="${RUN_DIR}/${RUN_TAG}_component_tracking.tsv"
  if [[ -s "$TREE_FILE" && -s "$TRACK_FILE" ]]; then
    m2b_log "STEP_B05: K-hierarchy Sankey"
    "${RSCRIPT_BIN}" "${SCRIPT_DIR}/STEP_B05_k_hierarchy_sankey.R" \
      --tree-file "${TREE_FILE}" \
      --tracking-file "${TRACK_FILE}" \
      --best-k "${BEST_K}" \
      --theme-file "${THEME_FILE}" \
      --out-dir "${MODULE2B_FIGURES}" \
      --file-prefix "${RUN_TAG}" \
      2>&1 || m2b_log "WARN: STEP_B05 had issues"
  else
    m2b_log "STEP_B05: SKIP (no merge tree files)"
  fi

  # ── STEP_B06: Relatedness summary ─────────────────────────────────────────
  m2b_log "STEP_B06: Relatedness summary"
  "${RSCRIPT_BIN}" "${SCRIPT_DIR}/STEP_B06_relatedness_summary.R" \
    --relatedness-file "${RELATEDNESS_RES}" \
    --ancestry-file "${RUN_DIR}/${RUN_TAG}_sample_ancestry.tsv" \
    --pruned-samples "${PRUNED_LIST}" \
    --natora-summary "${NATORA_SUMMARY}" \
    --best-k "${BEST_K}" \
    --theme-file "${THEME_FILE}" \
    --out-dir "${MODULE2B_FIGURES}" \
    --prefix "${RUN_TAG}_relatedness_summary" \
    2>&1 || m2b_log "WARN: STEP_B06 had issues"

  # ── STEP_B07: evalAdmix diagnostics ───────────────────────────────────────
  m2b_log "STEP_B07: evalAdmix diagnostics"
  "${RSCRIPT_BIN}" "${SCRIPT_DIR}/STEP_B07_evaladmix_diagnostic.R" \
    --eval-dir "${MODULE2B_EVALADMIX}" \
    --file-prefix "${RUN_TAG}" \
    --sample-file "${SFILE}" \
    --theme-file "${THEME_FILE}" \
    --out-dir "${MODULE2B_FIGURES}" \
    --k-min "${K_MIN}" --k-max "${K_MAX}" \
    --grid-only \
    2>&1 || m2b_log "WARN: STEP_B07 had issues"

  m2b_log "DONE: ${RUN_TAG}"
  echo ""

done

m2b_log "============================================"
m2b_log "All figures in: ${MODULE2B_FIGURES}"
m2b_log "============================================"
ls -lh "${MODULE2B_FIGURES}"/*.pdf 2>/dev/null | head -30 || true
