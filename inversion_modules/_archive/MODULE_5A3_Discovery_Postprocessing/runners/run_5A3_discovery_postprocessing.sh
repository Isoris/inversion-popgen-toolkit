#!/usr/bin/env bash
# =============================================================================
# run_5A3_discovery_postprocessing.sh — Runner for MODULE_5A3 (STEP11–13)
#
# Downstream annotation: theta overlap, regional PCA groups, combined panels.
#
# Usage:
#   bash run_5A3_discovery_postprocessing.sh               # all STEP11–13
#   bash run_5A3_discovery_postprocessing.sh --step 12     # STEP12 only
# =============================================================================

set -euo pipefail
source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MOD_DIR="$(dirname "${SCRIPT_DIR}")"
source "${MOD_DIR}/../../00_inversion_config.sh"

A3_STEPS="${MOD_DIR}/steps"

STEP_ONLY="${1:-}"

echo "================================================================"
echo "  MODULE_5A3 — Discovery Postprocessing (STEP11–13)"
echo "  Started: $(date)"
echo "================================================================"

# ── STEP11 ─────────────────────────────────────────────────────────────
if [[ -z "${STEP_ONLY}" || "${STEP_ONLY}" == "11" || "${STEP_ONLY}" == "--step" ]]; then
  [[ "${STEP_ONLY}" == "--step" ]] && STEP_ONLY="${2:-}"
  if [[ -z "${STEP_ONLY}" || "${STEP_ONLY}" == "11" ]]; then
    inv_log "STEP11: Overlap candidate regions with population theta"
    POP_TP="${BRIDGE_DIR}/population_mean_tP_windows.tsv.gz"
    "${RSCRIPT_BIN}" "${A3_STEPS}/STEP11_overlap_candidate_regions_with_theta_and_het.R" \
        "${CAND_FILE}" "${POP_TP}" "${BRIDGE_DIR}/inversion_localpca" \
        chrom WinStart WinStop tP
  fi
fi

# ── STEP12 ─────────────────────────────────────────────────────────────
if [[ -z "${STEP_ONLY}" || "${STEP_ONLY}" == "12" ]]; then
  inv_log "STEP12: Regional PCA (per candidate — provide CANDIDATE_ID and CHR)"
  echo "# Example:"
  echo "${RSCRIPT_BIN} ${A3_STEPS}/STEP12_candidate_region_pca_groups_and_plotC.R \\"
  echo "    ${CAND_FILE} \${CANDIDATE_ID} \\"
  echo "    ${DOSAGE_DIR}/\${CHR}.sites.tsv.gz \\"
  echo "    ${DOSAGE_DIR}/\${CHR}.dosage.tsv.gz \\"
  echo "    ${BRIDGE_DIR}/het_bridge.sample_tP_windows.tsv.gz \\"
  echo "    ${SAMPLES_IND} ${REGIONAL_DIR}/inversion_localpca 5"
fi

# ── STEP13 ─────────────────────────────────────────────────────────────
if [[ -z "${STEP_ONLY}" || "${STEP_ONLY}" == "13" ]]; then
  inv_log "STEP13: Combined panels"
  echo "# Example:"
  echo "${RSCRIPT_BIN} ${A3_STEPS}/STEP13_make_combined_panels_ABC.R \\"
  echo "    ${REGIONAL_DIR}/inversion_localpca \${CANDIDATE_ID} \\"
  echo "    ${COMBINED_DIR}/inversion_localpca"
fi

echo ""
echo "  MODULE_5A3 complete — $(date)"
