#!/bin/bash
# =============================================================================
# run_all.sh — MODULE_QC_ShelfDiagnosis driver
# =============================================================================
# Runs Q01 -> Q02 -> Q03 -> Q05 -> Q06 -> Q04 for one or all chromosomes.
#
# Usage:
#   bash run_all.sh C_gar_LG28
#   SHELF_START_MB=15 SHELF_END_MB=18 SMOOTH_WIN=5 bash run_all.sh C_gar_LG28
#   bash run_all.sh all
#
# Env overrides:
#   SKIP_Q01/02/03/04/05/06=1   skip individual steps
#   THETA_SCALE=win50000.step10000   theta scale for Q05
# =============================================================================
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/00_config.sh"

CHR="${1:-}"
[[ -z "${CHR}" ]] && qc_die "Usage: $0 <CHR|all>"

: "${SKIP_Q01:=0}"
: "${SKIP_Q02:=0}"
: "${SKIP_Q03:=0}"
: "${SKIP_Q04:=0}"
: "${SKIP_Q05:=0}"
: "${SKIP_Q06:=0}"
: "${SKIP_Q06_MULTI:=0}"
: "${SKIP_Q07:=0}"
: "${SKIP_Q08:=0}"
: "${THETA_SCALE:=win50000.step10000}"
: "${ANC_SCALES:=1x,5x,10x}"

run_chr() {
  local c="$1"
  qc_log "========= ${c} ========="
  [[ "${SKIP_Q01}" == "1" ]] || bash "${here}/STEP_Q01_snp_density.sh"       "${c}"
  [[ "${SKIP_Q02}" == "1" ]] || bash "${here}/STEP_Q02_beagle_uncertainty.sh" "${c}"
  [[ "${SKIP_Q03}" == "1" ]] || bash "${here}/STEP_Q03_coverage_tracks.sh"    "${c}"
  [[ "${SKIP_Q05}" == "1" ]] || bash "${here}/STEP_Q05_aggregate_theta.sh"    "${c}" "${THETA_SCALE}"
  [[ "${SKIP_Q06}" == "1" ]] || bash "${here}/STEP_Q06_ancestry_tracks.sh"    "${c}"
  [[ "${SKIP_Q06_MULTI}" == "1" ]] || bash "${here}/STEP_Q06_multiscale.sh"   "${c}" "${ANC_SCALES}"
  [[ "${SKIP_Q07}" == "1" ]] || bash "${here}/STEP_Q07_popstats.sh"           "${c}"
  [[ "${SKIP_Q04}" == "1" ]] || bash "${here}/STEP_Q04_plot_diagnostic.sh"    "${c}"
  # Q08 requires a shelf region
  if [[ "${SKIP_Q08}" != "1" && -n "${SHELF_START_MB:-}" && -n "${SHELF_END_MB:-}" ]]; then
    bash "${here}/STEP_Q08_shelf_heatmap.sh" "${c}"
  fi
}

if [[ "${CHR}" == "all" ]]; then
  while IFS= read -r c; do
    [[ -z "${c}" ]] && continue
    run_chr "${c}"
  done < "${BEAGLE_DIR}/chr.list"
else
  run_chr "${CHR}"
fi

qc_log "ALL DONE"
qc_log "  Tracks:  ${QC_TRACKS}"
qc_log "  Figures: ${QC_FIGS}"
qc_log "  Logs:    ${QC_LOGS}"
