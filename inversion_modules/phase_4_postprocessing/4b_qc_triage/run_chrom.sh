#!/bin/bash
# =============================================================================
# run_chrom.sh — single-chromosome driver
# =============================================================================
# Runs the shelf QC pipeline for ONE chromosome, producing all outputs:
#   Q01: SNP density tracks
#   Q02: BEAGLE uncertainty
#   Q03: coverage tracks (skipped if no mosdepth)
#   Q05: aggregated theta
#   Q06: ancestry tracks (single-scale + multi-scale precompute)
#   Q07: popstats invgt (k-means grouping + Fst/dXY via region_popstats)
#   Q04: diagnostic figure (BOTH WITH AND WITHOUT no-data strips)
#   Q08: shelf heatmap
#   Q09: gap characterization
#
# Skip flags via env: SKIP_Q01=1 SKIP_Q02=1 ... SKIP_Q09=1
# Shelf region env:   SHELF_START_MB=15 SHELF_END_MB=18
# Breakpoints env:    BP1_MB=15.115 BP2_MB=18.005
#
# Usage:
#   bash run_chrom.sh C_gar_LG28
#   SHELF_START_MB=15 SHELF_END_MB=18 BP1_MB=15.115 BP2_MB=18.005 \
#     bash run_chrom.sh C_gar_LG28
# =============================================================================
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/00_config.sh"

CHR="${1:-}"
[[ -z "${CHR}" ]] && qc_die "Usage: $0 <CHR>  (e.g. C_gar_LG28)"

# Shelf params default to whole chromosome if not given
SHELF_START_MB="${SHELF_START_MB:-0}"
SHELF_END_MB="${SHELF_END_MB:-1000}"
BP1_MB="${BP1_MB:-}"
BP2_MB="${BP2_MB:-}"
SMOOTH_WIN="${SMOOTH_WIN:-5}"
SNP_DENSITY_SCALE_KB="${SNP_DENSITY_SCALE_KB:-10}"

qc_log "=== run_chrom ${CHR} ==="
qc_log "shelf: ${SHELF_START_MB}-${SHELF_END_MB} Mb  breakpoints: ${BP1_MB:-none}/${BP2_MB:-none}"

run_step() {
  local step_name="$1"; shift
  local skip_var="SKIP_${step_name^^}"
  if [[ "${!skip_var:-}" == "1" ]]; then
    qc_log "SKIP ${step_name}"
    return 0
  fi
  qc_log "--- ${step_name} ---"
  "$@"
}

run_step Q01 bash "${here}/STEP_Q01_snp_density.sh" "${CHR}"
run_step Q02 bash "${here}/STEP_Q02_beagle_uncertainty.sh" "${CHR}"
run_step Q03 bash "${here}/STEP_Q03_coverage_tracks.sh" "${CHR}"
run_step Q05 bash "${here}/STEP_Q05_aggregate_theta.sh" "${CHR}"
run_step Q06 bash "${here}/STEP_Q06_precompute.sh" "${CHR}"
run_step Q06 bash "${here}/STEP_Q06_ancestry_tracks.sh" "${CHR}"
run_step Q06 bash "${here}/STEP_Q06_multiscale.sh" "${CHR}"
run_step Q07 bash "${here}/STEP_Q07_popstats.sh" "${CHR}"

# Q07b + Q07c: per-group ANGSD HWE (patched) + hobs_windower (Engine H, Merot)
run_step Q07B bash "${here}/STEP_Q07b_hobs_per_group.sh" "${CHR}"
run_step Q07C bash "${here}/STEP_Q07c_hobs_windower.sh" "${CHR}"

# Q04: render both variants (with + without no-data strips)
if [[ "${SKIP_Q04:-}" != "1" ]]; then
  qc_log "--- Q04 (two variants) ---"

  # With strips (diagnostic)
  OUT_SUFFIX="" Q04_NO_STRIPS="" \
  BREAKPOINT1_MB="${BP1_MB}" BREAKPOINT2_MB="${BP2_MB}" \
  SMOOTH_WIN="${SMOOTH_WIN}" SNP_DENSITY_SCALE_KB="${SNP_DENSITY_SCALE_KB}" \
    bash "${here}/STEP_Q04_plot_diagnostic.sh" "${CHR}"

  # Without strips (clean presentation)
  OUT_SUFFIX=".clean" Q04_NO_STRIPS="1" \
  BREAKPOINT1_MB="${BP1_MB}" BREAKPOINT2_MB="${BP2_MB}" \
  SMOOTH_WIN="${SMOOTH_WIN}" SNP_DENSITY_SCALE_KB="${SNP_DENSITY_SCALE_KB}" \
    bash "${here}/STEP_Q04_plot_diagnostic.sh" "${CHR}"
fi

run_step Q08 env SHELF_START_MB="${SHELF_START_MB}" SHELF_END_MB="${SHELF_END_MB}" \
              BP1_MB="${BP1_MB}" BP2_MB="${BP2_MB}" \
              bash "${here}/STEP_Q08_shelf_heatmap.sh" "${CHR}"

run_step Q09 bash "${here}/STEP_Q09_gap_characterization.sh" "${CHR}"

# Q10: bridge outputs into the four registries (sample / interval / results / evidence)
run_step Q10 env \
  SHELF_START_MB="${SHELF_START_MB}" SHELF_END_MB="${SHELF_END_MB}" \
  BP1_MB="${BP1_MB}" BP2_MB="${BP2_MB}" \
  SAMPLE_GROUP="${SAMPLE_GROUP:-all_226}" \
  METHOD_TAG="${METHOD_TAG:-phase_qc_shelf}" \
  CID_PREFIX="${CID_PREFIX:-inv}" \
  bash "${here}/STEP_Q10_register.sh" "${CHR}"

qc_log "=== run_chrom ${CHR} DONE ==="
qc_log "Figures:"
ls -la "${QC_FIGS}"/diagnostic."${CHR}"*.pdf 2>/dev/null | sed 's/^/  /'
ls -la "${QC_FIGS}"/shelf_heatmap."${CHR}"*.pdf 2>/dev/null | sed 's/^/  /'
qc_log "Gap table:"
ls -la "${QC_TRACKS}"/gap_features."${CHR}".tsv 2>/dev/null | sed 's/^/  /'
