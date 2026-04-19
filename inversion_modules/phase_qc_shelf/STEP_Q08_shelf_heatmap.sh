#!/bin/bash
# =============================================================================
# STEP_Q08_shelf_heatmap.sh
# =============================================================================
# Render sample x SNP genotype (dosage) heatmap for the shelf region.
#
# This implements Image 1 Panel C: samples ordered by PC1 (within the shelf),
# coloured by expected dosage (0=REF, 1=HET, 2=INV) in a block-block-block
# pattern if real inversion, noise if artifact.
#
# Produces a standalone PDF:
#   ${QC_FIGS}/shelf_heatmap.<CHR>.pdf
#
# Usage:
#   SHELF_START_MB=15 SHELF_END_MB=18 bash STEP_Q08_shelf_heatmap.sh C_gar_LG28
# =============================================================================
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/00_config.sh"

CHR="${1:-}"
[[ -z "${CHR}" ]] && qc_die "Usage: $0 <CHR>"
[[ -z "${SHELF_START_MB:-}" || -z "${SHELF_END_MB:-}" ]] && \
  qc_die "SHELF_START_MB and SHELF_END_MB required"

run_one() {
  local chr="$1"
  qc_log "Q08 ${chr}: shelf genotype heatmap"

  local beagle="${BEAGLE_DIR}/main_qcpass.${chr}.beagle.gz"
  [[ -f "${beagle}" ]] || { qc_log "SKIP ${chr}: no beagle"; return 0; }
  local pos="${BEAGLE_DIR}/main_qcpass.${chr}.pos.fixed"
  [[ -f "${pos}" ]]    || { qc_log "SKIP ${chr}: no pos file"; return 0; }
  local precomp="${PRECOMP_DIR}/${chr}.precomp.rds"
  [[ -f "${precomp}" ]] || precomp="${PRECOMP_DIR}/main_qcpass.${chr}.precomp.rds"
  [[ -f "${precomp}" ]] || { qc_log "SKIP ${chr}: no precomp"; return 0; }

  local out_pdf="${QC_FIGS}/shelf_heatmap.${chr}.pdf"
  local out_dir="${QC_FIGS}/shelf_heatmap.${chr}"
  local invgt_assign="${QC_TRACKS}/invgt_assignments.${chr}.tsv"

  local args=(
    --beagle  "${beagle}"
    --pos     "${pos}"
    --precomp "${precomp}"
    --chrom   "${chr}"
    --shelf_start_mb "${SHELF_START_MB}"
    --shelf_end_mb   "${SHELF_END_MB}"
    --sample_list    "${SAMPLE_LIST_POPSTATS}"
    --out     "${out_pdf}"
    --out_dir "${out_dir}"
  )
  [[ -f "${invgt_assign}" ]] && args+=( --invgt_assign "${invgt_assign}" )

  ${RSCRIPT_BIN} --vanilla "${here}/R/q08_shelf_heatmap.R" "${args[@]}" \
    2> "${QC_LOGS}/q08_${chr}.log"

  [[ -f "${out_pdf}" ]] && qc_log "  -> ${out_pdf}"
}

run_one "${CHR}"
qc_log "Q08 DONE."
