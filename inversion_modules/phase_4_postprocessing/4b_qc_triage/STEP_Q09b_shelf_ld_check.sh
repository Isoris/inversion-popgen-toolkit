#!/bin/bash
# =============================================================================
# STEP_Q09b_shelf_ld_check.sh
# =============================================================================
# For a chromosome with a shelf detected, test whether the shelf is a single
# inversion or multiple overlapping arrangements. Computes between-bin
# sample-arrangement correlation and emits a diagnostic PDF + summary TSV.
#
# Usage:
#   SHELF_START_MB=15 SHELF_END_MB=18 \
#     bash STEP_Q09b_shelf_ld_check.sh C_gar_LG28
# =============================================================================
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/00_config.sh"

CHR="${1:-}"
[[ -z "${CHR}" ]] && qc_die "Usage: $0 <CHR>"
: "${SHELF_START_MB:?SHELF_START_MB required}"
: "${SHELF_END_MB:?SHELF_END_MB required}"
N_BINS="${N_BINS:-20}"

qc_log "Q09b ${CHR}: shelf LD structure check"

beagle="${BEAGLE_DIR}/main_qcpass.${CHR}.beagle.gz"
pos="${BEAGLE_DIR}/main_qcpass.${CHR}.pos.fixed"
invgt="${QC_TRACKS}/invgt_assignments.${CHR}.tsv"

for f in "${beagle}" "${pos}" "${invgt}"; do
  [[ -f "${f}" ]] || { qc_log "SKIP ${CHR}: missing ${f}"; exit 0; }
done

out_pdf="${QC_FIGS}/shelf_ld_check.${CHR}.pdf"
out_tsv="${QC_TRACKS}/shelf_ld_summary.${CHR}.tsv"

${RSCRIPT_BIN} --vanilla "${here}/R/q09b_shelf_ld_check.R" \
  --beagle "${beagle}" --pos "${pos}" --invgt "${invgt}" \
  --chrom "${CHR}" \
  --shelf_start_mb "${SHELF_START_MB}" --shelf_end_mb "${SHELF_END_MB}" \
  --sample_list "${SAMPLE_LIST_POPSTATS}" \
  --out_pdf "${out_pdf}" --out_tsv "${out_tsv}" \
  --n_bins "${N_BINS}" \
  2> "${QC_LOGS}/q09b_${CHR}.log"

[[ -f "${out_pdf}" ]] && qc_log "  -> ${out_pdf}"
[[ -f "${out_tsv}" ]] && qc_log "  -> ${out_tsv}"
qc_log "Q09b ${CHR} DONE."
