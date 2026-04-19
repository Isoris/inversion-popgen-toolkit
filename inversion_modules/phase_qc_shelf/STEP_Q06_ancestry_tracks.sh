#!/bin/bash
# =============================================================================
# STEP_Q06_ancestry_tracks.sh
# =============================================================================
# Read the unified-ancestry Engine B caches (local_Q/<chr>.local_Q_summary.tsv.gz
# and optional local_Q/<chr>.local_Q_samples.tsv.gz) and emit two compact
# per-chromosome track TSVs for downstream use.
#
# Inputs (from unified_ancestry Engine B):
#   ${LOCAL_Q_DIR}/<CHR>.local_Q_summary.tsv.gz    (per-window means)
#   ${LOCAL_Q_DIR}/<CHR>.local_Q_samples.tsv.gz    (per-window x per-sample; optional)
#
# Outputs:
#   ${QC_TRACKS}/ancestry_window.<CHR>.tsv
#     cols: chrom  window_start_bp  window_end_bp  window_mid_mb
#           delta12  entropy  ena  maxQ_label  cv_delta12_across_samples
#
#   ${QC_TRACKS}/ancestry_sample.<CHR>.tsv.gz     (only if samples cache exists)
#     cols: sample  chrom  window_mid_bp  maxQ  maxQ_label  delta12  delta13  entropy  ena
#
# Usage:
#   bash STEP_Q06_ancestry_tracks.sh <CHR>
#   bash STEP_Q06_ancestry_tracks.sh all
#
# Env:
#   LOCAL_Q_DIR (default: ${BASE}/unified_ancestry/local_Q)
# =============================================================================
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/00_config.sh"

CHR="${1:-}"
: "${LOCAL_Q_DIR:=${BASE}/unified_ancestry/local_Q}"
[[ -z "${CHR}" ]] && qc_die "Usage: $0 <CHR|all>"

if [[ ! -d "${LOCAL_Q_DIR}" ]]; then
  qc_log "WARNING: LOCAL_Q_DIR not found: ${LOCAL_Q_DIR}"
  qc_log "  Set LOCAL_Q_DIR env var or check unified_ancestry install"
fi

run_one() {
  local chr="$1"
  local summary="${LOCAL_Q_DIR}/${chr}.local_Q_summary.tsv.gz"
  [[ -f "${summary}" ]] || summary="${LOCAL_Q_DIR}/${chr}.local_Q_summary.tsv"
  if [[ ! -f "${summary}" ]]; then
    qc_log "SKIP ${chr}: no local_Q summary in ${LOCAL_Q_DIR}"
    return 0
  fi
  local samples_file="${LOCAL_Q_DIR}/${chr}.local_Q_samples.tsv.gz"
  [[ -f "${samples_file}" ]] || samples_file="${LOCAL_Q_DIR}/${chr}.local_Q_samples.tsv"
  [[ -f "${samples_file}" ]] || samples_file=""

  local out_win="${QC_TRACKS}/ancestry_window.${chr}.tsv"
  local out_samp="${QC_TRACKS}/ancestry_sample.${chr}.tsv.gz"

  qc_log "Q06 ${chr}: summary=${summary}"
  [[ -n "${samples_file}" ]] && qc_log "Q06 ${chr}: samples=${samples_file}"

  ${RSCRIPT_BIN} --vanilla "${here}/R/q06_ancestry_tracks.R" \
      --summary "${summary}" \
      $( [[ -n "${samples_file}" ]] && echo "--samples ${samples_file}" ) \
      --chrom "${chr}" \
      --out_win "${out_win}" \
      --out_samp "${out_samp}" 2> "${QC_LOGS}/q06_${chr}.log"

  local n_bins
  n_bins=$(( $(wc -l < "${out_win}") - 1 ))
  qc_log "Q06 ${chr}: ${n_bins} window rows -> ${out_win}"
}

if [[ "${CHR}" == "all" ]]; then
  while IFS= read -r c; do
    [[ -z "${c}" ]] && continue
    run_one "${c}"
  done < "${BEAGLE_DIR}/chr.list"
else
  run_one "${CHR}"
fi

qc_log "Q06 DONE."
