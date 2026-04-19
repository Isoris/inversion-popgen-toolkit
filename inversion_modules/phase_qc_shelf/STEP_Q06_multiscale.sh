#!/bin/bash
# =============================================================================
# STEP_Q06_multiscale.sh
# =============================================================================
# Run Engine B (instant_q) at multiple window scales and emit stacked per-chrom
# ancestry tracks. Then consolidate them into a single multi-scale TSV so the
# Q04 diagnostic plot and the scrubber can show all scales at once.
#
# The point of multiple scales: small windows (fast, noisy, label-switching
# risk) vs large windows (slow, stable, low-resolution). If Δ12 is stable
# across all scales at a region, the signal is robust. If it only appears at
# the finest scale, it's likely noise.
#
# Input:  Engine B outputs at multiple scales, already computed:
#   ${LOCAL_Q_DIR}/scale_1x/<CHR>.local_Q_summary.tsv.gz    (per-window)
#   ${LOCAL_Q_DIR}/scale_5x/<CHR>.local_Q_summary.tsv.gz    (per 5 windows)
#   ${LOCAL_Q_DIR}/scale_10x/<CHR>.local_Q_summary.tsv.gz   (per 10 windows)
#
# If scale_1x is missing but a flat <CHR>.local_Q_summary.tsv.gz exists at
# ${LOCAL_Q_DIR}/, treat it as scale_1x for backward compatibility.
#
# Outputs:
#   ${QC_TRACKS}/ancestry_window_multiscale.<CHR>.tsv
#     cols: chrom scale window_start_bp window_end_bp window_mid_mb
#           delta12 entropy ena maxQ_label cv_delta12_across_samples
#   (long-format: one row per (scale, window))
#
# Usage:
#   bash STEP_Q06_multiscale.sh <CHR> [scales]
#   bash STEP_Q06_multiscale.sh C_gar_LG28 "1x,5x,10x"
# =============================================================================
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/00_config.sh"

CHR="${1:-}"
SCALES_SPEC="${2:-1x,5x,10x}"
[[ -z "${CHR}" ]] && qc_die "Usage: $0 <CHR|all> [scales=1x,5x,10x]"

resolve_summary() {
  # Given scale string (e.g. "5x"), find the corresponding summary file
  local chr="$1" scale="$2"
  # Priority order for location
  for p in \
    "${LOCAL_Q_DIR}/scale_${scale}/${chr}.local_Q_summary.tsv.gz" \
    "${LOCAL_Q_DIR}/scale_${scale}/${chr}.local_Q_summary.tsv" \
    "${LOCAL_Q_DIR}/${scale}/${chr}.local_Q_summary.tsv.gz" \
    "${LOCAL_Q_DIR}/${scale}/${chr}.local_Q_summary.tsv"; do
    if [[ -f "${p}" ]]; then echo "${p}"; return; fi
  done
  # Backward-compat: for 1x only, fall back to flat cache
  if [[ "${scale}" == "1x" ]]; then
    for p in "${LOCAL_Q_DIR}/${chr}.local_Q_summary.tsv.gz" \
             "${LOCAL_Q_DIR}/${chr}.local_Q_summary.tsv"; do
      if [[ -f "${p}" ]]; then echo "${p}"; return; fi
    done
  fi
  echo ""
}
resolve_samples() {
  local chr="$1" scale="$2"
  for p in \
    "${LOCAL_Q_DIR}/scale_${scale}/${chr}.local_Q_samples.tsv.gz" \
    "${LOCAL_Q_DIR}/scale_${scale}/${chr}.local_Q_samples.tsv" \
    "${LOCAL_Q_DIR}/${scale}/${chr}.local_Q_samples.tsv.gz" \
    "${LOCAL_Q_DIR}/${scale}/${chr}.local_Q_samples.tsv"; do
    if [[ -f "${p}" ]]; then echo "${p}"; return; fi
  done
  if [[ "${scale}" == "1x" ]]; then
    for p in "${LOCAL_Q_DIR}/${chr}.local_Q_samples.tsv.gz" \
             "${LOCAL_Q_DIR}/${chr}.local_Q_samples.tsv"; do
      if [[ -f "${p}" ]]; then echo "${p}"; return; fi
    done
  fi
  echo ""
}

run_one() {
  local chr="$1"
  local out_combined="${QC_TRACKS}/ancestry_window_multiscale.${chr}.tsv"
  qc_log "Q06-multi ${chr}: scales=${SCALES_SPEC}"

  # Parse scales
  IFS=',' read -ra SCALES_ARR <<< "${SCALES_SPEC}"

  # First pass: run per-scale Q06 -> per-scale track files
  local per_scale_tracks=()
  for scale in "${SCALES_ARR[@]}"; do
    local summary samples
    summary=$(resolve_summary "${chr}" "${scale}")
    samples=$(resolve_samples "${chr}" "${scale}")
    if [[ -z "${summary}" ]]; then
      qc_log "  scale=${scale}: no summary file, skipping"
      continue
    fi
    local out_win="${QC_TRACKS}/ancestry_window.${chr}.scale_${scale}.tsv"
    local out_samp="${QC_TRACKS}/ancestry_sample.${chr}.scale_${scale}.tsv.gz"

    local args=(
      --summary "${summary}"
      --chrom   "${chr}"
      --out_win "${out_win}"
      --out_samp "${out_samp}"
    )
    [[ -n "${samples}" ]] && args+=( --samples "${samples}" )

    ${RSCRIPT_BIN} --vanilla "${here}/R/q06_ancestry_tracks.R" "${args[@]}" \
        2>> "${QC_LOGS}/q06_multi_${chr}.log"

    per_scale_tracks+=( "${scale}:${out_win}" )
    qc_log "  scale=${scale}: ${out_win}"
  done

  if [[ ${#per_scale_tracks[@]} -eq 0 ]]; then
    qc_log "Q06-multi ${chr}: no scales resolved, no output"
    return 0
  fi

  # Second pass: concatenate into a long-format multiscale TSV with scale col
  local first=1
  for entry in "${per_scale_tracks[@]}"; do
    local scale="${entry%%:*}"
    local path="${entry#*:}"
    if [[ ${first} -eq 1 ]]; then
      # Keep header, prepend "scale"
      head -1 "${path}" | awk -v OFS='\t' '{print "scale", $0}' > "${out_combined}"
      tail -n +2 "${path}" | awk -v OFS='\t' -v s="${scale}" '{print s, $0}' >> "${out_combined}"
      first=0
    else
      tail -n +2 "${path}" | awk -v OFS='\t' -v s="${scale}" '{print s, $0}' >> "${out_combined}"
    fi
  done

  local nrows
  nrows=$(( $(wc -l < "${out_combined}") - 1 ))
  qc_log "Q06-multi ${chr}: combined ${nrows} rows -> ${out_combined}"
}

if [[ "${CHR}" == "all" ]]; then
  while IFS= read -r c; do
    [[ -z "${c}" ]] && continue
    run_one "${c}"
  done < "${BEAGLE_DIR}/chr.list"
else
  run_one "${CHR}"
fi

qc_log "Q06-multi DONE."
