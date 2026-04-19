#!/bin/bash
# =============================================================================
# STEP_Q02_beagle_uncertainty.sh
# =============================================================================
# Build per-chromosome tracks of BEAGLE genotype-posterior uncertainty, binned
# at BIN_MB. An "uncertain" site for a sample is one where the max of the three
# posteriors (P00, P01, P11) is below a threshold (default 0.9). In a clean
# ANGSD + BEAGLE call at decent coverage, uncertain rate should be small (<10%).
# A local spike in uncertain rate indicates the caller had no confidence —
# often a repeat/mapping issue or a coverage dropout.
#
# Input:  ${BEAGLE_DIR}/main_qcpass.<CHR>.beagle.gz
#         ${BEAGLE_DIR}/main_qcpass.<CHR>.pos.fixed (for genomic coords)
#
# Output: ${QC_TRACKS}/uncertainty.<CHR>.tsv
#   Columns: chrom bin_start_bp bin_end_bp bin_mid_mb n_sites
#            mean_max_post mean_uncertain_frac
#            cv_across_samples  n_samples_with_high_uncertain
#
#   - mean_max_post           : mean of (max posterior) across sites x samples
#   - mean_uncertain_frac     : fraction of (site x sample) with max_post < threshold
#   - cv_across_samples       : coefficient of variation of per-sample uncertain rate
#                               -> low CV: uniform; high CV: sample-specific
#   - n_samples_with_high_unc : count of samples whose uncertain rate > 2x the bin mean
#
# Usage:
#   bash STEP_Q02_beagle_uncertainty.sh <CHR>
#   bash STEP_Q02_beagle_uncertainty.sh all
# =============================================================================
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/00_config.sh"

CHR="${1:-}"
: "${UNCERTAINTY_THRESH:=0.9}"
[[ -z "${CHR}" ]] && qc_die "Usage: $0 <CHR|all>"

run_one() {
  local chr="$1"
  local beagle="${BEAGLE_DIR}/main_qcpass.${chr}.beagle.gz"
  local pos_file="${BEAGLE_DIR}/main_qcpass.${chr}.pos.fixed"
  [[ -f "${pos_file}" ]] || pos_file="${BEAGLE_DIR}/main_qcpass.${chr}.pos"
  [[ -f "${beagle}" ]] || { qc_log "SKIP ${chr}: no beagle.gz"; return 0; }
  [[ -f "${pos_file}" ]] || { qc_log "SKIP ${chr}: no pos file"; return 0; }

  local out="${QC_TRACKS}/uncertainty.${chr}.tsv"
  qc_log "Q02 ${chr}: streaming ${beagle} (thresh=${UNCERTAINTY_THRESH})"

  # Hand to R for efficient posterior parsing — pure awk gets unwieldy.
  ${RSCRIPT_BIN} --vanilla "${here}/R/q02_beagle_stream.R" \
      --beagle  "${beagle}" \
      --pos     "${pos_file}" \
      --chrom   "${chr}" \
      --bin_mb  "${BIN_MB}" \
      --thresh  "${UNCERTAINTY_THRESH}" \
      --out     "${out}" 2> "${QC_LOGS}/q02_${chr}.log"

  local n_bins
  n_bins=$(( $(wc -l < "${out}") - 1 ))
  qc_log "Q02 ${chr}: ${n_bins} bins written"
}

if [[ "${CHR}" == "all" ]]; then
  chr_list="${BEAGLE_DIR}/chr.list"
  while IFS= read -r c; do
    [[ -z "${c}" ]] && continue
    run_one "${c}"
  done < "${chr_list}"
else
  run_one "${CHR}"
fi

qc_log "Q02 DONE. Tracks in ${QC_TRACKS}/"
