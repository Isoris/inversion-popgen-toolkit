#!/bin/bash
# =============================================================================
# STEP_Q05_aggregate_theta.sh
# =============================================================================
# Aggregate per-sample ANGSD thetaStat .pestPG files for one chromosome into a
# compact long-format TSV. Used as input for the scrubber JSON export (Q06)
# and for downstream plots (genotype-group heterozygosity curves).
#
# Inputs:
#   ${HET_DIR}/02_heterozygosity/03_theta/*.win<WIN>.step<STEP>.pestPG        (main)
#   ${HET_DIR}/02_heterozygosity/03_theta/multiscale/*.win<WIN>.step<STEP>.pestPG
#
# Output:
#   ${QC_TRACKS}/theta.<CHR>.<SCALE>.tsv.gz
#     columns: sample chrom win_center_bp tP nSites theta_pi_persite
#
# Usage:
#   bash STEP_Q05_aggregate_theta.sh <CHR> [SCALE]
#   bash STEP_Q05_aggregate_theta.sh C_gar_LG28 win50000.step10000
#   bash STEP_Q05_aggregate_theta.sh all win500000.step500000
#
# Default SCALE is win500000.step500000 (main non-multiscale files).
# For the scrubber we usually want win50000.step10000 (finer resolution).
# =============================================================================
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/00_config.sh"

CHR="${1:-}"
SCALE="${2:-win500000.step500000}"
: "${HET_DIR:=${BASE}/het_roh}"
[[ -z "${CHR}" ]] && qc_die "Usage: $0 <CHR|all> [SCALE=win500000.step500000]"

# Select directory based on scale
if [[ "${SCALE}" == "win500000.step500000" ]]; then
  THETA_DIR="${HET_DIR}/02_heterozygosity/03_theta"
else
  THETA_DIR="${HET_DIR}/02_heterozygosity/03_theta/multiscale"
fi
[[ -d "${THETA_DIR}" ]] || qc_die "Theta dir not found: ${THETA_DIR}"

run_one() {
  local chr="$1"
  local out="${QC_TRACKS}/theta.${chr}.${SCALE}.tsv.gz"
  qc_log "Q05 ${chr}: scale=${SCALE}  dir=${THETA_DIR}"

  ${RSCRIPT_BIN} --vanilla "${here}/R/q05_aggregate_theta.R" \
      --theta_dir "${THETA_DIR}" \
      --scale     "${SCALE}" \
      --chrom     "${chr}" \
      --out       "${out}" 2> "${QC_LOGS}/q05_${chr}.log"

  local n_rows
  n_rows=$(( $(${ZCAT_BIN} "${out}" | wc -l) - 1 ))
  qc_log "Q05 ${chr}: ${n_rows} rows written -> ${out}"
}

if [[ "${CHR}" == "all" ]]; then
  while IFS= read -r c; do
    [[ -z "${c}" ]] && continue
    run_one "${c}"
  done < "${BEAGLE_DIR}/chr.list"
else
  run_one "${CHR}"
fi

qc_log "Q05 DONE."
