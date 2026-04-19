#!/bin/bash
# =============================================================================
# STEP_Q03_coverage_tracks.sh
# =============================================================================
# Build per-chromosome coverage tracks binned at BIN_MB, with per-sample
# coverage information collapsed into:
#   - mean coverage across all samples
#   - CV across samples (uniform drop vs sample-specific drop)
#   - count of samples with coverage < 50% of chromosome-median
#
# Signal: uniform coverage drop across all samples = repeat/mappability issue
#   (affects all reads equally). Sample-specific drop = batch effect or a
#   subset of samples carrying a deletion. Both interpretations matter for
#   diagnosing the Z plateau.
#
# Two modes:
#   (A) REUSE_MOSDEPTH=1 (default): parse existing ${MOSDEPTH_EXISTING}/*regions*.bed.gz
#       Expects mosdepth already ran with --by BIN_BP or similar.
#   (B) RUN_MOSDEPTH=1: run mosdepth on BAMs to generate the tracks.
#
# Output: ${QC_TRACKS}/coverage.<CHR>.tsv
#   Columns: chrom bin_start_bp bin_end_bp bin_mid_mb
#            mean_cov  median_cov  cv_across_samples  n_samples_low_cov
# =============================================================================
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/00_config.sh"

CHR="${1:-}"
: "${REUSE_MOSDEPTH:=1}"
: "${RUN_MOSDEPTH:=0}"
: "${MOSDEPTH_WIN_BP:=$(awk "BEGIN{printf \"%d\", ${BIN_MB}*1e6}")}"
[[ -z "${CHR}" ]] && qc_die "Usage: $0 <CHR|all>"

# If user wants fresh mosdepth runs
run_mosdepth_all() {
  local chr="$1"
  local outdir="${QC_OUT}/mosdepth"
  mkdir -p "${outdir}"
  local bams=( "${BAM_DIR}"/${BAM_PATTERN} )
  [[ -f "${bams[0]}" ]] || qc_die "No BAMs matched at ${BAM_DIR}/${BAM_PATTERN}"
  qc_log "Q03 ${chr}: running mosdepth on ${#bams[@]} BAMs (window=${MOSDEPTH_WIN_BP} bp)"

  for bam in "${bams[@]}"; do
    local stem
    stem=$(basename "${bam}" .bam)
    local prefix="${outdir}/${stem}.${chr}"
    if [[ -f "${prefix}.regions.bed.gz" ]]; then
      qc_log "  skip ${stem}: regions.bed.gz exists"
      continue
    fi
    ${MOSDEPTH_BIN} \
      --no-per-base \
      --by "${MOSDEPTH_WIN_BP}" \
      --chrom "${chr}" \
      --threads 2 \
      "${prefix}" "${bam}" \
      2>> "${QC_LOGS}/q03_mosdepth.${chr}.log"
  done
}

build_track() {
  local chr="$1"
  local search_dir="${MOSDEPTH_EXISTING}"
  [[ "${RUN_MOSDEPTH}" == "1" ]] && search_dir="${QC_OUT}/mosdepth"
  [[ -d "${search_dir}" ]] || qc_die "No mosdepth dir: ${search_dir}"

  # Collect all regions bed files that have this chrom
  local files=( )
  while IFS= read -r -d '' f; do
    files+=( "${f}" )
  done < <(find "${search_dir}" -maxdepth 2 -name "*regions*.bed.gz" -print0 2>/dev/null)

  [[ ${#files[@]} -gt 0 ]] || { qc_log "SKIP ${chr}: no regions.bed.gz in ${search_dir}"; return 0; }
  qc_log "Q03 ${chr}: collapsing ${#files[@]} per-sample files"

  local out="${QC_TRACKS}/coverage.${chr}.tsv"
  # Delegate the per-bin aggregation to R (cleaner for CV/median across samples)
  ${RSCRIPT_BIN} --vanilla "${here}/R/q03_coverage_collapse.R" \
      --files_glob "${search_dir}" \
      --chrom "${chr}" \
      --bin_mb "${BIN_MB}" \
      --out "${out}" 2> "${QC_LOGS}/q03_${chr}.log"

  local n_bins
  n_bins=$(( $(wc -l < "${out}") - 1 ))
  qc_log "Q03 ${chr}: ${n_bins} bins written"
}

if [[ "${RUN_MOSDEPTH}" == "1" ]]; then
  if [[ "${CHR}" == "all" ]]; then
    while IFS= read -r c; do
      [[ -z "${c}" ]] && continue
      run_mosdepth_all "${c}"
    done < "${BEAGLE_DIR}/chr.list"
  else
    run_mosdepth_all "${CHR}"
  fi
fi

if [[ "${CHR}" == "all" ]]; then
  while IFS= read -r c; do
    [[ -z "${c}" ]] && continue
    build_track "${c}"
  done < "${BEAGLE_DIR}/chr.list"
else
  build_track "${CHR}"
fi

qc_log "Q03 DONE. Tracks in ${QC_TRACKS}/"
