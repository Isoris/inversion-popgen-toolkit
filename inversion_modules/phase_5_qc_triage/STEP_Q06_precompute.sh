#!/bin/bash
# =============================================================================
# STEP_Q06_precompute.sh  v2  (multi-K × multi-scale)
# =============================================================================
# Run Engine B (instant_q) at any combination of BEAGLE resolution and K.
# Output: ${LOCAL_Q_DIR}/scale_<res>/K<NN>/<CHR>.local_Q_*.tsv.gz
#
# Also symlinks the canonical-K cache to the pre-v2 flat path so Q06 readers
# that don't know about K subdirs still work.
# =============================================================================
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/00_config.sh"

ANC_CFG="${UNIFIED_ANCESTRY_DIR}/00_ancestry_config.sh"
if [[ -f "${ANC_CFG}" ]]; then
  _preserve_bin="${INSTANT_Q_BIN:-}"
  _preserve_local_q="${LOCAL_Q_DIR:-}"
  # shellcheck disable=SC1090
  source "${ANC_CFG}"
  [[ -n "${_preserve_bin}"      ]] && INSTANT_Q_BIN="${_preserve_bin}"
  [[ -n "${_preserve_local_q}"  ]] && LOCAL_Q_DIR="${_preserve_local_q}"
  unset _preserve_bin _preserve_local_q
fi

CHR="${1:-}"
SCALES_SPEC="${2:-${Q06_SCALES:-thin,dense}}"
K_SPEC="${3:-${Q06_K_SWEEP:-8}}"
[[ -z "${CHR}" ]] && qc_die "Usage: $0 <CHR> [scales=thin,dense] [K_list=8]"

: "${IQ_BEAGLE_THIN:=${BASE}/popstruct_thin/04_beagle_byRF_majmin}"
: "${IQ_BEAGLE_DENSE:=${BASE}/inversion_localpca_v7/02_snps_beagle}"
: "${IQ_WINDOW_SIZE:=100}"
: "${IQ_WINDOW_STEP:=20}"
: "${IQ_EM_ITER:=100}"
: "${IQ_NCORES:=4}"
: "${INSTANT_Q_BIN:=${UNIFIED_ANCESTRY_DIR}/src/instant_q}"
: "${NGSADMIX_RUNS_DIR:=${BASE}/popstruct_thin/05_ngsadmix_global/runs_thin500}"
: "${NGSADMIX_PREFIX:=thin500_K}"
: "${NGSADMIX_SUFFIX:=_best}"
: "${CANONICAL_K:=8}"

[[ -x "${INSTANT_Q_BIN}" ]]     || qc_die "instant_q binary not executable: ${INSTANT_Q_BIN}"
[[ -d "${NGSADMIX_RUNS_DIR}" ]] || qc_die "NGSadmix runs dir not found: ${NGSADMIX_RUNS_DIR}"
mkdir -p "${LOCAL_Q_DIR}"

find_beagle() {
  local label="$1" chr="$2"
  case "${label}" in
    thin)  ls "${IQ_BEAGLE_THIN}"/*"${chr}"*thin*beagle.gz 2>/dev/null | head -1 ;;
    dense) ls "${IQ_BEAGLE_DENSE}"/*"${chr}"*beagle.gz 2>/dev/null | grep -v thin | head -1 ;;
    *) echo "" ;;
  esac
}

find_ngsadmix_for_k() {
  local K="$1" kk
  kk=$(printf "%02d" "${K}")
  local prefix="${NGSADMIX_RUNS_DIR}/${NGSADMIX_PREFIX}${kk}${NGSADMIX_SUFFIX}"
  local fopt="${prefix}.fopt.gz"
  local qopt="${prefix}.qopt"
  if [[ -f "${fopt}" && -f "${qopt}" ]]; then
    echo "${fopt}|${qopt}"
  else
    echo ""
  fi
}

symlink_canonical() {
  local outdir="$1" chr="$2" scale="$3" K="$4"
  if [[ "${K}" == "${CANONICAL_K}" ]]; then
    local parent="${LOCAL_Q_DIR}/scale_${scale}"
    local kk; kk=$(printf "%02d" "${K}")
    for suffix in summary samples meta; do
      local src="${outdir}/${chr}.local_Q_${suffix}.tsv.gz"
      local dst="${parent}/${chr}.local_Q_${suffix}.tsv.gz"
      if [[ -f "${src}" ]]; then
        ln -sfn "K${kk}/$(basename "${src}")" "${dst}"
      fi
    done
  fi
}

run_one() {
  local scale="$1" chr="$2" K="$3"
  local bgl; bgl=$(find_beagle "${scale}" "${chr}")
  if [[ -z "${bgl}" || ! -f "${bgl}" ]]; then
    qc_log "Q06-precomp ${chr}/${scale}/K${K}: no BEAGLE found — skip"
    return 0
  fi
  local fq_pair; fq_pair=$(find_ngsadmix_for_k "${K}")
  if [[ -z "${fq_pair}" ]]; then
    qc_log "Q06-precomp ${chr}/${scale}/K${K}: NGSadmix F/Q for K=${K} not found in ${NGSADMIX_RUNS_DIR}"
    return 0
  fi
  local fopt="${fq_pair%|*}"
  local qopt="${fq_pair#*|}"

  local kk; kk=$(printf "%02d" "${K}")
  local outdir="${LOCAL_Q_DIR}/scale_${scale}/K${kk}"
  mkdir -p "${outdir}"

  local summary="${outdir}/${chr}.local_Q_summary.tsv.gz"
  if [[ -f "${summary}" && "${FORCE:-0}" != "1" ]]; then
    qc_log "Q06-precomp ${chr}/${scale}/K${K}: cache present. Skip."
    symlink_canonical "${outdir}" "${chr}" "${scale}" "${K}"
    return 0
  fi

  qc_log "Q06-precomp ${chr}/${scale}/K${K}: beagle=$(basename "${bgl}") fopt=$(basename "${fopt}")"
  local t_start; t_start=$(date +%s)
  "${INSTANT_Q_BIN}" \
    --beagle "${bgl}" --fopt "${fopt}" --qinit "${qopt}" \
    --precompute --outdir "${outdir}" --chr "${chr}" \
    --window_size "${IQ_WINDOW_SIZE}" --window_step "${IQ_WINDOW_STEP}" \
    --em_iter "${IQ_EM_ITER}" --ncores "${IQ_NCORES}" --sample_output
  local t_end; t_end=$(date +%s)

  for f in "${outdir}/${chr}".local_Q_*.tsv; do
    [[ -f "${f}" ]] && gzip -f "${f}"
  done
  qc_log "Q06-precomp ${chr}/${scale}/K${K}: done in $((t_end - t_start))s"

  symlink_canonical "${outdir}" "${chr}" "${scale}" "${K}"
}

IFS=',' read -ra SCALES_ARR <<< "${SCALES_SPEC}"
IFS=',' read -ra K_ARR      <<< "${K_SPEC}"
qc_log "Q06-precomp ${CHR}: scales=[${SCALES_ARR[*]}] K=[${K_ARR[*]}]"
for scale in "${SCALES_ARR[@]}"; do
  for K in "${K_ARR[@]}"; do
    scale="${scale// /}"; K="${K// /}"
    [[ -z "${scale}" || -z "${K}" ]] && continue
    run_one "${scale}" "${CHR}" "${K}"
  done
done
qc_log "Q06-precomp ${CHR}: all (scale,K) pairs complete."
