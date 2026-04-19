#!/bin/bash
# =============================================================================
# STEP_Q06_precompute.sh
# =============================================================================
# Run Engine B (instant_q) for a chromosome at multiple BEAGLE resolutions
# and deposit the caches in scale-labeled subdirectories under ${LOCAL_Q_DIR}.
#
# Why: two distinct scientific views
#   - THIN BEAGLE (popstruct_thin/04_beagle_byRF_majmin/*.thin_10000.beagle.gz)
#     One SNP per ~10 kb. Engine B produces ~50-100 windows per chromosome at
#     the default 100-SNP window size. Coarse, low-noise, cohort-level view of
#     ancestry regime.
#   - DENSE BEAGLE (inversion_localpca_v7/02_snps_beagle/main_qcpass.*.beagle.gz)
#     Unthinned. Engine B produces thousands of windows matching the local-PCA
#     grid. Fine-resolution, higher-noise view — needed for seeing boundary
#     sharpness at inversion breakpoints.
#
# Engine B handles the F-matrix size mismatch internally: the F matrix was
# estimated from the thinned BEAGLE (~950k sites), but Engine B skips sites
# that aren't in the F matrix. So the dense BEAGLE contributes more positional
# coverage, the thin F provides the ancestry frequencies, and the intersection
# drives the local EM.
#
# Output layout (matches STEP_Q06_multiscale.sh reader):
#   ${LOCAL_Q_DIR}/scale_thin/<CHR>.local_Q_{summary,samples,meta}.tsv.gz
#   ${LOCAL_Q_DIR}/scale_dense/<CHR>.local_Q_{summary,samples,meta}.tsv.gz
#
# Usage:
#   bash STEP_Q06_precompute.sh <CHR>                 # both scales
#   bash STEP_Q06_precompute.sh <CHR> thin            # just thin
#   bash STEP_Q06_precompute.sh <CHR> dense           # just dense
#   bash STEP_Q06_precompute.sh <CHR> thin,dense      # explicit both
#
# Env overrides (via config.local.sh):
#   IQ_BEAGLE_THIN   (default: ${BASE}/popstruct_thin/04_beagle_byRF_majmin)
#   IQ_BEAGLE_DENSE  (default: ${BASE}/inversion_localpca_v7/02_snps_beagle)
#   IQ_WINDOW_SIZE   (default: 100)
#   IQ_WINDOW_STEP   (default: 20)
#   IQ_EM_ITER       (default: 100)
#   IQ_NCORES        (default: 4)
#   INSTANT_Q_BIN    (default: ${UNIFIED_ANCESTRY_DIR}/src/instant_q)
#   BEST_FOPT / BEST_QOPT  (from 00_ancestry_config.sh)
# =============================================================================
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/00_config.sh"

# Also source the ancestry config for BEST_FOPT / BEST_QOPT / IQ_WINDOW_SIZE etc.
# (These are not in phase_qc_shelf's own config since they belong to unified_ancestry)
ANC_CFG="${UNIFIED_ANCESTRY_DIR}/00_ancestry_config.sh"
if [[ -f "${ANC_CFG}" ]]; then
  # shellcheck disable=SC1090
  source "${ANC_CFG}"
fi

CHR="${1:-}"
SCALES_SPEC="${2:-thin,dense}"
[[ -z "${CHR}" ]] && qc_die "Usage: $0 <CHR> [thin|dense|thin,dense]"

# Resolve paths / binaries (with overrides honored)
: "${IQ_BEAGLE_THIN:=${BASE}/popstruct_thin/04_beagle_byRF_majmin}"
: "${IQ_BEAGLE_DENSE:=${BASE}/inversion_localpca_v7/02_snps_beagle}"
: "${IQ_WINDOW_SIZE:=100}"
: "${IQ_WINDOW_STEP:=20}"
: "${IQ_EM_ITER:=100}"
: "${IQ_NCORES:=4}"
: "${INSTANT_Q_BIN:=${UNIFIED_ANCESTRY_DIR}/src/instant_q}"

# Sanity
[[ -x "${INSTANT_Q_BIN}" ]]       || qc_die "instant_q binary not executable: ${INSTANT_Q_BIN}"
[[ -n "${BEST_FOPT:-}" && -f "${BEST_FOPT}"  ]] || qc_die "BEST_FOPT not found (check ${ANC_CFG}): ${BEST_FOPT:-unset}"
[[ -n "${BEST_QOPT:-}" && -f "${BEST_QOPT}"  ]] || qc_die "BEST_QOPT not found: ${BEST_QOPT:-unset}"

mkdir -p "${LOCAL_Q_DIR}"

# Map scale label -> BEAGLE resolver
find_beagle() {
  local label="$1" chr="$2" dir
  case "${label}" in
    thin)
      dir="${IQ_BEAGLE_THIN}"
      # e.g. catfish.C_gar_LG28.rf.txt.thin_10000.beagle.gz
      ls "${dir}"/*"${chr}"*thin*beagle.gz 2>/dev/null | head -1 ;;
    dense)
      dir="${IQ_BEAGLE_DENSE}"
      # e.g. main_qcpass.C_gar_LG28.beagle.gz
      ls "${dir}"/*"${chr}"*beagle.gz 2>/dev/null \
        | grep -v "thin" | head -1 ;;
    *)
      echo "" ;;
  esac
}

run_scale() {
  local scale="$1" chr="$2"
  local bgl; bgl=$(find_beagle "${scale}" "${chr}")

  if [[ -z "${bgl}" || ! -f "${bgl}" ]]; then
    qc_log "Q06-precomp ${chr}/${scale}: no BEAGLE found — skip"
    return 0
  fi

  local outdir="${LOCAL_Q_DIR}/scale_${scale}"
  mkdir -p "${outdir}"

  local summary="${outdir}/${chr}.local_Q_summary.tsv.gz"
  if [[ -f "${summary}" && "${FORCE:-0}" != "1" ]]; then
    qc_log "Q06-precomp ${chr}/${scale}: cache present (${summary}). Set FORCE=1 to rebuild. Skip."
    return 0
  fi

  qc_log "Q06-precomp ${chr}/${scale}: beagle=$(basename "${bgl}"), win=${IQ_WINDOW_SIZE}, step=${IQ_WINDOW_STEP}"

  local t_start; t_start=$(date +%s)
  "${INSTANT_Q_BIN}" \
    --beagle    "${bgl}" \
    --fopt      "${BEST_FOPT}" \
    --qinit     "${BEST_QOPT}" \
    --precompute \
    --outdir    "${outdir}" \
    --chr       "${chr}" \
    --window_size "${IQ_WINDOW_SIZE}" \
    --window_step "${IQ_WINDOW_STEP}" \
    --em_iter     "${IQ_EM_ITER}" \
    --ncores      "${IQ_NCORES}" \
    --sample_output
  local t_end; t_end=$(date +%s)

  # Gzip any TSVs the binary left behind
  for f in "${outdir}/${chr}".local_Q_*.tsv; do
    [[ -f "${f}" ]] && gzip -f "${f}"
  done

  qc_log "Q06-precomp ${chr}/${scale}: done in $((t_end - t_start))s"
  ls -la "${outdir}/${chr}".local_Q_* 2>/dev/null | awk '{print "  " $NF}' >&2
}

IFS=',' read -ra SCALES_ARR <<< "${SCALES_SPEC}"
for scale in "${SCALES_ARR[@]}"; do
  run_scale "${scale}" "${CHR}"
done

qc_log "Q06-precomp ${CHR}: all requested scales complete."
