#!/bin/bash
# =============================================================================
#  00_config.sh
#  -------------------------------------------------------------------
#  Central configuration and shared logging/preview helpers for phase_qc_shelf.
# =============================================================================
#
#  PURPOSE   Define paths, binary locations, runtime constants, and a small
#            set of shared logging helpers used by every STEP_Q* script.
#  SOURCED BY   Every STEP_Q*.sh, run_chrom.sh, run_all.sh
#  DEPENDS ON   config.local.sh (site-specific path overrides, optional)
#  STATUS    STABLE              LAST UPDATED   2026-04-20
#
# -------------------------------------------------------------------
#  How the helpers compose
# -------------------------------------------------------------------
#  Every step script should call:
#
#    source 00_config.sh                     # paths + helpers
#    qc_config_snapshot "Q01" BIN_MB ZCAT_BIN   # print env snapshot
#    qc_banner_open "Q01" "${chr}"           # open per-unit banner
#    qc_preview_file "input" "${path}"       # show head-3 of a file
#    ... work ...
#    qc_preview_file "output" "${out}"
#    qc_banner_close "Q01" "${chr}" "${sec}"
#
#  Log grep recipes (printed at end of every run via qc_log_grep_tips):
#    grep '🧾' log     -> config snapshots
#    grep '▶'  log     -> banners (one per unit)
#    grep '⚠️'  log     -> warnings
#    grep '❌'  log     -> errors
#    grep '📊' log     -> summaries
#    grep 'wall time' log -> timings
#    grep '📥\|📤' log -> inputs and outputs with their paths
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Paths (override in config.local.sh if the defaults don't match your site)
# -----------------------------------------------------------------------------
: "${BASE:=/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"
: "${BEAGLE_DIR:=${BASE}/inversion_localpca_v7/02_snps_beagle}"
: "${PRECOMP_DIR:=${BASE}/inversion_localpca_v7/06_mds_candidates/snake_regions_multiscale/precomp}"
: "${BAM_DIR:=${BASE}/MODULE_1_Reads}"
: "${BAM_GLOB:=*.bam}"
: "${MOSDEPTH_EXISTING:=${BASE}/het_roh/mosdepth_output}"
: "${HETDIR:=${BASE}/het_roh}"
: "${THETA_EXISTING:=${BASE}/het_roh/02_heterozygosity/03_theta}"
: "${THETA_MULTISCALE:=${BASE}/het_roh/02_heterozygosity/03_theta/multiscale}"
: "${UNIFIED_ANCESTRY_DIR:=${BASE}/inversion-popgen-toolkit/unified_ancestry}"
: "${LOCAL_Q_DIR:=${UNIFIED_ANCESTRY_DIR}/local_Q}"
: "${REGION_POPSTATS_BIN:=${UNIFIED_ANCESTRY_DIR}/engines/fst_dxy/region_popstats}"
: "${SAMPLES_IND:=${BASE}/het_roh/01_inputs_check/samples.ind}"

: "${QC_OUT:=${BASE}/inversion_localpca_v7/phase_qc_shelf_results}"
: "${QC_TRACKS:=${QC_OUT}/tracks}"
: "${QC_FIGS:=${QC_OUT}/figures}"
: "${QC_LOGS:=${QC_OUT}/logs}"
mkdir -p "${QC_TRACKS}" "${QC_FIGS}" "${QC_LOGS}"

# -----------------------------------------------------------------------------
# 2. Name translations (bridge older/parallel naming in other configs)
# -----------------------------------------------------------------------------
: "${SAMPLE_LIST:=${SAMPLES_IND}}"
: "${SAMPLE_LIST_POPSTATS:=${SAMPLES_IND}}"
: "${HET_DIR:=${HETDIR}}"

# -----------------------------------------------------------------------------
# 3. Binning constants + tool paths
# -----------------------------------------------------------------------------
: "${BIN_MB:=0.05}"
: "${RSCRIPT_BIN:=Rscript}"
: "${MOSDEPTH_BIN:=mosdepth}"
: "${SAMTOOLS_BIN:=samtools}"
: "${TABIX_BIN:=tabix}"
: "${BGZIP_BIN:=bgzip}"
: "${ZCAT_BIN:=zcat}"

# -----------------------------------------------------------------------------
# 4. Logging helpers
# -----------------------------------------------------------------------------
# qc_log "msg"     → timestamped info line
# qc_warn "msg"    → ⚠️-prefixed warning (greppable)
# qc_err "msg"     → ❌-prefixed error (greppable)
# qc_die "msg"     → ❌ error + exit 1
qc_log() {
  printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*"
}
qc_warn() {
  qc_log "⚠️  WARN: $*"
}
qc_err() {
  qc_log "❌ ERROR: $*"
}
qc_die() {
  qc_err "$*"
  exit 1
}

# -----------------------------------------------------------------------------
# 5. Banner helpers (per-unit open/close for a step)
# -----------------------------------------------------------------------------
# qc_banner_open  <step_label>  <unit_label>
#   prints a separator, a heading line, and a sub-heading with the unit.
# qc_banner_close <step_label>  <unit_label>  <seconds_elapsed>
#   prints a "done" marker with timing and a separator.
qc_banner_open() {
  local step="$1" unit="$2"
  qc_log "==============================================================="
  qc_log "▶  ${step}  ${0##*/}"
  qc_log "   unit: ${unit}"
  qc_log "==============================================================="
}
qc_banner_close() {
  local step="$1" unit="$2" secs="$3"
  qc_log "✅ ${step} ${unit}: done  (${secs}s wall time)"
  qc_log "==============================================================="
}

# -----------------------------------------------------------------------------
# 6. Config snapshot helper
# -----------------------------------------------------------------------------
# qc_config_snapshot <step_label> VAR1 VAR2 ...
#   prints the values of the named variables, plus standard BEAGLE_DIR /
#   QC_TRACKS / host / user lines. Call once at the top of each step.
qc_config_snapshot() {
  local step="$1"
  shift
  qc_log "🧾 ${step} config snapshot:"
  # Always-printed baseline
  qc_log "    BEAGLE_DIR  = ${BEAGLE_DIR:-<unset>}"
  qc_log "    QC_TRACKS   = ${QC_TRACKS:-<unset>}"
  qc_log "    QC_FIGS     = ${QC_FIGS:-<unset>}"
  qc_log "    QC_LOGS     = ${QC_LOGS:-<unset>}"
  # Step-specific extras
  for var in "$@"; do
    local val="${!var:-<unset>}"
    qc_log "    $(printf '%-11s' "${var}")= ${val}"
  done
  qc_log "    host        = $(hostname -s 2>/dev/null || echo unknown)"
  qc_log "    user        = ${USER:-unknown}"
  qc_log "    pwd         = $(pwd)"
  qc_log "    started     = $(date -u +%Y-%m-%dT%H:%M:%SZ)"
}

# -----------------------------------------------------------------------------
# 7. File preview helper
# -----------------------------------------------------------------------------
# qc_preview_file <label> <path> [n_lines=3]
#   emits a log line with the path + total line count, then the head-N lines
#   boxed with ASCII rails for clean grep/display.
qc_preview_file() {
  local label="$1"
  local path="$2"
  local n="${3:-3}"
  local arrow
  case "${label}" in
    input*|in)   arrow="📥 input " ;;
    output*|out) arrow="📤 output" ;;
    *)           arrow="📄 ${label}" ;;
  esac
  if [[ ! -f "${path}" ]]; then
    qc_warn "${arrow}  ${path}  (missing!)"
    return 1
  fi
  local reader="cat"
  [[ "${path}" == *.gz ]] && reader="${ZCAT_BIN}"
  local nlines
  nlines=$(${reader} "${path}" | wc -l)
  local size_kb
  size_kb=$(du -k "${path}" 2>/dev/null | awk '{print $1}')
  qc_log "${arrow}  ${path}  (${nlines} lines, ${size_kb:-?} kb)"
  qc_log "    |-- first ${n} lines --------------------------------"
  ${reader} "${path}" | head -n "${n}" | while IFS= read -r ln; do
    qc_log "    | ${ln}"
  done
  qc_log "    +----------------------------------------------------"
}

# -----------------------------------------------------------------------------
# 8. Log grep-tips helper (printed at end of long runs so logs are self-docu-
#    menting). Call this once from run_all.sh / run_chrom.sh after all steps.
# -----------------------------------------------------------------------------
qc_log_grep_tips() {
  qc_log "📚 log grep recipes:"
  qc_log "    grep '🧾' <log>       # config snapshots"
  qc_log "    grep '▶'  <log>       # per-unit banners"
  qc_log "    grep '⚠️'  <log>       # warnings"
  qc_log "    grep '❌' <log>       # errors"
  qc_log "    grep '📊' <log>       # summaries"
  qc_log "    grep '📥\\|📤' <log>   # inputs and outputs with paths"
  qc_log "    grep 'wall time' <log>  # per-unit timings"
  qc_log "    grep '✅\\|❌' <log>  # outcome (success or fail)"
}

# -----------------------------------------------------------------------------
# 9. Sanity check the environment at source time
# -----------------------------------------------------------------------------
_check_qc_config() {
  local missing=0
  if [[ ! -d "${BEAGLE_DIR}" ]]; then
    qc_warn "BEAGLE_DIR not a directory: ${BEAGLE_DIR}"
    missing=$((missing + 1))
  fi
  if [[ ! -d "${PRECOMP_DIR}" ]]; then
    qc_warn "PRECOMP_DIR not a directory: ${PRECOMP_DIR}"
    missing=$((missing + 1))
  fi
  if [[ ! -r "${SAMPLES_IND}" ]]; then
    qc_warn "SAMPLES_IND not readable: ${SAMPLES_IND}"
    missing=$((missing + 1))
  fi
  if [[ "${missing}" -gt 0 ]]; then
    qc_warn "${missing} paths unresolved. Step scripts may fail."
  fi
}
_check_qc_config

# -----------------------------------------------------------------------------
# Optional: site-local override file
# -----------------------------------------------------------------------------
local_conf="$(dirname "${BASH_SOURCE[0]}")/config.local.sh"
[[ -f "${local_conf}" ]] && source "${local_conf}"
