#!/bin/bash
# =============================================================================
# install.sh — MODULE_QC_ShelfDiagnosis one-step setup
# =============================================================================
# Run once after git-cloning. Checks that required inputs exist at the paths in
# 00_config.sh, compiles Engine B if its source exists and its binary doesn't,
# and writes a local overlay config (config.local.sh — gitignored) that the
# user can edit without polluting the repo.
#
# Usage:
#   cd ${BASE}/inversion-popgen-toolkit/Modules/MODULE_QC_ShelfDiagnosis
#   bash install.sh
#   bash install.sh --compile-only        # just (re)compile the C++ engine
#   bash install.sh --check-only          # just verify paths, don't build
# =============================================================================
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---- Parse args --------------------------------------------------------------
COMPILE_ONLY=0
CHECK_ONLY=0
for a in "$@"; do
  case "$a" in
    --compile-only) COMPILE_ONLY=1 ;;
    --check-only)   CHECK_ONLY=1 ;;
    -h|--help)
      sed -n '3,25p' "$0" | sed 's/^# \{0,1\}//'
      exit 0 ;;
    *) echo "unknown arg: $a" >&2; exit 2 ;;
  esac
done

# ---- Source config -----------------------------------------------------------
# Default 00_config.sh ships with the repo. If config.local.sh exists, it
# overrides individual variables (because env vars set in config.local.sh are
# honoured by the `: "${VAR:=default}"` pattern in 00_config.sh).
if [[ -f "${here}/config.local.sh" ]]; then
  # shellcheck disable=SC1091
  source "${here}/config.local.sh"
fi
# shellcheck disable=SC1091
source "${here}/00_config.sh"

echo "== MODULE_QC_ShelfDiagnosis install =="
echo "  BASE              = ${BASE}"
echo "  BEAGLE_DIR        = ${BEAGLE_DIR}"
echo "  PRECOMP_DIR       = ${PRECOMP_DIR}"
echo "  HET_DIR           = ${HET_DIR}"
echo "  UNIFIED_ANCESTRY  = ${UNIFIED_ANCESTRY_DIR}"
echo "  LOCAL_Q_DIR       = ${LOCAL_Q_DIR}"
echo "  QC_OUT            = ${QC_OUT}"
echo ""

# ---- Path checks -------------------------------------------------------------
missing=0
check_path() {
  local label="$1" path="$2" hard="${3:-0}"
  if [[ -e "${path}" ]]; then
    printf "  [OK]    %-20s %s\n" "${label}" "${path}"
  else
    printf "  [MISS]  %-20s %s\n" "${label}" "${path}"
    if [[ "${hard}" == "1" ]]; then
      missing=$((missing + 1))
    fi
  fi
}

echo "Required inputs:"
check_path "BEAGLE_DIR"    "${BEAGLE_DIR}"        1
check_path "PRECOMP_DIR"   "${PRECOMP_DIR}"       1
check_path "chr.list"      "${BEAGLE_DIR}/chr.list" 1
echo ""
echo "Optional inputs (used by specific Q-steps):"
check_path "HET_DIR"            "${HET_DIR}"              0
check_path "THETA_MAIN_DIR"     "${THETA_MAIN_DIR}"       0
check_path "THETA_MULTI_DIR"    "${THETA_MULTI_DIR}"      0
check_path "MOSDEPTH_EXISTING"  "${MOSDEPTH_EXISTING}"    0
check_path "UNIFIED_ANCESTRY_DIR" "${UNIFIED_ANCESTRY_DIR}" 0
check_path "LOCAL_Q_DIR"        "${LOCAL_Q_DIR}"          0
echo ""

if [[ ${missing} -gt 0 ]]; then
  echo "ERROR: ${missing} required input(s) missing."
  echo "Edit config.local.sh (or 00_config.sh) to point to the correct locations."
  [[ "${CHECK_ONLY}" == "1" ]] || exit 1
fi

# ---- Compile Engine B (instant_q) if source exists and binary missing -------
INSTANT_Q_SRC="${UNIFIED_ANCESTRY_DIR}/src/instant_q.cpp"
INSTANT_Q_BIN="${UNIFIED_ANCESTRY_DIR}/src/instant_q"

if [[ "${CHECK_ONLY}" == "1" ]]; then
  echo ""
  echo "Engine B status:"
  if [[ -x "${INSTANT_Q_BIN}" ]]; then
    echo "  [OK]   instant_q built: ${INSTANT_Q_BIN}"
  elif [[ -f "${INSTANT_Q_SRC}" ]]; then
    echo "  [TODO] instant_q.cpp exists but not compiled. Run: bash install.sh --compile-only"
  else
    echo "  [MISS] no instant_q.cpp at ${INSTANT_Q_SRC}"
    echo "         Engine B tracks (Q06) will be skipped automatically."
  fi
else
  if [[ -f "${INSTANT_Q_SRC}" && ( ! -x "${INSTANT_Q_BIN}" || "${COMPILE_ONLY}" == "1" ) ]]; then
    echo ""
    echo "Compiling Engine B (instant_q)..."
    if [[ -f "${UNIFIED_ANCESTRY_DIR}/src/Makefile" ]]; then
      ( cd "${UNIFIED_ANCESTRY_DIR}/src" && make )
    else
      echo "  No Makefile, trying direct g++..."
      g++ -O3 -march=native -fopenmp -std=c++17 \
          -o "${INSTANT_Q_BIN}" "${INSTANT_Q_SRC}" -lz
    fi
    if [[ -x "${INSTANT_Q_BIN}" ]]; then
      echo "  [OK] Built: ${INSTANT_Q_BIN}"
    else
      echo "  [FAIL] Compile produced no binary. Check compiler/OpenMP/zlib."
      exit 1
    fi
  elif [[ -x "${INSTANT_Q_BIN}" ]]; then
    echo ""
    echo "Engine B binary present: ${INSTANT_Q_BIN}"
  else
    echo ""
    echo "Engine B source not found at ${INSTANT_Q_SRC} — Q06 tracks will be skipped."
  fi
fi

# ---- Compile Engine F (region_popstats) -------------------------------------
POPSTATS_DIR="${UNIFIED_ANCESTRY_DIR}/engines/fst_dxy"
POPSTATS_SRC="${POPSTATS_DIR}/region_popstats.c"

if [[ "${CHECK_ONLY}" == "1" ]]; then
  echo ""
  echo "Engine F status:"
  if [[ -x "${REGION_POPSTATS_BIN}" ]]; then
    echo "  [OK]   region_popstats built: ${REGION_POPSTATS_BIN}"
  elif [[ -f "${POPSTATS_SRC}" ]]; then
    echo "  [TODO] region_popstats.c exists but not compiled."
  else
    echo "  [MISS] no region_popstats.c at ${POPSTATS_SRC}"
    echo "         Engine F tracks (Q07) will be skipped automatically."
  fi
else
  if [[ -f "${POPSTATS_SRC}" && ( ! -x "${REGION_POPSTATS_BIN}" || "${COMPILE_ONLY}" == "1" ) ]]; then
    echo ""
    echo "Compiling Engine F (region_popstats)..."
    if [[ -f "${POPSTATS_DIR}/Makefile" ]]; then
      ( cd "${POPSTATS_DIR}" && make )
    else
      echo "  No Makefile, trying direct gcc..."
      gcc -O3 -march=native -fopenmp \
          -o "${REGION_POPSTATS_BIN}" "${POPSTATS_SRC}" -lz -lm
    fi
    if [[ -x "${REGION_POPSTATS_BIN}" ]]; then
      echo "  [OK] Built: ${REGION_POPSTATS_BIN}"
    else
      echo "  [FAIL] Compile produced no binary. Check compiler/OpenMP/zlib/libm."
      exit 1
    fi
  elif [[ -x "${REGION_POPSTATS_BIN}" ]]; then
    echo ""
    echo "Engine F binary present: ${REGION_POPSTATS_BIN}"
  else
    echo ""
    echo "Engine F source not found at ${POPSTATS_SRC} — Q07 tracks will be skipped."
  fi
fi

# ---- Create config.local.sh skeleton if absent ------------------------------
if [[ ! -f "${here}/config.local.sh" && "${CHECK_ONLY}" != "1" && "${COMPILE_ONLY}" != "1" ]]; then
  cat > "${here}/config.local.sh" <<EOF
# =============================================================================
# config.local.sh — per-user overrides for 00_config.sh (gitignored)
# =============================================================================
# Uncomment and edit any variable to override the default in 00_config.sh.
# This file is sourced before 00_config.sh so your overrides win.
# =============================================================================

# export BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
# export BEAGLE_DIR="\${BASE}/inversion_localpca_v7/02_snps_beagle"
# export PRECOMP_DIR="\${BASE}/inversion_localpca_v7/06_mds_candidates/snake_regions_multiscale/precomp"
# export HET_DIR="\${BASE}/het_roh"
# export MOSDEPTH_EXISTING="\${BASE}/het_roh/mosdepth_output"
# export UNIFIED_ANCESTRY_DIR="\${BASE}/unified_ancestry"
# export LOCAL_Q_DIR="\${UNIFIED_ANCESTRY_DIR}/local_Q"
# export QC_OUT="\${BASE}/MODULE_QC_ShelfDiagnosis/results"

# Per-run settings
# export BIN_MB=0.05
# export UNCERTAINTY_THRESH=0.9
# export THETA_SCALE=win50000.step10000
# export ANC_SCALES=1x,5x,10x
# export SMOOTH_WIN=5
EOF
  echo ""
  echo "Created ${here}/config.local.sh — edit to customize paths."
fi

# ---- Output directories ------------------------------------------------------
echo ""
echo "Creating output directories..."
mkdir -p "${QC_TRACKS}" "${QC_FIGS}" "${QC_LOGS}"
echo "  ${QC_TRACKS}"
echo "  ${QC_FIGS}"
echo "  ${QC_LOGS}"

echo ""
echo "== Install complete =="
echo "Next steps:"
echo "  1) Edit config.local.sh if any path is wrong."
echo "  2) Test a single chromosome:"
echo "       SHELF_START_MB=15 SHELF_END_MB=18 bash run_all.sh C_gar_LG28"
echo "  3) Full 28-chr run (SLURM):"
echo "       sbatch slurm/array_28chrom.sh"
