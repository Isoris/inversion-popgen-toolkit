#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=127
#SBATCH --mem=237G
#SBATCH -t 2:00:00
#SBATCH -J inv_C01b
#SBATCH -o logs/inv_C01b_%j.out
#SBATCH -e logs/inv_C01b_%j.err

# =============================================================================
# LAUNCH_C01b_snake1.slurm
#
# Multi-step launcher for Snake 1 pipeline:
#
#   --step all          Run cores → merge → diag (default)
#   --step cores        Run cores only (saves snake1_cores_*.rds)
#   --step merge        Run fuzzy merge only (reads cores RDS)
#   --step diag         Run diagnostic plots only
#   --step diag_precomp Run C01a precomp diagnostics only
#
# Combine steps with comma:
#   --step merge,diag   Run merge then diagnostics
#   --step cores,merge  Run cores then merge (no diag)
#
# Tag creates a subfolder for the output:
#   --step merge --tag relaxed_v2
#
# Usage:
#   sbatch LAUNCH_C01b_snake1.slurm                          # full pipeline
#   sbatch LAUNCH_C01b_snake1.slurm --step cores             # cores only
#   sbatch LAUNCH_C01b_snake1.slurm --step merge             # merge only
#   sbatch LAUNCH_C01b_snake1.slurm --step merge,diag        # merge + plots
#   sbatch LAUNCH_C01b_snake1.slurm --step diag              # plots only
#   sbatch LAUNCH_C01b_snake1.slurm --step diag_precomp      # C01a plots
#   sbatch LAUNCH_C01b_snake1.slurm --step merge --tag v2    # tagged merge
# =============================================================================

set -euo pipefail
source ~/.bashrc
mamba activate assembly

CONFIG="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_codebase_v8.5/00_inversion_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }

set -a
source "${CONFIG}"
set +a

inv_init_dirs

A2_SNAKES="${DISCOVERY2DIR}/snakes"
PRECOMP_DIR="${MDS_DIR}/snake_regions_multiscale/precomp"
OUTDIR="${MDS_DIR}/snake_regions_multiscale"

# ── Parse arguments ──
STEP="all"
TAG=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --step)  STEP="$2"; shift 2 ;;
    --tag)   TAG="$2"; shift 2 ;;
    *)       shift ;;
  esac
done

TAG_ARGS=""
EFFECTIVE_OUTDIR="${OUTDIR}"
if [[ -n "${TAG}" ]]; then
  TAG_ARGS="--tag ${TAG}"
  EFFECTIVE_OUTDIR="${OUTDIR}/${TAG}"
fi

DIAG_C01A="${EFFECTIVE_OUTDIR}/diag_C01a"
DIAG_C01B="${EFFECTIVE_OUTDIR}/diag_C01b"

echo "================================================================"
echo "  C01b: Snake 1 — step: ${STEP}"
echo "  Tag: ${TAG:-default}"
echo "  Precomp: ${PRECOMP_DIR}"
echo "  Output: ${EFFECTIVE_OUTDIR}"
echo "  Diag C01a: ${DIAG_C01A}"
echo "  Diag C01b: ${DIAG_C01B}"
echo "  JobID: ${SLURM_JOB_ID}"
echo "  Started: $(date)"
echo "================================================================"

# ── Helper functions ──
run_cores() {
  echo ""
  echo "[LAUNCH] ── Phase: CORES (S/M/L) ──"
  "${RSCRIPT_BIN}" "${A2_SNAKES}/STEP_C01b_1_cores.R" \
    "${PRECOMP_DIR}" "${OUTDIR}" ${TAG_ARGS}

  # Verify RDS files were saved
  local n_rds=$(ls "${EFFECTIVE_OUTDIR}"/snake1_cores_*.rds 2>/dev/null | wc -l)
  echo "[LAUNCH] Core RDS files saved: ${n_rds}"
  if [[ "${n_rds}" -eq 0 ]]; then
    echo "[WARN] No core RDS files found — merge step will fail"
  fi
}

run_merge() {
  echo ""
  echo "[LAUNCH] ── Phase: FUZZY MERGE (A/B/C) ──"

  # Check cores exist
  local n_rds=$(ls "${EFFECTIVE_OUTDIR}"/snake1_cores_*.rds 2>/dev/null | wc -l)
  if [[ "${n_rds}" -eq 0 ]]; then
    echo "[ERROR] No core RDS files in ${EFFECTIVE_OUTDIR}/"
    echo "[ERROR] Run --step cores first"
    return 1
  fi
  echo "[LAUNCH] Found ${n_rds} core RDS files"

  # v8.5: pass landscape dir if it exists
  local LANDSCAPE_OPTS=""
  local LANDSCAPE_DIR="${MDS_DIR}/snake_regions_multiscale/landscape"
  if [[ -d "${LANDSCAPE_DIR}" ]]; then
    LANDSCAPE_OPTS="--landscape_dir ${LANDSCAPE_DIR}"
    echo "[LAUNCH] Landscape data found: ${LANDSCAPE_DIR}"
  else
    echo "[LAUNCH] No landscape dir — running without 01C integration"
  fi

  "${RSCRIPT_BIN}" "${A2_SNAKES}/STEP_C01b_2_merge.R" \
    "${PRECOMP_DIR}" "${OUTDIR}" ${TAG_ARGS} ${LANDSCAPE_OPTS}
}

run_diag() {
  echo ""
  echo "[LAUNCH] ── Phase: DIAGNOSTICS (C01b) ──"
  "${RSCRIPT_BIN}" "${A2_SNAKES}/STEP_C01b_diag.R" \
    "${PRECOMP_DIR}" "${EFFECTIVE_OUTDIR}" "${DIAG_C01B}"
}

run_diag_precomp() {
  echo ""
  echo "[LAUNCH] -- Phase: DIAGNOSTICS (C01a precomp) --"
  "${RSCRIPT_BIN}" "${A2_SNAKES}/STEP_C01a_diag.R" \
    "${PRECOMP_DIR}" "${DIAG_C01A}"
}

run_network_diag() {
  echo ""
  echo "[LAUNCH] -- Phase: NETWORK DIAGNOSTICS (C01a2) --"
  local NETWORK_DIR="${EFFECTIVE_OUTDIR}/diag_network"
  local NET_ARGS=""
  # Auto-detect ancestry Q-matrix and relatedness file if they exist
  local Q_FILE="${MDS_DIR}/../07_round1_structure/best_seed_K8.qopt"
  local REL_FILE="${MDS_DIR}/../05_ngsrelate/catfish_226_for_natora.txt"
  if [[ -f "${Q_FILE}" ]]; then
    NET_ARGS="${NET_ARGS} --ancestry ${Q_FILE}"
  fi
  if [[ -f "${REL_FILE}" ]]; then
    NET_ARGS="${NET_ARGS} --relatedness ${REL_FILE}"
  fi
  "${RSCRIPT_BIN}" "${A2_SNAKES}/STEP_C01a2_network_diag.R" \
    "${PRECOMP_DIR}" "${NETWORK_DIR}" ${NET_ARGS}
}

# -- Dispatch --
# "all" expands to cores,merge,diag
if [[ "${STEP}" == "all" ]]; then
  STEP="cores,merge,diag"
fi

# Split comma-separated steps and run in order
IFS=',' read -ra STEPS <<< "${STEP}"
for s in "${STEPS[@]}"; do
  case "${s}" in
    cores)        run_cores ;;
    merge)        run_merge ;;
    diag)         run_diag ;;
    diag_precomp) run_diag_precomp ;;
    network_diag) run_network_diag ;;
    *)
      echo "[ERROR] Unknown step '${s}'. Use: cores, merge, diag, diag_precomp, network_diag, all"
      exit 1
      ;;
  esac
done

echo ""
echo "================================================================"
echo "  C01b steps '${STEP}' complete — $(date)"
echo "  Output: ${EFFECTIVE_OUTDIR}/"
echo "================================================================"
