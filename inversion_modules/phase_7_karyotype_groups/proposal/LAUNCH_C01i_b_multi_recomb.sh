#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -A lt200308
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -t 3:00:00
#SBATCH -J inv_4b2_recomb
#SBATCH -o logs/inv_4b2_recomb_%j.out
#SBATCH -e logs/inv_4b2_recomb_%j.err

# =============================================================================
# LAUNCH_C01i_b_multi_recomb.sh
#
# Phase 4b.2 — recombinant call gate.  Combines Regime signal (from C01j) with
# GHSL signal (from C04b). Emits per-candidate recombinant_map Tier-2 block;
# classifies each sample as RECOMBINANT_GC / RECOMBINANT_DCO / RECOMBINANT /
# NOT_RECOMB per the chat-12 gate.
#
# Called by: orchestrator/run_phase4b.sh (depends on 4b.1 decompose).
#
# Upstream:    C01i_decompose out (internal_dynamics block)
#              C01j regime_memberships.tsv.gz
#              C04b GHSL per-candidate output (optional but recommended)
# Downstream:  C01i_d_seal (reads recomb_map for final class resolution)
#
# Chat-16 sanity guards: MAX_N_SEGMENTS_HARD=10, MAX_N_SEGMENTS_SANE=6,
# MIN_SEGMENT_BP=10_000. Watch for loud REFUSING warnings.
# =============================================================================

set -euo pipefail
source ~/.bashrc
mamba activate assembly

CONFIG="${CONFIG:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_modules/00_inversion_config.sh}"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
set -a
source "${CONFIG}"
set +a

BRIDGE_SH="${BRIDGE_SH:-${BASE}/utils/pipeline_bridge.sh}"
[[ -f "${BRIDGE_SH}" ]] && source "${BRIDGE_SH}"

inv_init_dirs
mkdir -p logs

# ── Args ─────────────────────────────────────────────────────────────────
CANDIDATES="${CANDIDATES:-${SCORING_DIR}/pass1/candidate_scores.tsv.gz}"
DECOMP_OUT="${DECOMP_OUT:-${MDS_DIR}/snake_regions_multiscale/phase4b_out/decompose}"
RECOMB_OUT="${RECOMB_OUT:-${MDS_DIR}/snake_regions_multiscale/phase4b_out/multi_recomb}"
REGIME_MEMB="${REGIME_MEMB:-${REGIME_DIR}/regime_memberships.tsv.gz}"
GHSL_PER_CAND_DIR="${GHSL_PER_CAND_DIR:-${GHSL_DIR}/per_candidate}"
GC_DIR="${GC_DIR:-}"
TIER_MAX="${TIER_MAX:-3}"
MIN_DEV_FRAC="${MIN_DEV_FRAC:-0.15}"

mkdir -p "${RECOMB_OUT}"

[[ -f "${CANDIDATES}" ]] || { echo "[C01i_b] ERROR: missing ${CANDIDATES}" >&2; exit 2; }
[[ -d "${DECOMP_OUT}" ]]  || { echo "[C01i_b] ERROR: decomp dir missing — ${DECOMP_OUT}" >&2; exit 2; }
[[ -f "${REGIME_MEMB}" ]] || {
  echo "[C01i_b] ERROR: regime memberships required — ${REGIME_MEMB}" >&2
  echo "[C01i_b]         C01j regime engine must run first." >&2
  exit 2;
}

GHSL_OPTS=""
if [[ -d "${GHSL_PER_CAND_DIR}" ]]; then
  GHSL_OPTS="--ghsl_dir ${GHSL_PER_CAND_DIR}"
else
  echo "[C01i_b] WARN: no per-candidate GHSL dir — R-only fallback gate" >&2
fi

GC_OPTS=""
if [[ -n "${GC_DIR}" && -d "${GC_DIR}" ]]; then
  GC_OPTS="--gc_dir ${GC_DIR}"
fi

echo "================================================================"
echo "  4b.2 / C01i multi_recomb"
echo "  Candidates:   ${CANDIDATES}"
echo "  Decomp dir:   ${DECOMP_OUT}"
echo "  Regime memb:  ${REGIME_MEMB}"
echo "  GHSL dir:     ${GHSL_PER_CAND_DIR:-(skipped)}"
echo "  GC dir:       ${GC_DIR:-(skipped)}"
echo "  Out dir:      ${RECOMB_OUT}"
echo "  Tier max:     ${TIER_MAX}"
echo "  Min dev frac: ${MIN_DEV_FRAC}"
echo "  JobID:        ${SLURM_JOB_ID}"
echo "  Started:      $(date)"
echo "================================================================"

"${RSCRIPT_BIN}" \
  "${BASE}/inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_b_multi_recomb.R" \
  --candidates   "${CANDIDATES}" \
  --regime_memb  "${REGIME_MEMB}" \
  --decomp_dir   "${DECOMP_OUT}" \
  --outdir       "${RECOMB_OUT}" \
  --tier_max     "${TIER_MAX}" \
  --min_dev_frac "${MIN_DEV_FRAC}" \
  ${GHSL_OPTS} ${GC_OPTS}

echo "================================================================"
echo "  4b.2 complete — $(date)"
echo "  Recomb out: ${RECOMB_OUT}/"
echo "  Next (4b.4): sbatch LAUNCH_C01i_d_seal.sh  (after 4b.3 also done)"
echo "================================================================"
