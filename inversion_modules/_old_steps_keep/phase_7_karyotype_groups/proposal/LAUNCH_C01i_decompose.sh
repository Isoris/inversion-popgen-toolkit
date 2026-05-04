#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -A lt200308
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -t 4:00:00
#SBATCH -J inv_4b1_decomp
#SBATCH -o logs/inv_4b1_decomp_%j.out
#SBATCH -e logs/inv_4b1_decomp_%j.err

# =============================================================================
# LAUNCH_C01i_decompose.sh
#
# Phase 4b.1 — per-candidate decomposition. For each Tier ≤ TIER_MAX
# candidate, cluster samples into HOM_REF / HET / HOM_INV using PCA +
# phasing + flashlight seeds. Writes internal_dynamics Tier-2 JSON block.
#
# Called by: orchestrator/run_phase4b.sh (currently via --wrap; this is a
# standalone-launcher alternative).
#
# Upstream:    candidate_scores.tsv.gz (from C01d pass-1)
# Downstream:  4b.2 multi_recomb (reads decomp output)
#              4b.4 seal (reads decomp + recomb + nested)
#
# NOTE ON SCRIPT PATH: run_phase4b.sh expects the script under
# ${MODULE_DIR}/phase4b_rewrite/R/. The current tree (post pass 15) has
# it under phase_7_karyotype_groups/proposal/. This launcher uses the
# current-tree path.
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
TIER_MAX="${TIER_MAX:-3}"
K_DECOMP="${K_DECOMP:-3}"
STEP03_SEEDS_DIR="${STEP03_SEEDS_DIR:-}"

mkdir -p "${DECOMP_OUT}"

[[ -f "${CANDIDATES}" ]] || {
  echo "[C01i_decomp] ERROR: candidate table not found — ${CANDIDATES}" >&2
  echo "[C01i_decomp]            Run LAUNCH_C01d_scoring_pass1.sh first." >&2
  exit 2;
}

SEED_OPTS=""
if [[ -n "${STEP03_SEEDS_DIR}" ]]; then
  [[ -d "${STEP03_SEEDS_DIR}" ]] || {
    echo "[C01i_decomp] ERROR: step03 seeds dir missing — ${STEP03_SEEDS_DIR}" >&2; exit 2;
  }
  SEED_OPTS="--step03_seeds_dir ${STEP03_SEEDS_DIR}"
fi

echo "================================================================"
echo "  4b.1 / C01i decompose"
echo "  Candidates:   ${CANDIDATES}"
echo "  Out dir:      ${DECOMP_OUT}"
echo "  Tier max:     ${TIER_MAX}"
echo "  K decomp:     ${K_DECOMP}"
echo "  Seeds dir:    ${STEP03_SEEDS_DIR:-(none)}"
echo "  SAMPLE_GROUP: ${SAMPLE_GROUP:-all_226}"
echo "  JobID:        ${SLURM_JOB_ID}"
echo "  Started:      $(date)"
echo "================================================================"

"${RSCRIPT_BIN}" \
  "${BASE}/inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_decompose.R" \
  --candidates "${CANDIDATES}" \
  --outdir     "${DECOMP_OUT}" \
  --tier_max   "${TIER_MAX}" \
  --k_decomp   "${K_DECOMP}" \
  ${SEED_OPTS}

echo "================================================================"
echo "  4b.1 complete — $(date)"
echo "  Decompose out: ${DECOMP_OUT}/"
echo "  Next (4b.2):   sbatch LAUNCH_C01i_b_multi_recomb.sh"
echo "================================================================"
