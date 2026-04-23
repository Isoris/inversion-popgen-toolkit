#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -A lt200308
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -t 2:00:00
#SBATCH -J inv_4b3_nested
#SBATCH -o logs/inv_4b3_nested_%j.out
#SBATCH -e logs/inv_4b3_nested_%j.err

# =============================================================================
# LAUNCH_C01i_c_nested_comp.sh
#
# Phase 4b.3 — nested ancestry composition. Reads per-candidate Engine B
# Q-matrix cache and checks for sub-candidate ancestry structure inside
# each candidate interval. Parallel with 4b.1 decompose (no dependency).
#
# Called by: orchestrator/run_phase4b.sh.
#
# ⚠ STUB WARNING ⚠
# The Python script this launcher calls —
# STEP_C01i_c_nested_composition.py — is referenced by run_phase4b.sh but
# is NOT present in the current inversion_modules tree. This launcher is
# written in the correct shape so that when the script is added, it will
# Just Work. Until then, it will fail loudly on the script-exists check
# below. That is the correct behaviour: a missing phase-4b.3 should fail
# the DAG, not silently skip.
#
# Upstream:    candidate_scores.tsv.gz (from C01d pass-1)
#              ${LOCAL_Q_DIR}/K*/  (Engine B q-matrix cache from instant_q)
# Downstream:  C01i_d_seal (reads nested output for final class resolution)
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
Q_CACHE_DIR="${Q_CACHE_DIR:-${LOCAL_Q_DIR}}"
NESTED_OUT="${NESTED_OUT:-${MDS_DIR}/snake_regions_multiscale/phase4b_out/nested_comp}"
TIER_MAX="${TIER_MAX:-3}"

mkdir -p "${NESTED_OUT}"

# ── STUB GUARD — fail loudly if the script is missing ────────────────────
NESTED_SCRIPT="${BASE}/inversion_modules/phase_4_postprocessing/4d_group_proposal/STEP_C01i_c_nested_composition.py"
if [[ ! -f "${NESTED_SCRIPT}" ]]; then
  cat >&2 <<EOF
================================================================
  4b.3 / C01i nested_composition — SCRIPT MISSING
================================================================
  Expected at: ${NESTED_SCRIPT}

  This phase is referenced by run_phase4b.sh but the Python script
  has not been added to the tree yet. Options:
    (a) Add the script and re-run this launcher.
    (b) Skip 4b.3 for this run: the seal step (4b.4) will write
        unknown_no_engine_b stubs for each candidate.
        Set SKIP_4B3=1 in the orchestrator to signal this explicitly.

  Failing the job so the DAG blocks instead of silently skipping.
================================================================
EOF
  exit 2
fi

[[ -f "${CANDIDATES}" ]] || { echo "[C01i_c] ERROR: missing ${CANDIDATES}" >&2; exit 2; }
[[ -d "${Q_CACHE_DIR}" ]] || {
  echo "[C01i_c] WARN: Q cache dir missing — ${Q_CACHE_DIR}" >&2
  echo "[C01i_c]        nested_comp will write unknown_no_engine_b stubs." >&2
}

echo "================================================================"
echo "  4b.3 / C01i nested_composition"
echo "  Candidates:   ${CANDIDATES}"
echo "  Q cache:      ${Q_CACHE_DIR}"
echo "  Out dir:      ${NESTED_OUT}"
echo "  Tier max:     ${TIER_MAX}"
echo "  SAMPLE_GROUP: ${SAMPLE_GROUP:-all_226}"
echo "  JobID:        ${SLURM_JOB_ID}"
echo "  Started:      $(date)"
echo "================================================================"

python3 "${NESTED_SCRIPT}" \
  --candidates  "${CANDIDATES}" \
  --q_cache_dir "${Q_CACHE_DIR}" \
  --outdir      "${NESTED_OUT}" \
  --tier_max    "${TIER_MAX}"

echo "================================================================"
echo "  4b.3 complete — $(date)"
echo "  Nested out: ${NESTED_OUT}/"
echo "================================================================"
