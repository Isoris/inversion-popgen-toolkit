#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -A lt200308
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -t 1:00:00
#SBATCH -J inv_4b4_seal
#SBATCH -o logs/inv_4b4_seal_%j.out
#SBATCH -e logs/inv_4b4_seal_%j.err

# =============================================================================
# LAUNCH_C01i_d_seal.sh
#
# Phase 4b.4 — seal. Reads the three upstream Tier-2 blocks, resolves per-
# sample FINAL class per the PHASE4B_REWRITE_ARCHITECTURE §3 rules, and
# registers HOM_REF / HET / HOM_INV / RECOMBINANT groups (plus optional
# GC / DCO subgroups, HOM_STD alias) in sample_registry.
#
# Sets q6_group_validation = UNCERTAIN with a quality_flag. C01f may
# promote above UNCERTAIN except when composite_flag = likely_composite.
#
# Called by: orchestrator/run_phase4b.sh  (depends on 4b.1, 4b.2, 4b.3).
#
# Upstream:    4b.1 decompose (internal_dynamics block)
#              4b.2 multi_recomb (recombinant_map block)
#              4b.3 nested_comp (internal_ancestry_composition block)
# Downstream:  4c hypothesis tests (C01f reads sealed groups from registry)
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
[[ -f "${BRIDGE_SH}" ]] || { echo "Missing bridge (need SAMPLE_REGISTRY FK)" >&2; exit 1; }
source "${BRIDGE_SH}"

inv_init_dirs
mkdir -p logs

# ── Args ─────────────────────────────────────────────────────────────────
CANDIDATES="${CANDIDATES:-${SCORING_DIR}/pass1/candidate_scores.tsv.gz}"
DECOMP_OUT="${DECOMP_OUT:-${MDS_DIR}/snake_regions_multiscale/phase4b_out/decompose}"
RECOMB_OUT="${RECOMB_OUT:-${MDS_DIR}/snake_regions_multiscale/phase4b_out/multi_recomb}"
NESTED_OUT="${NESTED_OUT:-${MDS_DIR}/snake_regions_multiscale/phase4b_out/nested_comp}"
SEAL_OUT="${SEAL_OUT:-${MDS_DIR}/snake_regions_multiscale/phase4b_out/seal}"
TIER_MAX="${TIER_MAX:-3}"

mkdir -p "${SEAL_OUT}"

[[ -f "${CANDIDATES}" ]] || { echo "[C01i_d] ERROR: missing ${CANDIDATES}" >&2; exit 2; }
for d in "${DECOMP_OUT}" "${RECOMB_OUT}"; do
  [[ -d "$d" ]] || { echo "[C01i_d] ERROR: upstream dir missing — $d" >&2; exit 2; }
done
if [[ ! -d "${NESTED_OUT}" ]]; then
  echo "[C01i_d] WARN: nested_comp dir missing — ${NESTED_OUT}" >&2
  echo "[C01i_d]        seal will write unknown_no_engine_b stubs." >&2
fi

# ── SAMPLE_GROUP sanity — seal writes into sample_registry with FK ──────
SAMPLE_GROUP="${SAMPLE_GROUP:-all_226}"
if [[ -n "${SAMPLE_REGISTRY:-}" && -f "${SAMPLE_REGISTRY}/groups/${SAMPLE_GROUP}.txt" ]]; then
  echo "[C01i_d] Sample group: ${SAMPLE_GROUP} (registered)"
else
  echo "[C01i_d] WARN: SAMPLE_GROUP=${SAMPLE_GROUP} not registered" >&2
fi
export SAMPLE_GROUP

echo "================================================================"
echo "  4b.4 / C01i seal"
echo "  Candidates:   ${CANDIDATES}"
echo "  Decomp dir:   ${DECOMP_OUT}"
echo "  Recomb dir:   ${RECOMB_OUT}"
echo "  Nested dir:   ${NESTED_OUT}"
echo "  Out dir:      ${SEAL_OUT}"
echo "  Tier max:     ${TIER_MAX}"
echo "  SAMPLE_GROUP: ${SAMPLE_GROUP}"
echo "  JobID:        ${SLURM_JOB_ID}"
echo "  Started:      $(date)"
echo "================================================================"

"${RSCRIPT_BIN}" \
  "${BASE}/inversion_modules/phase_4_postprocessing/4d_group_proposal/STEP_C01i_d_seal.R" \
  --candidates "${CANDIDATES}" \
  --decomp_dir "${DECOMP_OUT}" \
  --recomb_dir "${RECOMB_OUT}" \
  --nested_dir "${NESTED_OUT}" \
  --outdir     "${SEAL_OUT}" \
  --tier_max   "${TIER_MAX}"

echo "================================================================"
echo "  4b.4 complete — $(date)"
echo "  Seal out: ${SEAL_OUT}/"
echo "  Groups registered in sample_registry with q6_group_validation=UNCERTAIN"
echo "  Next (4c): sbatch LAUNCH_C01f_hypothesis.sh"
echo "================================================================"
