#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -A lt200308
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -t 2:00:00
#SBATCH -J inv_4e_C01d_p2
#SBATCH -o logs/inv_4e_C01d_p2_%j.out
#SBATCH -e logs/inv_4e_C01d_p2_%j.err

# =============================================================================
# LAUNCH_C01d_scoring_pass2.sh
#
# Phase 4e — candidate scoring PASS 2. Re-runs C01d scoring with D11
# (boundary concordance) and hypothesis-test outputs now available, so the
# final tier/score table reflects all evidence.
#
# Called by: orchestrator/run_phase4.sh  (depends on 4a boundary + 4c hyp
# + 4d group_cheats).
#
# Upstream:    inv_detect_v9.3 scoring tables (DETECTOR_DIR)
#              STEP_C01a precomp RDS
#              4a C01g boundary catalog (BOUNDARY_DIR)
#              4c C01f hypothesis tests (HYP_DIR)
#              STEP_C01b_1 cores
# Downstream:  LAUNCH_characterize_classify.sh
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

# ── Parse orchestrator-style args ────────────────────────────────────────
CHROMS_CSV=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --chroms) CHROMS_CSV="$2"; shift 2 ;;
    --array=*) shift ;;
    --array)   shift 2 ;;
    *) echo "[C01d_p2] unknown arg (ignored): $1" >&2; shift ;;
  esac
done

# ── Resolve dirs ─────────────────────────────────────────────────────────
DETECTOR_DIR="${DETECTOR_DIR:-${MDS_DIR}/snake_regions_multiscale/candidate_detection_out}"
OUTDIR_RUN="${OUTDIR_RUN:-${SCORING_DIR}/pass2}"
PRECOMP_DIR_RUN="${PRECOMP_DIR_RUN:-${MDS_DIR}/snake_regions_multiscale/precomp}"
CORES_DIR_RUN="${CORES_DIR_RUN:-${MDS_DIR}/snake_regions_multiscale}"
BOUNDARY_DIR_RUN="${BOUNDARY_DIR_RUN:-${MDS_DIR}/snake_regions_multiscale/boundary_catalog}"
HYP_DIR_RUN="${HYP_DIR_RUN:-${MDS_DIR}/snake_regions_multiscale/hypothesis_tests}"

mkdir -p "${OUTDIR_RUN}"

for d in "${DETECTOR_DIR}" "${BOUNDARY_DIR_RUN}" "${HYP_DIR_RUN}"; do
  [[ -d "$d" ]] || { echo "[C01d_p2] ERROR: missing dir — $d" >&2; exit 2; }
done

echo "================================================================"
echo "  4e / C01d pass-2 — candidate scoring (final pass)"
echo "  Detector dir: ${DETECTOR_DIR}"
echo "  Boundary dir: ${BOUNDARY_DIR_RUN}"
echo "  Hyp dir:      ${HYP_DIR_RUN}"
echo "  Precomp dir:  ${PRECOMP_DIR_RUN}"
echo "  Cores dir:    ${CORES_DIR_RUN}"
echo "  Out dir:      ${OUTDIR_RUN}"
echo "  Chroms:       ${CHROMS_CSV:-all}"
echo "  SAMPLE_GROUP: ${SAMPLE_GROUP:-all_226}"
echo "  JobID:        ${SLURM_JOB_ID}"
echo "  Started:      $(date)"
echo "================================================================"

# Pass 2 — WITH --boundary_dir and --hyp_dir (the pass-1/pass-2 distinction)
"${RSCRIPT_BIN}" \
  "${BASE}/inversion_modules/phase_4_catalog/STEP_C01d_candidate_scoring_wired_25_v934_registry.R" \
  "${DETECTOR_DIR}" "${OUTDIR_RUN}" \
  --precomp_dir  "${PRECOMP_DIR_RUN}" \
  --cores_dir    "${CORES_DIR_RUN}" \
  --boundary_dir "${BOUNDARY_DIR_RUN}" \
  --hyp_dir      "${HYP_DIR_RUN}"

echo "================================================================"
echo "  C01d pass-2 complete — $(date)"
echo "  Final candidate table: ${OUTDIR_RUN}/candidate_scores.tsv.gz"
echo "  Next: sbatch LAUNCH_characterize_classify.sh"
echo "================================================================"
