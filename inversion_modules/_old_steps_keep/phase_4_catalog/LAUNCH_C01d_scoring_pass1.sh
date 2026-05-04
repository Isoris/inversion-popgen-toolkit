#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -A lt200308
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -t 2:00:00
#SBATCH -J inv_4a_C01d_p1
#SBATCH -o logs/inv_4a_C01d_p1_%j.out
#SBATCH -e logs/inv_4a_C01d_p1_%j.err

# =============================================================================
# LAUNCH_C01d_scoring_pass1.sh
#
# Phase 4a — candidate scoring PASS 1 (cold pass, no boundary/hypothesis).
# Consumes the v9.3 inv_detect scoring tables and emits the Tier-1/2/3
# candidate list. This pass runs BEFORE 4c hypothesis tests so it does NOT
# take --boundary_dir or --hyp_dir. Pass 2 (4e) re-scores with both.
#
# Called by: orchestrator/run_phase4.sh
#
# Upstream:    inv_detect_v9.3 scoring tables (DETECTOR_DIR)
#              STEP_C01a precomp RDS
#              STEP_C01b_1 cores output
# Downstream:  Phase 4b (decomposition reads candidate list)
#              LAUNCH_C01d_scoring_pass2.sh (4e, pass-2)
#
# CLI accepted from orchestrator:
#   --array=<n>            (consumed by sbatch itself, no-op here)
#   --chroms LG01,LG02,... (filters to specific chromosomes via CHROMS_CSV)
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

# ── Parse orchestrator args (--chroms only used for provenance echo) ─────
CHROMS_CSV=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --chroms) CHROMS_CSV="$2"; shift 2 ;;
    --array=*) shift ;;                  # no-op (sbatch consumed it)
    --array)   shift 2 ;;                # no-op (sbatch consumed it)
    *) echo "[C01d_p1] unknown arg (ignored): $1" >&2; shift ;;
  esac
done

# ── Resolve dirs ─────────────────────────────────────────────────────────
DETECTOR_DIR="${DETECTOR_DIR:-${MDS_DIR}/snake_regions_multiscale/candidate_detection_out}"
OUTDIR_RUN="${OUTDIR_RUN:-${SCORING_DIR}/pass1}"
PRECOMP_DIR_RUN="${PRECOMP_DIR_RUN:-${MDS_DIR}/snake_regions_multiscale/precomp}"
CORES_DIR_RUN="${CORES_DIR_RUN:-${MDS_DIR}/snake_regions_multiscale}"

mkdir -p "${OUTDIR_RUN}"

[[ -d "${DETECTOR_DIR}" ]] || {
  echo "[C01d_p1] ERROR: detector dir missing — ${DETECTOR_DIR}" >&2
  echo "[C01d_p1]          Run phase_2d LAUNCH_run_all.slurm first." >&2
  exit 2;
}

echo "================================================================"
echo "  4a / C01d pass-1 — candidate scoring (cold pass)"
echo "  Detector dir: ${DETECTOR_DIR}"
echo "  Out dir:      ${OUTDIR_RUN}"
echo "  Precomp dir:  ${PRECOMP_DIR_RUN}"
echo "  Cores dir:    ${CORES_DIR_RUN}"
echo "  Chroms:       ${CHROMS_CSV:-all}"
echo "  SAMPLE_GROUP: ${SAMPLE_GROUP:-all_226}"
echo "  JobID:        ${SLURM_JOB_ID}"
echo "  Started:      $(date)"
echo "================================================================"

# Pass 1 — no --boundary_dir / --hyp_dir (those come from 4a-boundary and 4c)
"${RSCRIPT_BIN}" \
  "${BASE}/inversion_modules/phase_4_catalog/STEP_C01d_candidate_scoring_wired_25_v934_registry.R" \
  "${DETECTOR_DIR}" "${OUTDIR_RUN}" \
  --precomp_dir "${PRECOMP_DIR_RUN}" \
  --cores_dir   "${CORES_DIR_RUN}"

echo "================================================================"
echo "  C01d pass-1 complete — $(date)"
echo "  Candidate list: ${OUTDIR_RUN}/candidate_scores.tsv.gz"
echo "  Next (4b):      run_phase4b.sh  (decomposition + recombinant + seal)"
echo "================================================================"
