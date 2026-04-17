#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -A lt200308
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -t 2:00:00
#SBATCH -J inv_4a_C01g
#SBATCH -o logs/inv_4a_C01g_%j.out
#SBATCH -e logs/inv_4a_C01g_%j.err

# =============================================================================
# LAUNCH_C01g_boundary.sh
#
# Phase 4a — boundary catalog. Builds the per-chromosome boundary catalog
# that D11 in C01d (pass-2) consults for concordance scoring. Reads Engine B
# for Cheat 4 (ancestry transition evidence at boundaries) via load_bridge.R.
#
# Called by: orchestrator/run_phase4.sh  (runs in parallel with C01d pass-1)
#
# Upstream:    STEP_C01a precomp RDS
#              PHASE_01C_block_detect landscape output
#              inv_detect_v9.3 scoring tables (staircase)
#              STEP_C01b_1 cores output
#              SV prior (from STEP_C00)
# Downstream:  C01d pass-2 (4e), consumed via --boundary_dir
#
# CLI accepted from orchestrator:
#   --array=<n>            (consumed by sbatch itself, no-op here)
#   --chroms LG01,LG02,... (passed to script as --chrom if single)
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

# ── Parse orchestrator args ──────────────────────────────────────────────
CHROMS_CSV=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --chroms) CHROMS_CSV="$2"; shift 2 ;;
    --array=*) shift ;;
    --array)   shift 2 ;;
    *) echo "[C01g] unknown arg (ignored): $1" >&2; shift ;;
  esac
done

# If a single chrom was given via --chroms, pass it to the script.
CHROM_OPTS=""
if [[ -n "${CHROMS_CSV}" ]] && [[ "${CHROMS_CSV}" != *,* ]]; then
  CHROM_OPTS="--chrom ${CHROMS_CSV}"
fi

# ── Resolve dirs ─────────────────────────────────────────────────────────
PRECOMP_DIR_RUN="${PRECOMP_DIR_RUN:-${MDS_DIR}/snake_regions_multiscale/precomp}"
LANDSCAPE_DIR_RUN="${LANDSCAPE_DIR_RUN:-${MDS_DIR}/snake_regions_multiscale/landscape}"
STAIRCASE_DIR_RUN="${STAIRCASE_DIR_RUN:-${MDS_DIR}/snake_regions_multiscale/candidate_detection_out}"
CORES_DIR_RUN="${CORES_DIR_RUN:-${MDS_DIR}/snake_regions_multiscale}"
SV_PRIOR_DIR_RUN="${SV_PRIOR_DIR_RUN:-${MDS_DIR}/snake_regions_multiscale/sv_prior}"
OUTDIR_RUN="${OUTDIR_RUN:-${MDS_DIR}/snake_regions_multiscale/boundary_catalog}"

mkdir -p "${OUTDIR_RUN}"

for d in "${PRECOMP_DIR_RUN}" "${LANDSCAPE_DIR_RUN}" "${STAIRCASE_DIR_RUN}"; do
  [[ -d "$d" ]] || { echo "[C01g] ERROR: missing dir — $d" >&2; exit 2; }
done

# ── Optional: BAM manifest + mosdepth + repeatmasker for full D-scoring ──
BAM_OPTS=""
[[ -n "${BAM_MANIFEST:-}" && -f "${BAM_MANIFEST}" ]] && BAM_OPTS="${BAM_OPTS} --bam_manifest ${BAM_MANIFEST}"
[[ -n "${MOSDEPTH_DIR:-}" && -d "${MOSDEPTH_DIR}" ]] && BAM_OPTS="${BAM_OPTS} --mosdepth_dir ${MOSDEPTH_DIR}"
[[ -n "${REPEATMASKER_FILE:-}" && -f "${REPEATMASKER_FILE}" ]] && \
  BAM_OPTS="${BAM_OPTS} --repeatmasker ${REPEATMASKER_FILE}"

echo "================================================================"
echo "  4a / C01g — boundary catalog"
echo "  Precomp dir:   ${PRECOMP_DIR_RUN}"
echo "  Landscape dir: ${LANDSCAPE_DIR_RUN}"
echo "  Staircase dir: ${STAIRCASE_DIR_RUN}"
echo "  Cores dir:     ${CORES_DIR_RUN}"
echo "  SV prior dir:  ${SV_PRIOR_DIR_RUN}"
echo "  Out dir:       ${OUTDIR_RUN}"
echo "  Chroms:        ${CHROMS_CSV:-all}"
echo "  SAMPLE_GROUP:  ${SAMPLE_GROUP:-all_226}"
echo "  Bonus opts:    ${BAM_OPTS:-(none)}"
echo "  JobID:         ${SLURM_JOB_ID}"
echo "  Started:       $(date)"
echo "================================================================"

"${RSCRIPT_BIN}" \
  "${BASE}/inversion_modules/phase_4_postprocessing/4a_existence_layers/STEP_C01g_boundary_catalog_wired_4_8_10_11_17_21_v934_registry.R" \
  --precomp     "${PRECOMP_DIR_RUN}" \
  --landscape   "${LANDSCAPE_DIR_RUN}" \
  --staircase   "${STAIRCASE_DIR_RUN}" \
  --cores       "${CORES_DIR_RUN}" \
  --sv_prior    "${SV_PRIOR_DIR_RUN}" \
  --outdir      "${OUTDIR_RUN}" \
  ${CHROM_OPTS} ${BAM_OPTS}

echo "================================================================"
echo "  C01g complete — $(date)"
echo "  Boundary catalog: ${OUTDIR_RUN}/"
echo "  Consumed by: C01d pass-2 (4e) via --boundary_dir"
echo "================================================================"
