#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -A lt200308
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=237G
#SBATCH -t 1:00:00
#SBATCH -J inv_01C_block
#SBATCH -o logs/inv_01C_block_%j.out
#SBATCH -e logs/inv_01C_block_%j.err

# =============================================================================
# LAUNCH_PHASE_01C_block_detect.slurm
#
# Phase 2 / 2c — block detection + boundary classification + blue-cross
# diagnosis + block concordance. Run AFTER C01a precompute, BEFORE C01b merge.
#
# Upstream:    STEP_C01a precompute (PRECOMP_DIR)
# Downstream:  STEP_C01b_1_seeded_regions (reads landscape/)
#
# Usage:
#   sbatch LAUNCH_PHASE_01C_block_detect.slurm                 # all chromosomes
#   sbatch LAUNCH_PHASE_01C_block_detect.slurm C_gar_LG28      # single chrom
#
# Sizing: parallel across chroms via mclapply; 28 CPU / 237G matches the
# archive template. Walltime ~30–60 min full-genome.
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

# ── Resolve in/out dirs (consistent with C01a output layout) ─────────────
PRECOMP_DIR_RUN="${PRECOMP_DIR_RUN:-${MDS_DIR}/snake_regions_multiscale/precomp}"
LANDSCAPE_DIR_RUN="${LANDSCAPE_DIR_RUN:-${MDS_DIR}/snake_regions_multiscale/landscape}"

if [[ ! -d "${PRECOMP_DIR_RUN}" ]]; then
  echo "[01C] ERROR: precomp dir missing — ${PRECOMP_DIR_RUN}" >&2
  echo "[01C]        Run LAUNCH_STEP_C01a_precompute.slurm first." >&2
  exit 2
fi

# ── Optional single-chrom arg ────────────────────────────────────────────
CHR_ARG="${1:-}"
CHR_OPTS=""
if [[ -n "${CHR_ARG}" ]]; then
  CHR_OPTS="--chrom ${CHR_ARG}"
fi

# ── Optional gap/AGP annotations (blue-cross) ────────────────────────────
GAP_OPTS=""
AGP_GAPS="${AGP_GAPS:-${BASE}/00-samples/fClaHyb_Gar_LG.gaps.bed}"
AGP_SCAFFOLD="${AGP_SCAFFOLD:-${BASE}/00-samples/fClaHyb_Gar_LG.agp}"
[[ -f "${AGP_GAPS}" ]]     && GAP_OPTS="${GAP_OPTS} --gaps ${AGP_GAPS}"
[[ -f "${AGP_SCAFFOLD}" ]] && GAP_OPTS="${GAP_OPTS} --agp  ${AGP_SCAFFOLD}"

echo "================================================================"
echo "  01C: Block detection + boundary classification"
echo "  Chrom:     ${CHR_ARG:-all}"
echo "  Precomp:   ${PRECOMP_DIR_RUN}"
echo "  Out:       ${LANDSCAPE_DIR_RUN}"
echo "  Mode:      ${POP_MODE:-hatchery}"
echo "  AGP opts:  ${GAP_OPTS:-(none)}"
echo "  JobID:     ${SLURM_JOB_ID}"
echo "  Started:   $(date)"
echo "================================================================"

"${RSCRIPT_BIN}" \
  "${BASE}/inversion_modules/phase_2_discovery/2c_precomp/PHASE_01C_block_detect.R" \
  "${PRECOMP_DIR_RUN}" "${LANDSCAPE_DIR_RUN}" \
  --mode "${POP_MODE:-hatchery}" \
  ${CHR_OPTS} ${GAP_OPTS}

echo "================================================================"
echo "  01C complete — $(date)"
echo "  Landscape: ${LANDSCAPE_DIR_RUN}/"
echo "  Next:      sbatch LAUNCH_STEP_C01b_1_seeded_regions.slurm"
echo "================================================================"
