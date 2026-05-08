#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -A lt200308
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=237G
#SBATCH -t 1:00:00
#SBATCH -J inv_C01b1
#SBATCH -o logs/inv_C01b1_%j.out
#SBATCH -e logs/inv_C01b1_%j.err

# =============================================================================
# LAUNCH_STEP_C01b_1_seeded_regions.slurm
#
# Phase 2 / 2c — seeded region detection (cores pass). Takes C01a precomp and
# 01C landscape/block registry; emits seeded_regions_windows_<chr>.tsv.gz.
#
# Upstream:    STEP_C01a_precompute.R  (precomp RDS)
#              PHASE_01C_block_detect.R (landscape/ block_registry)
# Downstream:  STEP_C01d candidate scoring (phase 4a)
#
# Usage:
#   sbatch LAUNCH_STEP_C01b_1_seeded_regions.slurm                   # all chroms
#   sbatch LAUNCH_STEP_C01b_1_seeded_regions.slurm C_gar_LG12        # single chrom
#
# Optional env: TAG="my_label"  to tag this run's output
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

# ── Resolve dirs ─────────────────────────────────────────────────────────
PRECOMP_DIR_RUN="${PRECOMP_DIR_RUN:-${MDS_DIR}/snake_regions_multiscale/precomp}"
OUTDIR_RUN="${OUTDIR_RUN:-${MDS_DIR}/snake_regions_multiscale}"

[[ -d "${PRECOMP_DIR_RUN}" ]] || {
  echo "[C01b1] ERROR: precomp dir missing — ${PRECOMP_DIR_RUN}" >&2; exit 2;
}

# ── Optional single-chrom arg / tag ──────────────────────────────────────
CHR_ARG="${1:-}"
CHR_OPTS=""
if [[ -n "${CHR_ARG}" ]]; then
  CHR_OPTS="--chr ${CHR_ARG}"
fi

TAG_OPTS=""
if [[ -n "${TAG:-}" ]]; then
  TAG_OPTS="--tag ${TAG}"
fi

echo "================================================================"
echo "  C01b_1: Seeded regions (cores pass)"
echo "  Chrom:   ${CHR_ARG:-all}"
echo "  Precomp: ${PRECOMP_DIR_RUN}"
echo "  Out:     ${OUTDIR_RUN}"
echo "  Tag:     ${TAG:-(none)}"
echo "  JobID:   ${SLURM_JOB_ID}"
echo "  Started: $(date)"
echo "================================================================"

"${RSCRIPT_BIN}" \
  "${BASE}/inversion_modules/phase_2_discovery/2c_precomp/STEP_C01b_1_seeded_regions.R" \
  "${PRECOMP_DIR_RUN}" "${OUTDIR_RUN}" \
  ${CHR_OPTS} ${TAG_OPTS}

echo "================================================================"
echo "  C01b_1 complete — $(date)"
echo "  Seeded regions: ${OUTDIR_RUN}/seeded_regions_windows_<chr>.tsv.gz"
echo "  Next:           sbatch phase_2d/LAUNCH_run_all.slurm"
echo "================================================================"
