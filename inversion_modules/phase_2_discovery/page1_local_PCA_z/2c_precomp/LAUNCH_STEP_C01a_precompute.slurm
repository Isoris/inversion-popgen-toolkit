#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -A lt200308
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=237G
#SBATCH -t 2:00:00
#SBATCH -J inv_C01a
#SBATCH -o logs/inv_C01a_%j.out
#SBATCH -e logs/inv_C01a_%j.err

# =============================================================================
# LAUNCH_STEP_C01a_precompute.slurm
#
# Phase 2 / 2c — per-chromosome precompute. Builds similarity matrices,
# NN distances, robust z-scores, inv-likeness, AND stamps localQ_* columns
# from the ancestry precompute + SV overlap columns from STEP_C00 output.
#
# Run AFTER:
#   1. LAUNCH_instant_q_precompute.slurm   (at least K=8; ancestry precompute)
#   2. LAUNCH_STEP_C00_sv_prior.slurm      (SV evidence prior)
#
# Chat-16 note: merge_local_Q_into_invlikeness() in this script was the
# silent-no-op target. SAMPLE_GROUP=all_226 (env default) must match the
# filename convention the ancestry launcher wrote with, e.g.
#   ${LOCAL_Q_DIR}/K08/C_gar_LG12.all_226.local_Q_summary.tsv.gz
# If that file doesn't exist, the merge is a no-op and localQ_* columns
# come out as NA.
#
# Usage:
#   sbatch LAUNCH_STEP_C01a_precompute.slurm
#
# Sizing: 28 cores for parallel per-chromosome precompute via mclapply
# (one chrom per core; 28 chroms maps 1:1). 237G is conservative head-
# room for the full-genome precompute. Walltime ~30–60 min.
# =============================================================================

set -euo pipefail
source ~/.bashrc
mamba activate assembly

# ── Load configs ─────────────────────────────────────────────────────────
CONFIG="${CONFIG:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_modules/00_inversion_config.sh}"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
set -a
source "${CONFIG}"
set +a

# Pipeline bridge gives SAMPLE_GROUP, RESULTS_REGISTRY_DIR, SAMPLE_REGISTRY,
# LOCAL_Q_DIR — all needed for the chat-16 localQ merge to land correctly.
BRIDGE_SH="${BRIDGE_SH:-${BASE}/utils/pipeline_bridge.sh}"
[[ -f "${BRIDGE_SH}" ]] || { echo "Missing bridge: ${BRIDGE_SH}" >&2; exit 1; }
source "${BRIDGE_SH}"

inv_init_dirs

mkdir -p logs

# ── Resolve upstream paths ───────────────────────────────────────────────
# STEP_C00 writes to ${MDS_DIR}/snake_regions_multiscale/sv_prior (hard-
# coded in the script). If you moved it, adjust here or pass SV_PRIOR_DIR_RUN.
SV_PRIOR_DIR_RUN="${SV_PRIOR_DIR_RUN:-${MDS_DIR}/snake_regions_multiscale/sv_prior}"
if [[ ! -d "${SV_PRIOR_DIR_RUN}" ]]; then
  echo "[C01a] WARN: SV prior dir missing — ${SV_PRIOR_DIR_RUN}" >&2
  echo "[C01a]       SV annotation columns will be 0/NA. Run STEP_C00 first." >&2
fi

# STEP_C01a positional args: <step10_outprefix> <outdir>
# step10 prefix = where STEP10 (window MDS) wrote its per-chrom outputs
STEP10_PREFIX="${STEP10_PREFIX:-${MDS_DIR}/inversion_localpca}"
OUTDIR="${PRECOMP_OUTDIR:-${MDS_DIR}/snake_regions_multiscale}"

# ── SAMPLE_GROUP sanity — matches ancestry launcher's filename convention ─
SAMPLE_GROUP="${SAMPLE_GROUP:-all_226}"
if [[ -n "${SAMPLE_REGISTRY:-}" && -f "${SAMPLE_REGISTRY}/groups/${SAMPLE_GROUP}.txt" ]]; then
  echo "[C01a] Sample group: ${SAMPLE_GROUP} (registered)"
else
  echo "[C01a] WARN: SAMPLE_GROUP=${SAMPLE_GROUP} not registered in sample_registry" >&2
  echo "[C01a]       FK-enforced put_* calls downstream will fail." >&2
fi
export SAMPLE_GROUP

# ── localQ cache existence check (the chat-16 silent-bug surface) ────────
# If localQ files for SAMPLE_GROUP don't exist, merge_local_Q_into_invlikeness()
# no-ops and localQ_* columns come out NA. Warn early.
CANONICAL_KK=$(printf "K%02d" "${CANONICAL_K:-8}")
LOCALQ_PROBE="${LOCAL_Q_DIR}/${CANONICAL_KK}"
if [[ -d "${LOCALQ_PROBE}" ]]; then
  N_LOCALQ=$(find "${LOCALQ_PROBE}" -name "*.${SAMPLE_GROUP}.local_Q_summary.tsv.gz" 2>/dev/null | wc -l)
  echo "[C01a] localQ cache (${CANONICAL_KK}, ${SAMPLE_GROUP}): ${N_LOCALQ} chrom files"
  if [[ "${N_LOCALQ}" -eq 0 ]]; then
    echo "[C01a] WARN: no localQ files found — localQ_* columns will be NA" >&2
    echo "[C01a]       Run LAUNCH_instant_q_precompute.slurm first." >&2
  fi
else
  echo "[C01a] WARN: LOCAL_Q_DIR/${CANONICAL_KK} missing — ${LOCALQ_PROBE}" >&2
fi

echo "================================================================"
echo "  C01a: per-chrom precompute (parallel, ${SLURM_CPUS_PER_TASK} cores)"
echo "  STEP10 prefix:    ${STEP10_PREFIX}"
echo "  Out dir:          ${OUTDIR}"
echo "  SV prior dir:     ${SV_PRIOR_DIR_RUN}"
echo "  Dosage dir:       ${DOSAGE_DIR}"
echo "  SAMPLE_GROUP:     ${SAMPLE_GROUP}"
echo "  LOCAL_Q_DIR:      ${LOCAL_Q_DIR}"
echo "  RESULTS_REGISTRY: ${RESULTS_REGISTRY_DIR}"
echo "  JobID:            ${SLURM_JOB_ID}"
echo "  Started:          $(date)"
echo "================================================================"

"${RSCRIPT_BIN}" \
  "${BASE}/inversion_modules/phase_2_discovery/2c_precomp/STEP_C01a_precompute.R" \
  "${STEP10_PREFIX}" \
  "${OUTDIR}" \
  --dosage_dir   "${DOSAGE_DIR}" \
  --sv_prior_dir "${SV_PRIOR_DIR_RUN}" \
  --mode         "${POP_MODE:-hatchery}"

echo "================================================================"
echo "  C01a complete — $(date)"
echo "  Precomp RDS:  ${OUTDIR}/precomp/"
echo ""
echo "  Sanity check — confirm localQ columns got stamped:"
echo "    Rscript -e 'x <- readRDS(\"${OUTDIR}/precomp/C_gar_LG12.precomp.rds\");"
echo "                print(grep(\"^localQ_\", names(x), value = TRUE));"
echo "                cat(\"non-NA:\", sum(!is.na(x\$localQ_delta12_K08)), \"/\", nrow(x), \"\\n\")'"
echo "  Expected: 3 localQ_* columns, non-NA count == nrow."
echo ""
echo "  Next: sbatch LAUNCH_01C_block_detect.slurm"
echo "================================================================"
