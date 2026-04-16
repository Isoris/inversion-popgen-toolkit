#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64GB
#SBATCH -t 12:00:00
#SBATCH -J inv_10v2_s1
#SBATCH -o logs/inv_10v2_s1.%A_%a.out
#SBATCH -e logs/inv_10v2_s1.%A_%a.err

# =============================================================================
# LAUNCH_STEP10v2_parallel.slurm  (v8.3-parallel)
#
# Stage 1: one focal chromosome per array task
# =============================================================================

set -euo pipefail
source ~/.bashrc
mamba activate assembly

CONFIG="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_codebase_v8.5/00_inversion_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }

set -a
source "${CONFIG}"
set +a

inv_init_dirs

CHROM_LIST="${1:?ERROR: provide chrom.list as argv1}"
MDS_MODE_ARG="${2:-chunked_4x}"
MDS_SEED_ARG="${3:-42}"

A2_STEPS="${DISCOVERY2DIR}/steps"
mkdir -p "${MDS_DIR}/tmp" "${INVDIR}/logs"

IDX="${SLURM_ARRAY_TASK_ID}"
LINE=$((IDX + 1))

CHROM="$(sed -n "${LINE}p" "${CHROM_LIST}")"
if [[ -z "${CHROM}" ]]; then
  echo "[ERROR] Could not read chromosome from ${CHROM_LIST} line ${LINE}" >&2
  exit 1
fi

RDS_DIR="${INVDIR}/05b_dense_registry"
[[ -d "${RDS_DIR}" ]] && [[ -n "$(ls "${RDS_DIR}"/*.rds 2>/dev/null)" ]] || RDS_DIR="${LPCA_DIR}"

echo "================================================================"
echo "  STEP10v2 Stage 1 — ${CHROM} (mode=${MDS_MODE_ARG})"
echo "  Array task: ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_TASK_MAX:-?}"
echo "  RDS source: ${RDS_DIR}"
echo "  JobID: ${SLURM_JOB_ID}"
echo "  Started: $(date)"
echo "================================================================"

"${RSCRIPT_BIN}" "${A2_STEPS}/STEP10v2_stage1_perchr_mds.R" \
  --rds_dir "${RDS_DIR}" \
  --outdir "${MDS_DIR}" \
  --outprefix "inversion_localpca" \
  --focal_chr "${CHROM}" \
  --mds_mode "${MDS_MODE_ARG}" \
  --npc "${NPC}" \
  --mds_dims 20 \
  --z_thresh 3 \
  --seed "${MDS_SEED_ARG}"

echo "[DONE] STEP10v2 Stage 1 — ${CHROM} — $(date)"
