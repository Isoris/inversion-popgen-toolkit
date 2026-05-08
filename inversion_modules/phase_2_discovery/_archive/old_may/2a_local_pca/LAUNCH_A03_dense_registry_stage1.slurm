#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48GB
#SBATCH -t 6:00:00
#SBATCH -J inv_09b_s1
#SBATCH -o logs/inv_09b_s1.%A_%a.out
#SBATCH -e logs/inv_09b_s1.%A_%a.err

# =============================================================================
# LAUNCH_A03_dense_registry_stage1.slurm  (v8.5)
#
# Stage 1: one chromosome per array task
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
WINSIZE="${2:-100}"
STEP_SIZE="${3:-10}"
NPC_ARG="${4:-${NPC}}"

A2_STEPS="${DISCOVERY2DIR}/steps"
REGISTRY_DIR="${INVDIR}/05b_dense_registry"
mkdir -p "${REGISTRY_DIR}/tmp" "${INVDIR}/logs"

IDX="${SLURM_ARRAY_TASK_ID}"
LINE=$((IDX + 1))

CHROM="$(sed -n "${LINE}p" "${CHROM_LIST}")"
if [[ -z "${CHROM}" ]]; then
  echo "[ERROR] Could not read chromosome from ${CHROM_LIST} line ${LINE}" >&2
  exit 1
fi

echo "================================================================"
echo "  STEP09b Stage 1 — ${CHROM}"
echo "  Array task: ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_TASK_MAX:-?}"
echo "  JobID: ${SLURM_JOB_ID}"
echo "  Started: $(date)"
echo "================================================================"

"${RSCRIPT_BIN}" "${A2_STEPS}/STEP_A03_dense_registry_stage1.R" \
  --dosage_dir "${DOSAGE_DIR}" \
  --outdir "${REGISTRY_DIR}" \
  --chrom "${CHROM}" \
  --winsize "${WINSIZE}" \
  --step "${STEP_SIZE}" \
  --npc "${NPC_ARG}"

echo "[DONE] STEP09b Stage 1 — ${CHROM} — $(date)"
