#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16GB
#SBATCH -t 00:30:00
#SBATCH -J inv_09b_s2
#SBATCH -o logs/inv_09b_s2.%j.out
#SBATCH -e logs/inv_09b_s2.%j.err

# =============================================================================
# LAUNCH_A03_dense_registry_stage2.slurm  (v8.5)
#
# Stage 2: merge per-chromosome temp outputs into master registry
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

A2_STEPS="${DISCOVERY2DIR}/steps"
REGISTRY_DIR="${INVDIR}/05b_dense_registry"
mkdir -p "${REGISTRY_DIR}" "${INVDIR}/logs"

echo "================================================================"
echo "  STEP09b Stage 2 — Master registry merge"
echo "  JobID: ${SLURM_JOB_ID}"
echo "  Started: $(date)"
echo "================================================================"

"${RSCRIPT_BIN}" "${A2_STEPS}/STEP_A03_dense_registry_stage2.R" \
  --outdir "${REGISTRY_DIR}"

echo ""
echo "[DONE] STEP09b Stage 2 — Master registry complete — $(date)"
echo "  Registry: ${REGISTRY_DIR}/windows_master.tsv.gz"
echo "  Summary:  ${REGISTRY_DIR}/windows_master_summary.tsv"
