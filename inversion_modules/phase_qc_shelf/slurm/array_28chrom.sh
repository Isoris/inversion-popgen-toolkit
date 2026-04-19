#!/bin/bash
#SBATCH --job-name=QC_shelf
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --array=1-28%8
#SBATCH --output=logs/qc_shelf_%A_%a.out
#SBATCH --error=logs/qc_shelf_%A_%a.err
# =============================================================================
# SLURM driver for MODULE_QC_ShelfDiagnosis
# Runs one chrom per array task. LG01..LG28.
#
# Submit:
#   cd MODULE_QC_ShelfDiagnosis
#   sbatch slurm/array_28chrom.sh
#
# Shelf region can be passed via env (only used by Q04):
#   sbatch --export=ALL,SHELF_START_MB=15,SHELF_END_MB=18 slurm/array_28chrom.sh
# =============================================================================
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${here}/00_config.sh"

# Pick chromosome from chr.list by array index
CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${BEAGLE_DIR}/chr.list")
[[ -z "${CHR}" ]] && qc_die "No chrom at line ${SLURM_ARRAY_TASK_ID}"

qc_log "SLURM task ${SLURM_ARRAY_TASK_ID} -> ${CHR}"

# Activate R environment if needed (adjust for your HPC)
# source ~/.bashrc
# conda activate assembly

bash "${here}/run_all.sh" "${CHR}"
