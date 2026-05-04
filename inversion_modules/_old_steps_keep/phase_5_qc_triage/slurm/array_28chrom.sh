#!/bin/bash
#SBATCH --job-name=phase_qc_shelf
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --array=1-28%8
#SBATCH --output=logs/phase_qc_shelf_%A_%a.out
#SBATCH --error=logs/phase_qc_shelf_%A_%a.err
# =============================================================================
# slurm/array_28chrom.sh — run the full Q01→Q10 pipeline for every chromosome
# =============================================================================
# One array task per chromosome, throttled at 8 concurrent.
# CHR list comes from ${BEAGLE_DIR}/chr.list (one chromosome per line).
#
# Per-chromosome runs:
#   Q01 snp density → Q02 BEAGLE uncertainty → Q03 coverage →
#   Q05 theta → Q06 ancestry (precompute + tracks + multiscale) →
#   Q07 popstats + invgt grouping → Q04 diagnostic (both variants) →
#   Q08 shelf heatmap → Q09 gap characterization → Q10 registry ingest
#
# Optional per-chrom shelf coords via shelf_coords.tsv:
#   chrom<TAB>shelf_start_mb<TAB>shelf_end_mb<TAB>bp1_mb<TAB>bp2_mb
#
# Submit:
#   cd /scratch/.../phase_qc_shelf
#   mkdir -p logs
#   sbatch slurm/array_28chrom.sh
#
# After all complete:
#   bash STEP_Q09_gap_characterization.sh ALL
#   bash STEP_Q10_register.sh ALL
#   bash scripts/registry_query.sh summary
# =============================================================================
set -euo pipefail

MOD_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${MOD_ROOT}"
source ./00_config.sh

CHR_LIST="${BEAGLE_DIR}/chr.list"
[[ -f "${CHR_LIST}" ]] || { echo "chr.list not found at ${CHR_LIST}"; exit 1; }

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${CHR_LIST}")
[[ -z "${CHR}" ]] && { echo "No chrom at line ${SLURM_ARRAY_TASK_ID}"; exit 1; }

SHELF_COORDS_TSV="${SHELF_COORDS_TSV:-${MOD_ROOT}/shelf_coords.tsv}"
if [[ -f "${SHELF_COORDS_TSV}" ]]; then
  row=$(awk -v c="${CHR}" '$1==c' "${SHELF_COORDS_TSV}" | head -1)
  if [[ -n "${row}" ]]; then
    SHELF_START_MB=$(echo "${row}" | awk '{print $2}')
    SHELF_END_MB=$(echo   "${row}" | awk '{print $3}')
    BP1_MB=$(echo         "${row}" | awk '{print $4}')
    BP2_MB=$(echo         "${row}" | awk '{print $5}')
    echo "[slurm] ${CHR}: ${SHELF_START_MB}-${SHELF_END_MB} Mb  BP:${BP1_MB}/${BP2_MB}"
  fi
fi

: "${SHELF_START_MB:=0}"
: "${SHELF_END_MB:=1000}"
: "${BP1_MB:=}"
: "${BP2_MB:=}"
: "${SMOOTH_WIN:=5}"
: "${SNP_DENSITY_SCALE_KB:=10}"
: "${SAMPLE_GROUP:=all_226}"
: "${METHOD_TAG:=phase_qc_shelf}"

export SHELF_START_MB SHELF_END_MB BP1_MB BP2_MB
export SMOOTH_WIN SNP_DENSITY_SCALE_KB
export SAMPLE_GROUP METHOD_TAG

echo "[slurm] ${CHR}: starting at $(date -u +%Y-%m-%dT%H:%M:%SZ)"
bash run_chrom.sh "${CHR}"
echo "[slurm] ${CHR}: done at $(date -u +%Y-%m-%dT%H:%M:%SZ)"
