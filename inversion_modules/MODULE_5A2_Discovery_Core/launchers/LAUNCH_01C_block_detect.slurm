#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=237G
#SBATCH -t 1:00:00
#SBATCH -J inv_01C
#SBATCH -o logs/inv_01C_%j.out
#SBATCH -e logs/inv_01C_%j.err

# =============================================================================
# LAUNCH_01C_block_detect.slurm
#
# PHASE_01C: Block detection + boundary classification + blue-cross diagnosis
# + block concordance. Run AFTER C01a precompute, BEFORE C01b merge.
#
# Outputs consumed by C01b_2 merge (--landscape_dir) and downstream.
#
# Usage:
#   sbatch LAUNCH_01C_block_detect.slurm                 # all chromosomes
#   sbatch LAUNCH_01C_block_detect.slurm C_gar_LG28      # single chromosome
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

A2_SNAKES="${DISCOVERY2DIR}/snakes"
PRECOMP_DIR="${MDS_DIR}/snake_regions_multiscale/precomp"
LANDSCAPE_DIR="${MDS_DIR}/snake_regions_multiscale/landscape"

CHROM_ARG="${1:-}"
CHROM_OPTS=""
if [[ -n "${CHROM_ARG}" ]]; then
  CHROM_OPTS="--chrom ${CHROM_ARG}"
fi

# Optional: assembly gaps + AGP for blue-cross annotation
GAP_OPTS=""
AGP_FILE="${BASE}/00-samples/fClaHyb_Gar_LG.gaps.bed"
AGP_SCAFFOLD="${BASE}/00-samples/fClaHyb_Gar_LG.agp"
if [[ -f "${AGP_FILE}" ]]; then GAP_OPTS="${GAP_OPTS} --gaps ${AGP_FILE}"; fi
if [[ -f "${AGP_SCAFFOLD}" ]]; then GAP_OPTS="${GAP_OPTS} --agp ${AGP_SCAFFOLD}"; fi

echo "================================================================"
echo "  01C: Block detection + boundary classification"
echo "  Chrom: ${CHROM_ARG:-all}"
echo "  Precomp: ${PRECOMP_DIR}"
echo "  Output: ${LANDSCAPE_DIR}"
echo "  JobID: ${SLURM_JOB_ID}"
echo "  Started: $(date)"
echo "================================================================"

"${RSCRIPT_BIN}" "${A2_SNAKES}/PHASE_01C_block_detect.R" \
    "${PRECOMP_DIR}" "${LANDSCAPE_DIR}" \
    --mode hatchery \
    ${CHROM_OPTS} ${GAP_OPTS}

echo "================================================================"
echo "  01C complete — $(date)"
echo "  Landscape: ${LANDSCAPE_DIR}/"
echo "  Next: sbatch LAUNCH_C01b_snake1.slurm --step merge"
echo "================================================================"
