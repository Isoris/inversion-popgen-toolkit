#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=237G
#SBATCH -t 1:00:00
#SBATCH -J inv_C01c
#SBATCH -o logs/inv_C01c_%j.out
#SBATCH -e logs/inv_C01c_%j.err

# =============================================================================
# LAUNCH_C01c_triangles.slurm
#
# PHASE_02C: Triangle regime detection — insulation boundaries, squareness
# check, hatchery mode. Runs AFTER merge (C01b_2). Reads precomp.
#
# Usage:
#   sbatch LAUNCH_C01c_triangles.slurm               # all chromosomes
#   sbatch LAUNCH_C01c_triangles.slurm C_gar_LG28    # single chromosome
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
TRIANGLE_DIR="${MDS_DIR}/snake_regions_multiscale/triangles"

CHROM_ARG="${1:-}"
CHROM_OPTS=""
if [[ -n "${CHROM_ARG}" ]]; then
  CHROM_OPTS="--chrom ${CHROM_ARG}"
fi

# Auto-detect ancestry Q-matrix and relatedness for overlay
ANC_OPTS=""
Q_FILE="${MDS_DIR}/../07_round1_structure/best_seed_K8.qopt"
REL_FILE="${MDS_DIR}/../05_ngsrelate/catfish_226_for_natora.txt"
SAMPLES_FILE="${HETDIR}/01_inputs_check/samples.ind"
if [[ -f "${Q_FILE}" ]]; then ANC_OPTS="${ANC_OPTS} --ancestry ${Q_FILE}"; fi
if [[ -f "${REL_FILE}" ]]; then ANC_OPTS="${ANC_OPTS} --relatedness ${REL_FILE}"; fi
if [[ -f "${SAMPLES_FILE}" ]]; then ANC_OPTS="${ANC_OPTS} --samples ${SAMPLES_FILE}"; fi

echo "================================================================"
echo "  C01c: Triangle regimes (hatchery mode)"
echo "  Chrom: ${CHROM_ARG:-all}"
echo "  Precomp: ${PRECOMP_DIR}"
echo "  Output: ${TRIANGLE_DIR}"
echo "  JobID: ${SLURM_JOB_ID}"
echo "  Started: $(date)"
echo "================================================================"

"${RSCRIPT_BIN}" "${A2_SNAKES}/STEP_C01c_triangle_regimes.R" \
    "${PRECOMP_DIR}" "${TRIANGLE_DIR}" \
    --mode hatchery \
    ${CHROM_OPTS} ${ANC_OPTS}

echo "================================================================"
echo "  C01c complete — $(date)"
echo "  Output: ${TRIANGLE_DIR}/"
echo "  Next: sbatch LAUNCH_C01d_scoring.slurm"
echo "================================================================"
