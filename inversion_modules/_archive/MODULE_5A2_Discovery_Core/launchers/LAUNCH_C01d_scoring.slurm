#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=237G
#SBATCH -t 0:30:00
#SBATCH -J inv_C01d
#SBATCH -o logs/inv_C01d_%j.out
#SBATCH -e logs/inv_C01d_%j.err

# =============================================================================
# LAUNCH_C01d_scoring.slurm
#
# PHASE_03A: Candidate scoring — 10 dimensions + PA aggregation.
# Runs AFTER triangles (C01c). Reads triangle output + optionally cores
# and network dirs for additional dimensions.
#
# Usage:
#   sbatch LAUNCH_C01d_scoring.slurm
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
TRIANGLE_DIR="${MDS_DIR}/snake_regions_multiscale/triangles"
SCORING_DIR="${MDS_DIR}/snake_regions_multiscale/scoring"
CORES_DIR="${MDS_DIR}/snake_regions_multiscale"
NETWORK_DIR="${MDS_DIR}/snake_regions_multiscale/diag_network"

CORES_OPTS=""
if [[ -d "${CORES_DIR}" ]]; then CORES_OPTS="--cores_dir ${CORES_DIR}"; fi
NET_OPTS=""
if [[ -d "${NETWORK_DIR}" ]]; then NET_OPTS="--network_dir ${NETWORK_DIR}"; fi

echo "================================================================"
echo "  C01d: Candidate scoring (10D + PA aggregation)"
echo "  Triangles: ${TRIANGLE_DIR}"
echo "  Output: ${SCORING_DIR}"
echo "  JobID: ${SLURM_JOB_ID}"
echo "  Started: $(date)"
echo "================================================================"

"${RSCRIPT_BIN}" "${A2_SNAKES}/STEP_C01d_candidate_scoring.R" \
    "${TRIANGLE_DIR}" "${SCORING_DIR}" \
    ${CORES_OPTS} ${NET_OPTS}

echo "================================================================"
echo "  C01d complete — $(date)"
echo "  Output: ${SCORING_DIR}/"
echo "  Next: PHASE_04 characterization (C01i, C01j, etc.)"
echo "================================================================"
