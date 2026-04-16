#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=237G
#SBATCH -t 1:00:00
#SBATCH -J inv_C01a
#SBATCH -o logs/inv_C01a_%j.out
#SBATCH -e logs/inv_C01a_%j.err

# =============================================================================
# LAUNCH_C01a_snake1_precompute.slurm
#
# Precompute step — run once. Loads .mds.rds, builds per-chromosome
# similarity matrices, NN distances, robust z-scores, inv-likeness.
# Saves .precomp.rds files for fast C01b reruns.
#
# Uses 28 cores for parallel per-chromosome precompute via mclapply.
#
# Usage:
#   sbatch LAUNCH_C01a_snake1_precompute.slurm
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

echo "================================================================"
echo "  C01a: Snake 1 precompute (parallel, ${SLURM_CPUS_PER_TASK} cores)"
echo "  JobID: ${SLURM_JOB_ID}"
echo "  Started: $(date)"
echo "================================================================"

"${RSCRIPT_BIN}" "${A2_SNAKES}/STEP_C01a_snake1_precompute.R" \
    "${MDS_DIR}/inversion_localpca" \
    "${MDS_DIR}/snake_regions_multiscale"

echo "================================================================"
echo "  C01a complete — $(date)"
echo "  Precomp: ${MDS_DIR}/snake_regions_multiscale/precomp/"
echo "  Next: sbatch LAUNCH_C01b_snake1_run.slurm"
echo "================================================================"

