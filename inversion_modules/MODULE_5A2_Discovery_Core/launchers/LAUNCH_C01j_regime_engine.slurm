#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=237G
#SBATCH -t 2:00:00
#SBATCH -J inv_C01j
#SBATCH -o logs/inv_C01j_%j.out
#SBATCH -e logs/inv_C01j_%j.err

# =============================================================================
# LAUNCH_C01j_regime_engine.slurm
#
# PHASE_04C: Regime compatibility engine — compatibility grouping, sample
# types, regime PA, combined all-layers strip plot.
#
# Requires: C01d scoring output + Clair3 VCFs + C01i coseg output.
#
# Usage:
#   sbatch LAUNCH_C01j_regime_engine.slurm
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
SCORING_DIR="${MDS_DIR}/snake_regions_multiscale/scoring"
REGIME_DIR="${MDS_DIR}/snake_regions_multiscale/regime_engine"
COSEG_DIR="${MDS_DIR}/snake_regions_multiscale/decomposition"
VCF_DIR="${BASE}/clair3_pipeline/05_final_catalogs"
SAMPLES_FILE="${HETDIR}/01_inputs_check/samples.ind"

SCORES_FILE="${SCORING_DIR}/candidate_scores.tsv.gz"
if [[ ! -f "${SCORES_FILE}" ]]; then
  echo "[ERROR] Missing scores: ${SCORES_FILE}" >&2
  echo "[ERROR] Run LAUNCH_C01d_scoring.slurm first" >&2
  exit 1
fi

COSEG_OPTS=""
if [[ -d "${COSEG_DIR}" ]]; then COSEG_OPTS="--coseg_dir ${COSEG_DIR}"; fi
SAMPLE_OPTS=""
if [[ -f "${SAMPLES_FILE}" ]]; then SAMPLE_OPTS="--samples ${SAMPLES_FILE}"; fi

echo "================================================================"
echo "  C01j: Regime compatibility engine"
echo "  Scores: ${SCORES_FILE}"
echo "  VCF dir: ${VCF_DIR}"
echo "  Output: ${REGIME_DIR}"
echo "  JobID: ${SLURM_JOB_ID}"
echo "  Started: $(date)"
echo "================================================================"

"${RSCRIPT_BIN}" "${A2_SNAKES}/STEP_C01j_regime_compatibility_engine.R" \
    --scores "${SCORES_FILE}" \
    --vcf_dir "${VCF_DIR}" \
    --outdir "${REGIME_DIR}" \
    ${COSEG_OPTS} ${SAMPLE_OPTS}

echo "================================================================"
echo "  C01j complete — $(date)"
echo "  Output: ${REGIME_DIR}/"
echo "  Next: PHASE_05 confirmation (C01e figures, C01k annotated simmat)"
echo "================================================================"
