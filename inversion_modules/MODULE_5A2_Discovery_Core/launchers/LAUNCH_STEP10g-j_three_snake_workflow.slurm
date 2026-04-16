#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH -J three_snake
#SBATCH -o logs/STEP10g-j_three_snake_%j.out
#SBATCH -e logs/STEP10g-j_three_snake_%j.err

# =============================================================================
# LAUNCH_STEP10g-j_three_snake_workflow.slurm
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

SCRIPT_DIR="${DISCOVERY2DIR}/launchers"
A2_SNAKES="${DISCOVERY2DIR}/snakes"
A2_HMM="${DISCOVERY2DIR}/hmm"

SKIP_HMM="${1:-no}"

SNAKE_DIR="${INVDIR}/06_mds_candidates/snake_regions"
SNAKE2_DIR="${INVDIR}/06_mds_candidates/snake2_community"
SNAKE3_DIR="${INVDIR}/06_mds_candidates/snake3_ghsl"
CONSENSUS_DIR="${INVDIR}/06_mds_candidates/three_snake_consensus"
HMM_DIR="${INVDIR}/06_mds_candidates/regime_hmm"
PHASED_DIR="${INVDIR}/phased_summaries"

mkdir -p "${SNAKE2_DIR}" "${SNAKE3_DIR}" "${CONSENSUS_DIR}" "${HMM_DIR}" "${INVDIR}/logs"

echo "================================================================"
echo "  3-Snake Workflow (STEP10g → STEP10j)"
echo "  Started: $(date)"
echo "================================================================"

inv_log "STEP10g: Snake 2 — community continuity"
"${RSCRIPT_BIN}" "${A2_SNAKES}/STEP10g_snake2_community_continuity.R" \
    "${MDS_DIR}/inversion_localpca" "${DOSAGE_DIR}" "${SAMPLES_IND}" "${SNAKE2_DIR}"

if [[ -d "${PHASED_DIR}" ]]; then
  SNAKE3_V2="${DISCOVERY2DIR}/ghsl/STEP10h_snake3_ghsl_haplotype_contrast_v2.R"
  SNAKE3_V1="${A2_SNAKES}/STEP10h_snake3_ghsl_haplotype_contrast.R"

  if [[ -f "${SNAKE3_V2}" ]]; then
    inv_log "STEP10h: Snake 3 v2 — GHSL haplotype contrast (adaptive thresholds)"
    "${RSCRIPT_BIN}" "${SNAKE3_V2}" \
        "${MDS_DIR}/inversion_localpca" "${PHASED_DIR}" "${SAMPLES_IND}" "${SNAKE3_DIR}"
  else
    inv_log "STEP10h: Snake 3 v1 — GHSL-inspired haplotype contrast"
    "${RSCRIPT_BIN}" "${SNAKE3_V1}" \
        "${MDS_DIR}/inversion_localpca" "${PHASED_DIR}" "${SAMPLES_IND}" "${SNAKE3_DIR}"
  fi
else
  inv_log "STEP10h: SKIP — phased summary dir not found: ${PHASED_DIR}"
  inv_log "Snake 3 will emit FAIL for all windows (consensus still works with 2 snakes)"
fi

inv_log "STEP10i: 3-snake overlap consensus"
"${RSCRIPT_BIN}" "${A2_SNAKES}/STEP10i_three_snake_consensus.R" \
    "${MDS_DIR}/inversion_localpca" "${SNAKE_DIR}" "${SNAKE2_DIR}" \
    "${SNAKE3_DIR}" "${CONSENSUS_DIR}"

if [[ "${SKIP_HMM}" != "skip_hmm" ]]; then
  inv_log "STEP10j: Regime HMM"
  "${RSCRIPT_BIN}" "${A2_HMM}/STEP10j_regime_hmm.R" \
      "${MDS_DIR}/inversion_localpca" "${SNAKE_DIR}" "${SNAKE2_DIR}" \
      "${SNAKE3_DIR}" "${CONSENSUS_DIR}" "${HMM_DIR}"
else
  inv_log "STEP10j: SKIP (skip_hmm requested)"
fi

echo "================================================================"
echo "  3-Snake Workflow complete — $(date)"
echo "================================================================"
echo ""
echo "Outputs:"
echo "  Snake 2:    ${SNAKE2_DIR}/"
echo "  Snake 3:    ${SNAKE3_DIR}/"
echo "  Consensus:  ${CONSENSUS_DIR}/"
echo "  HMM:        ${HMM_DIR}/"
