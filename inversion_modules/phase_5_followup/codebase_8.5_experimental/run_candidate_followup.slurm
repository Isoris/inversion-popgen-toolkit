#!/bin/bash
#SBATCH --job-name=inv_followup
#SBATCH --output=logs/inv_followup_%A_%a.out
#SBATCH --error=logs/inv_followup_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=compute

set -euo pipefail

source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_DIR="${SCRIPT_DIR}"
CONFIG="${MODULE_DIR}/config_inversion_followup.R"
RSCRIPT_BIN="/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript"

mkdir -p logs

[[ -x "${RSCRIPT_BIN}" ]] || { echo "[ERROR] Rscript not executable: ${RSCRIPT_BIN}" >&2; exit 1; }
[[ -s "${CONFIG}" ]] || { echo "[ERROR] Missing config file: ${CONFIG}" >&2; exit 1; }

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    CID="${SLURM_ARRAY_TASK_ID}"
elif [[ -n "${1:-}" ]]; then
    CID="$1"
else
    CID="all"
fi

MODE="${2:-full}"

echo "=========================================="
echo "Inversion candidate follow-up pipeline v6"
echo "Candidate: ${CID}"
echo "Mode: ${MODE}"
echo "Start: $(date)"
echo "Rscript: ${RSCRIPT_BIN}"
echo "Module dir: ${MODULE_DIR}"
echo "Config: ${CONFIG}"
echo "=========================================="

run_step() {
  local label="$1"
  local script="$2"
  shift 2
  echo "[${label}] Running ${script}..."
  "${RSCRIPT_BIN}" "${MODULE_DIR}/${script}" "$@"
}

if [[ "${MODE}" == "full" || "${MODE}" == "core" ]]; then
  run_step STEP20 STEP20_candidate_window_analysis.R           "${CONFIG}" "${CID}"
  run_step STEP21 STEP21_candidate_state_assignment.R          "${CONFIG}" "${CID}"
  run_step STEP22 STEP22_candidate_subclustering.R             "${CONFIG}" "${CID}"
  run_step STEP23 STEP23_candidate_ancestry_association.R      "${CONFIG}" best "${CID}"
  run_step STEP24 STEP24_candidate_het_and_marker_support.R    "${CONFIG}" "${CID}"
  run_step STEP25 STEP25_candidate_interpretation_table.R      "${CONFIG}" "${CID}"
fi

if [[ "${MODE}" == "full" || "${MODE}" == "extended" ]]; then
  run_step STEP29 STEP29_candidate_coherence_and_polarity.R      "${CONFIG}" "${CID}"
  run_step STEP30 STEP30_candidate_window_trajectories.R         "${CONFIG}" "${CID}"
  run_step STEP32 STEP32_candidate_multiscale_windows.R          "${CONFIG}" "${CID}"
  run_step STEP33 STEP33_candidate_within_stripe_analysis.R      "${CONFIG}" "${CID}"
  run_step STEP34 STEP34_candidate_anchor_resolution.R           "${CONFIG}" "${CID}"
  run_step STEP36 STEP36_candidate_clair3_local_signatures.R     "${CONFIG}" "${CID}"
  run_step STEP25B STEP25B_candidate_interpretation_table_v2.R   "${CONFIG}" "${CID}"
fi

if [[ "${MODE}" == "full" || "${MODE}" == "plots" ]]; then
  run_step STEP26 STEP26_candidate_figure_catalogue.R            "${CONFIG}" "${CID}" intermediate
  run_step STEP27 STEP27_candidate_faceted_local_pca.R           "${CONFIG}" "${CID}"
  run_step STEP28 STEP28_candidate_marker_heatmap.R              "${CONFIG}" "${CID}"
  run_step STEP31 STEP31_candidate_diagnostic_figures.R          "${CONFIG}" "${CID}"
  run_step STEP35 STEP35_candidate_harmonized_heatmaps.R         "${CONFIG}" "${CID}"
fi

echo "=========================================="
echo "Pipeline complete for candidate: ${CID}"
echo "Mode: ${MODE}"
echo "End: $(date)"
echo "=========================================="
