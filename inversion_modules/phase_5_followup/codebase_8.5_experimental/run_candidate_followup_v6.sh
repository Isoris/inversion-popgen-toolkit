#!/bin/bash
set -euo pipefail

source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source shared config if available; fall back to local derivation
if [[ -f "${SCRIPT_DIR}/../00_inversion_config.sh" ]]; then
  source "${SCRIPT_DIR}/../00_inversion_config.sh"
fi

CONFIG="${SCRIPT_DIR}/config_inversion_followup.R"
RSCRIPT_BIN="${RSCRIPT_BIN:-/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript}"
DISCOVERYDIR="${DISCOVERYDIR:-${SCRIPT_DIR}/../MODULE_5A2_Discovery_Core}"

CID="${1:-all}"
MODE="${2:-full}"

[[ -x "${RSCRIPT_BIN}" ]] || { echo "[ERROR] Rscript not executable: ${RSCRIPT_BIN}" >&2; exit 1; }
[[ -s "${CONFIG}" ]] || { echo "[ERROR] Missing config file: ${CONFIG}" >&2; exit 1; }

echo "================================================================"
echo "  MODULE_5B Inversion Follow-up Pipeline v7.4"
echo "  Module dir   : ${SCRIPT_DIR}"
echo "  Discovery dir: ${DISCOVERYDIR}"
echo "  Config       : ${CONFIG}"
echo "  Rscript      : ${RSCRIPT_BIN}"
echo "  Candidate    : ${CID}"
echo "  Mode         : ${MODE}"
echo "  Started      : $(date)"
echo "================================================================"

run() {
  local label="$1" script="$2"; shift 2
  echo ""
  echo "──── ${label} ──────────────────────────────────────────"
  if [[ -f "${SCRIPT_DIR}/${script}" ]]; then
    "${RSCRIPT_BIN}" "${SCRIPT_DIR}/${script}" "$@"
  else
    echo "[SKIP] ${script} not found"
  fi
}

if [[ "${MODE}" == "full" || "${MODE}" == "core" ]]; then
  # Window sample-belonging analysis (discovery-level, goes before STEP20)
  run STEP10c "${DISCOVERYDIR}/steps/STEP_D02_window_sample_belonging.R" "${CONFIG}" "${CID}"
  # Layer comparison QC gate (automatic after STEP10c)
  run STEP10d "${DISCOVERYDIR}/steps/STEP_D03_layer_comparison_qc.R" "${CONFIG}" "${CID}"
  # New sample-structure-first engine (v3)
  run STEP20  current_followup/STEP20_candidate_window_profile_engine.R  "${CONFIG}" "${CID}"
  run STEP20b current_followup/STEP20b_anchor_geometry_module.R          "${CONFIG}" "${CID}"
  run STEP20c current_followup/STEP20c_gradient_marker_experiments.R     "${CONFIG}" "${CID}"
  # Existing follow-up (consumes STEP12 outputs)
  run STEP21 current_followup/STEP21_candidate_state_assignment.R        "${CONFIG}" "${CID}"
  run STEP22 current_followup/STEP22_candidate_subclustering.R           "${CONFIG}" "${CID}"
  run STEP23 current_followup/STEP23_candidate_ancestry_association.R    "${CONFIG}" best "${CID}"
  run STEP24 current_followup/STEP24_candidate_het_and_marker_support.R  "${CONFIG}" "${CID}"
  run STEP25 current_followup/STEP25B_candidate_interpretation_table_v2.R "${CONFIG}" "${CID}"
fi

if [[ "${MODE}" == "full" || "${MODE}" == "extended" ]]; then
  run STEP29 current_followup/STEP29_candidate_coherence_and_polarity.R      "${CONFIG}" "${CID}"
  run STEP30 current_followup/STEP30_candidate_window_trajectories.R         "${CONFIG}" "${CID}"
  run STEP32 current_followup/STEP32_candidate_multiscale_windows.R          "${CONFIG}" "${CID}"
  run STEP33 current_followup/STEP33_candidate_within_stripe_analysis.R      "${CONFIG}" "${CID}"
  run STEP34 current_followup/STEP34_candidate_anchor_resolution.R           "${CONFIG}" "${CID}"
  run STEP36 current_followup/STEP36_candidate_clair3_local_signatures.R     "${CONFIG}" "${CID}"
  # Breakpoint support from DELLY2 + Manta SV calls (requires --delly_dir / --manta_dir)
  DELLY_SV_DIR="${BASE}/delly_sv_discovery/per_sample_merged"
  MANTA_SV_DIR="${BASE}/manta_sv_discovery/per_sample_merged"
  if [[ -d "${DELLY_SV_DIR}" || -d "${MANTA_SV_DIR}" ]]; then
    run STEP37 current_followup/STEP37_breakpoint_support_from_sv.R "${CONFIG}" "${CID}" \
      --delly_dir "${DELLY_SV_DIR}" --manta_dir "${MANTA_SV_DIR}"
  else
    echo "[SKIP] STEP37: No DELLY/Manta SV directories found"
  fi
fi

if [[ "${MODE}" == "full" || "${MODE}" == "plots" ]]; then
  run STEP26 current_followup/STEP26_candidate_figure_catalogue.R            "${CONFIG}" "${CID}" intermediate
  run STEP27 current_followup/STEP27_candidate_faceted_local_pca.R           "${CONFIG}" "${CID}"
  run STEP28 current_followup/STEP28_candidate_marker_heatmap.R              "${CONFIG}" "${CID}"
  run STEP31 current_followup/STEP31_candidate_diagnostic_figures.R          "${CONFIG}" "${CID}"
  run STEP35 current_followup/STEP35_candidate_harmonized_heatmaps.R         "${CONFIG}" "${CID}"
  # Composite inversion figure (requires STEP37 breakpoint support + snake tracks)
  run STEP38 current_followup/STEP38_inversion_composite_figure.R            "${CONFIG}" "${CID}"
  # Breakpoint visualization (population overview + focal highlighted)
  DELLY_SV_DIR="${BASE}/delly_sv_discovery/per_sample_merged"
  run STEP39 current_followup/STEP39_breakpoint_visualization.R              "${CONFIG}" --sv_dir "${DELLY_SV_DIR}" --cid "${CID}"
  # Internal coherence diagnostic (composite candidate detection)
  run STEP40 current_followup/STEP40_candidate_internal_coherence.R          "${CONFIG}" "${CID}"
  # Membership trajectories (H1 stable sub-variants vs H2 positional changes)
  run STEP41 current_followup/STEP41_candidate_membership_trajectories.R     "${CONFIG}" "${CID}"
fi

echo ""
echo "================================================================"
echo "  Pipeline v7.4 complete — Candidate: ${CID}"
echo "  Finished : $(date)"
echo "================================================================"
