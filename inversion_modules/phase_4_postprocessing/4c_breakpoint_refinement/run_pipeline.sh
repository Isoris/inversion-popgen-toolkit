#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh  (v1.1 — registry-wired)
#
# Runs the inversion breakpoint pipeline for one candidate (or all candidates
# registered in interval_registry). Chains the four core scripts 01 → 02 →
# 03 → 04 in order. Optionally runs 05 (four-encoding diagnostic) if
# RUN_05=1 is set.
#
# Registry wiring (chat-18):
#   - Sources utils/pipeline_bridge.sh to get REGISTRY_BRIDGE into env.
#   - Candidates come from registries/data/interval_registry/candidate_intervals.tsv
#     (not a flat CANDIDATE_TABLE path).
#   - Each script reads/writes via reg$* — no FOLLOWUP_DIR needed.
#   - config.R is now OPTIONAL and used for parameter overrides only
#     (BP01_PARAMS_OVERRIDE, BP02_PARAMS_OVERRIDE, etc.), not for paths.
#
# 06 (stream graph) and 07 (panel D pileup) are standalone visualizations;
# call them directly with --candidate <cid> when wanted.
#
# Usage:
#   run_pipeline.sh                       # all candidates
#   run_pipeline.sh LG28_1                # just LG28_1
#   run_pipeline.sh LG28_1 my_overrides.R # LG28_1 with parameter overrides
#
# Environment variables:
#   BASE             repo root (auto-detected if SCRIPT_DIR is inside it)
#   REGISTRY_BRIDGE  path to utils/registry_bridge.R (auto-detected)
#   DOSAGE_DIR       where <chrom>.dosage.tsv.gz / <chrom>.sites.tsv.gz live
#                    (defaults to $BASE/popstruct_thin/04_beagle_byRF_majmin)
#   RUN_05=1         also run 05 four-encoding diagnostic
#   LOG_DIR          log output dir (default ./logs/breakpoint_pipeline)
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- Locate BASE (repo root) ---
# The scripts live in <BASE>/inversion_modules/phase_4_postprocessing/4c_breakpoint_refinement/ or
# similar. Walk up until we find utils/pipeline_bridge.sh, or use $BASE if set.
if [[ -z "${BASE:-}" ]]; then
  _dir="${SCRIPT_DIR}"
  for _ in 1 2 3 4 5; do
    if [[ -f "${_dir}/utils/pipeline_bridge.sh" ]]; then
      BASE="${_dir}"
      break
    fi
    _dir="$(dirname "${_dir}")"
  done
fi
if [[ -z "${BASE:-}" ]] || [[ ! -f "${BASE}/utils/pipeline_bridge.sh" ]]; then
  echo "ERROR: cannot locate BASE (repo root containing utils/pipeline_bridge.sh)." >&2
  echo "Set BASE=/path/to/inversion-popgen-toolkit and retry." >&2
  exit 2
fi
export BASE

# --- Source the bash bridge to populate REGISTRY_BRIDGE etc. ---
# shellcheck disable=SC1091
source "${BASE}/utils/pipeline_bridge.sh"

LOG_DIR="${LOG_DIR:-./logs/breakpoint_pipeline}"
mkdir -p "${LOG_DIR}"

CID="${1:-all}"
CONFIG="${2:-}"

# Validate config if given (it's now optional)
if [[ -n "${CONFIG}" && ! -f "${CONFIG}" ]]; then
  echo "ERROR: config file not found: ${CONFIG}" >&2
  exit 2
fi

TS=$(date +%Y%m%d_%H%M%S)
LOG="${LOG_DIR}/pipeline_${CID}_${TS}.log"
{
  echo "[pipeline] cid=${CID}"
  echo "[pipeline] config=${CONFIG:-<none>}"
  echo "[pipeline] BASE=${BASE}"
  echo "[pipeline] REGISTRY_BRIDGE=${REGISTRY_BRIDGE:-<unset>}"
  echo "[pipeline] log: ${LOG}"
} | tee -a "${LOG}"

run_step() {
  local step_name="$1"
  local script_path="$2"
  echo "[pipeline] >>> ${step_name} ..." | tee -a "${LOG}"
  local start_time
  start_time=$(date +%s)

  # Build arg list: optional --config for overrides; then cid (or "all")
  local -a script_args=()
  if [[ -n "${CONFIG}" ]]; then
    script_args+=("--config" "${CONFIG}")
  fi
  script_args+=("${CID}")

  if ! Rscript "${script_path}" "${script_args[@]}" >> "${LOG}" 2>&1; then
    echo "[pipeline] !!! ${step_name} FAILED — see ${LOG}" | tee -a "${LOG}"
    tail -n 30 "${LOG}"
    exit 1
  fi
  local elapsed=$(( $(date +%s) - start_time ))
  echo "[pipeline] <<< ${step_name} OK (${elapsed}s)" | tee -a "${LOG}"
}

# Core pipeline
run_step "01 dosage_signal"       "${SCRIPT_DIR}/01_dosage_signal.R"
run_step "02 ancestral_fragments" "${SCRIPT_DIR}/02_ancestral_fragments.R"
run_step "03 consensus_merge"     "${SCRIPT_DIR}/03_consensus_merge.R"
run_step "04 diagnostic_figure"   "${SCRIPT_DIR}/04_diagnostic_figure.R"

# Optional diagnostic — set RUN_05=1 to include
if [[ "${RUN_05:-0}" == "1" ]]; then
  run_step "05 four_encoding (optional)" "${SCRIPT_DIR}/05_four_encoding_diagnostic.R"
fi

{
  echo "[pipeline] ============================================="
  echo "[pipeline] DONE for cid=${CID}"
  echo "[pipeline] Outputs per candidate are now in the registry:"
  echo "[pipeline]   registries/data/evidence_registry/per_candidate/<cid>/"
  echo "[pipeline]     structured/dosage_blocks.json                (from 01)"
  echo "[pipeline]     structured/ancestral_fragments_summary.json  (from 02)"
  echo "[pipeline]     structured/boundary_refined_left.json        (from 03)  <-- MAIN RESULT"
  echo "[pipeline]     structured/boundary_refined_right.json       (from 03)  <-- MAIN RESULT"
  echo "[pipeline]     structured/breakpoints_per_method.json       (from 03)"
  echo "[pipeline]     figures/breakpoint_diagnostic.pdf            (from 04)"
  if [[ "${RUN_05:-0}" == "1" ]]; then
    echo "[pipeline]     structured/encoding_robustness.json          (from 05)"
    echo "[pipeline]     figures/encoding/<label>_heatmaps_panel.pdf  (from 05)"
  fi
  echo "[pipeline]     raw/dosage_informative_markers.tsv.gz        (from 01)"
  echo "[pipeline]     raw/ancestral_fragments_per_sample.tsv.gz    (from 02)"
  echo "[pipeline] "
  echo "[pipeline] Quick query for refined breakpoint, left side:"
  echo "[pipeline]   jq '.data.final_bp' registries/data/evidence_registry/per_candidate/${CID}/structured/boundary_refined_left.json"
  echo "[pipeline] ============================================="
} | tee -a "${LOG}"
