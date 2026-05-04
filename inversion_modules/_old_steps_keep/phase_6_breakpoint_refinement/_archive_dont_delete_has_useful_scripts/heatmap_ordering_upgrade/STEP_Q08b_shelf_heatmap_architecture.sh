#!/bin/bash
# =============================================================================
# STEP_Q08b_shelf_heatmap_architecture.sh
# =============================================================================
# Architecture-aware version of STEP_Q08_shelf_heatmap.sh.
#
# Backward compatible: if no STEP40 outputs exist, produces identical results
# to the original Q08 script (A/B/C panels, unchanged). When STEP40 outputs
# are present, adds:
#
#   - Architecture label overlay on every panel title
#   - Classifier feature line under each panel
#   - Subgroup annotation column on each heatmap
#   - D1/D2/D3 separate-first locally-harmonized panels (for composite /
#     unresolved candidates only)
#
# Auto-discovery logic: the script looks for STEP40 outputs in the standard
# follow-up directory and passes them along only if found. Explicit overrides
# are possible via environment variables:
#
#   Q08B_ARCH_FILE     explicit path to candidate_architecture_class.tsv
#   Q08B_SUBGROUP_FILE explicit path to candidate_sample_subgroups.tsv
#   Q08B_POLARITY_FILE explicit path to candidate_marker_polarity.tsv
#   Q08B_FOLLOWUP_DIR  explicit follow-up root (default: ${FOLLOWUP_DIR} from
#                      00_config.sh, if that variable exists)
#
# Usage:
#   SHELF_START_MB=15 SHELF_END_MB=18 bash STEP_Q08b_shelf_heatmap_architecture.sh C_gar_LG28
#
# For the catalog-wide run (all chromosomes), loop externally as before.
# =============================================================================
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/00_config.sh"

CHR="${1:-}"
[[ -z "${CHR}" ]] && qc_die "Usage: $0 <CHR>"
[[ -z "${SHELF_START_MB:-}" || -z "${SHELF_END_MB:-}" ]] && \
  qc_die "SHELF_START_MB and SHELF_END_MB required"

# ---------------------------------------------------------------------------
# Discover optional STEP40 follow-up outputs
# ---------------------------------------------------------------------------
# The STEP40 classifier writes per-candidate outputs to:
#   <FOLLOWUP_DIR>/<CHR>.candidate_<cid>/candidate_architecture_class.tsv
#   <FOLLOWUP_DIR>/<CHR>.candidate_<cid>/candidate_sample_subgroups.tsv
#   <FOLLOWUP_DIR>/<CHR>.candidate_<cid>/candidate_marker_polarity.tsv
#
# And a catalog at:
#   <FOLLOWUP_DIR>/architecture_catalog.tsv
#
# For Q08 we are typically dealing with a *single* candidate per chromosome
# (the shelf region). We first honor explicit overrides, then look under
# FOLLOWUP_DIR if that variable is defined, then fall back to searching a
# conventional location.
# ---------------------------------------------------------------------------

discover_followup() {
  local chr="$1"
  local root=""
  if [[ -n "${Q08B_FOLLOWUP_DIR:-}" ]]; then
    root="${Q08B_FOLLOWUP_DIR}"
  elif [[ -n "${FOLLOWUP_DIR:-}" ]]; then
    root="${FOLLOWUP_DIR}"
  elif [[ -n "${QC_ROOT:-}" && -d "${QC_ROOT}/followup" ]]; then
    root="${QC_ROOT}/followup"
  fi

  if [[ -z "${root}" || ! -d "${root}" ]]; then
    echo ""
    return 0
  fi

  # Prefer the candidate subfolder that matches the current chromosome.
  # Candidates are named <chrom>.candidate_<cid>; with multiple, take the
  # one with the architecture_class.tsv present.
  local match
  match=$(find "${root}" -maxdepth 2 -type d -name "${chr}.candidate_*" 2>/dev/null | head -n 1)
  echo "${match}"
}

CANDIDATE_DIR="$(discover_followup "${CHR}")"
ARCH_FILE="${Q08B_ARCH_FILE:-}"
SUBG_FILE="${Q08B_SUBGROUP_FILE:-}"
POL_FILE="${Q08B_POLARITY_FILE:-}"

if [[ -n "${CANDIDATE_DIR}" ]]; then
  qc_log "Q08b: discovered candidate dir = ${CANDIDATE_DIR}"
  [[ -z "${ARCH_FILE}" && -f "${CANDIDATE_DIR}/candidate_architecture_class.tsv" ]] && \
    ARCH_FILE="${CANDIDATE_DIR}/candidate_architecture_class.tsv"
  [[ -z "${SUBG_FILE}" && -f "${CANDIDATE_DIR}/candidate_sample_subgroups.tsv" ]] && \
    SUBG_FILE="${CANDIDATE_DIR}/candidate_sample_subgroups.tsv"
  [[ -z "${POL_FILE}"  && -f "${CANDIDATE_DIR}/candidate_marker_polarity.tsv" ]] && \
    POL_FILE="${CANDIDATE_DIR}/candidate_marker_polarity.tsv"
fi

run_one() {
  local chr="$1"
  qc_log "Q08b ${chr}: architecture-aware shelf genotype heatmap"

  local beagle="${BEAGLE_DIR}/main_qcpass.${chr}.beagle.gz"
  [[ -f "${beagle}" ]] || { qc_log "SKIP ${chr}: no beagle"; return 0; }
  local pos="${BEAGLE_DIR}/main_qcpass.${chr}.pos.fixed"
  [[ -f "${pos}" ]]    || { qc_log "SKIP ${chr}: no pos file"; return 0; }
  local precomp="${PRECOMP_DIR}/${chr}.precomp.rds"
  [[ -f "${precomp}" ]] || precomp="${PRECOMP_DIR}/main_qcpass.${chr}.precomp.rds"
  [[ -f "${precomp}" ]] || { qc_log "SKIP ${chr}: no precomp"; return 0; }

  local out_pdf="${QC_FIGS}/shelf_heatmap_arch.${chr}.pdf"
  local out_dir="${QC_FIGS}/shelf_heatmap_arch.${chr}"
  local invgt_assign="${QC_TRACKS}/invgt_assignments.${chr}.tsv"

  local args=(
    --beagle  "${beagle}"
    --pos     "${pos}"
    --precomp "${precomp}"
    --chrom   "${chr}"
    --shelf_start_mb "${SHELF_START_MB}"
    --shelf_end_mb   "${SHELF_END_MB}"
    --sample_list    "${SAMPLE_LIST_POPSTATS}"
    --out     "${out_pdf}"
    --out_dir "${out_dir}"
  )
  [[ -n "${BP1_MB:-}" ]] && args+=( --breakpoint1_mb "${BP1_MB}" )
  [[ -n "${BP2_MB:-}" ]] && args+=( --breakpoint2_mb "${BP2_MB}" )
  [[ -f "${invgt_assign}" ]] && args+=( --invgt_assign "${invgt_assign}" )

  # Architecture-aware extras
  if [[ -n "${ARCH_FILE}" && -f "${ARCH_FILE}" ]]; then
    args+=( --arch_class "${ARCH_FILE}" )
    qc_log "  -> arch class file: ${ARCH_FILE}"
  fi
  if [[ -n "${SUBG_FILE}" && -f "${SUBG_FILE}" ]]; then
    args+=( --subgroup_file "${SUBG_FILE}" )
    qc_log "  -> subgroup file : ${SUBG_FILE}"
  fi
  if [[ -n "${POL_FILE}" && -f "${POL_FILE}" ]]; then
    args+=( --polarity_file "${POL_FILE}" )
    qc_log "  -> polarity file : ${POL_FILE}"
  fi
  if [[ -z "${ARCH_FILE}" && -z "${SUBG_FILE}" && -z "${POL_FILE}" ]]; then
    qc_log "  (no STEP40 outputs found — running in legacy mode)"
  fi

  ${RSCRIPT_BIN} --vanilla "${here}/R/q08b_shelf_heatmap_architecture.R" "${args[@]}" \
    2> "${QC_LOGS}/q08b_${chr}.log"

  [[ -f "${out_pdf}" ]] && qc_log "  -> ${out_pdf}"
}

run_one "${CHR}"
qc_log "Q08b DONE."
