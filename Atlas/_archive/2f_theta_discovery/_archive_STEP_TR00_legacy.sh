#!/usr/bin/env bash
# =============================================================================
# STEP_TR00_aggregate_theta.sh — aggregate per-sample pestPG → per-chrom tsv.gz
# =============================================================================
# This is the prerequisite for STEP_TR01. Per-sample θπ already exists from
# MODULE_3 het/ROH (02_run_heterozygosity.sh produces *.win50000.step10000.pestPG
# per sample). This script aggregates them into per-chromosome long-format TSVs
# the next step pivots into wide samples × windows arrays.
#
# This wraps the existing STEP10b_parse_pestPG_to_sample_theta_windows.R logic.
# If that script is the canonical aggregator already, this script can simply
# call it with the right arguments. If it doesn't exist (or is in an old folder
# location), this script implements the parsing inline.
#
# Output:
#   ${THETA_TSV_DIR}/theta.<CHROM>.win50000.step10000.tsv.gz
# Columns:
#   sample  chrom  win_start  win_end  win_center  tP  nSites
#
# Usage:
#   ./STEP_TR00_aggregate_theta.sh                 # all chromosomes
#   ./STEP_TR00_aggregate_theta.sh C_gar_LG28      # one chromosome
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_theta_config.sh"

CHROM_FILTER="${1:-}"

mkdir -p "${THETA_TSV_DIR}"

# Find the existing aggregator (if any). Check three plausible locations
# in priority order:
AGGREGATOR=""
for candidate in \
  "${CODEBASE}/phase_2_discovery/2f_theta_discovery/_lib/parse_pestPG.R" \
  "${CODEBASE}/MODULE_5A2_Discovery_Core/snakes/STEP10b_parse_pestPG_to_sample_theta_windows.R" \
  "${BASE}/het_roh/scripts/STEP10b_parse_pestPG_to_sample_theta_windows.R"
do
  if [[ -f "${candidate}" ]]; then
    AGGREGATOR="${candidate}"
    break
  fi
done

if [[ -z "${AGGREGATOR}" ]]; then
  echo "[ERROR] No pestPG aggregator script found." >&2
  echo "  Checked:" >&2
  echo "    ${CODEBASE}/phase_2_discovery/2f_theta_discovery/_lib/parse_pestPG.R" >&2
  echo "    ${CODEBASE}/MODULE_5A2_Discovery_Core/snakes/STEP10b_parse_pestPG_to_sample_theta_windows.R" >&2
  echo "    ${BASE}/het_roh/scripts/STEP10b_parse_pestPG_to_sample_theta_windows.R" >&2
  echo "  Either copy the existing one to one of those paths, or implement" >&2
  echo "  inline parsing here. The .pestPG header is:" >&2
  echo "    #(indexStart,indexStop)(firstPos,lastPos)(WinStart,WinStop) Chr WinCenter tW tP tF tH tL Tajima fuf fud fayh zeng nSites" >&2
  echo "  Use WinStart/WinStop/WinCenter (not firstPos/lastPos) and the tP column." >&2
  exit 1
fi

echo "[INFO] Using aggregator: ${AGGREGATOR}"

# Build chromosome list for this invocation
if [[ -n "${CHROM_FILTER}" ]]; then
  CHROMS_TO_RUN=("${CHROM_FILTER}")
else
  CHROMS_TO_RUN=("${CHROM_LIST[@]}")
fi

for CHROM in "${CHROMS_TO_RUN[@]}"; do
  OUT_TSV="${THETA_TSV_DIR}/theta.${CHROM}.${SCALE_LABEL}.tsv.gz"

  if [[ -s "${OUT_TSV}" ]]; then
    echo "[SKIP] ${CHROM} — output already exists: ${OUT_TSV}"
    continue
  fi

  echo "[RUN]  ${CHROM} → ${OUT_TSV}"

  "${RSCRIPT}" "${AGGREGATOR}" \
    "${PESTPG_DIR}" \
    "${SAMPLES_IND}" \
    "${THETA_TSV_DIR}/theta.${CHROM}" \
    "${SCALE_LABEL}" \
    "${CHROM}" \
    2> "${LOG_DIR}/STEP_TR00.${CHROM}.log"

  # Aggregator writes <prefix>.sample_tP_windows.tsv.gz; rename to canonical
  if [[ -f "${THETA_TSV_DIR}/theta.${CHROM}.sample_tP_windows.tsv.gz" ]]; then
    mv "${THETA_TSV_DIR}/theta.${CHROM}.sample_tP_windows.tsv.gz" "${OUT_TSV}"
  fi

  if [[ ! -s "${OUT_TSV}" ]]; then
    echo "[FAIL] ${CHROM} — aggregator did not produce expected output. See log:" >&2
    echo "  ${LOG_DIR}/STEP_TR00.${CHROM}.log" >&2
    exit 2
  fi
done

echo "[DONE] STEP_TR00 — ${#CHROMS_TO_RUN[@]} chromosome(s) aggregated."
echo "       Outputs at: ${THETA_TSV_DIR}/theta.*.${SCALE_LABEL}.tsv.gz"
