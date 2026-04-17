#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -A lt200308
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH -t 1:00:00
#SBATCH -J inv_4e_char
#SBATCH -o logs/inv_4e_char_%j.out
#SBATCH -e logs/inv_4e_char_%j.err

# =============================================================================
# LAUNCH_characterize_classify.sh
#
# Phase 4e — final characterization + Q1–Q7 classification. Reads per-
# candidate keys from evidence_registry and emits the final status / class
# table for the manuscript. Two Rscript calls in sequence:
#
#   1. run_characterize.R        <registry_dir> <outdir>
#       — aggregates per-candidate evidence into the characterize table.
#   2. compute_candidate_status.R <registry_dir> <outdir>
#       — applies the 7-question classification logic (CLASS / UNCLASS /
#         UNCERTAIN / REJECTED) from the 41-key registry spec.
#
# Called by: orchestrator/run_phase4.sh (depends on C01d pass-2).
#
# Upstream:    4e C01d pass-2 (final candidate table)
#              evidence_registry keys.tsv (populated across phases)
# Downstream:  manuscript figures (figrid).
# =============================================================================

set -euo pipefail
source ~/.bashrc
mamba activate assembly

CONFIG="${CONFIG:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_modules/00_inversion_config.sh}"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
set -a
source "${CONFIG}"
set +a

BRIDGE_SH="${BRIDGE_SH:-${BASE}/utils/pipeline_bridge.sh}"
[[ -f "${BRIDGE_SH}" ]] || { echo "Missing bridge (need evidence_registry)" >&2; exit 1; }
source "${BRIDGE_SH}"

inv_init_dirs
mkdir -p logs

# ── Parse orchestrator-style args ────────────────────────────────────────
CHROMS_CSV=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --chroms) CHROMS_CSV="$2"; shift 2 ;;
    --array=*) shift ;;
    --array)   shift 2 ;;
    *) echo "[char] unknown arg (ignored): $1" >&2; shift ;;
  esac
done

# ── Args ─────────────────────────────────────────────────────────────────
REGISTRY_DIR_RUN="${REGISTRY_DIR_RUN:-${BASE}/registries/data/evidence_registry}"
OUTDIR_RUN="${OUTDIR_RUN:-${MDS_DIR}/snake_regions_multiscale/phase4e_final}"
mkdir -p "${OUTDIR_RUN}"

[[ -d "${REGISTRY_DIR_RUN}" ]] || {
  echo "[char] ERROR: evidence_registry dir missing — ${REGISTRY_DIR_RUN}" >&2; exit 2;
}

CHAR_SCRIPT="${BASE}/inversion_modules/phase_4_postprocessing/4e_final_classification/run_characterize.R"
STATUS_SCRIPT="${BASE}/inversion_modules/phase_4_postprocessing/4e_final_classification/compute_candidate_status.R"

for s in "${CHAR_SCRIPT}" "${STATUS_SCRIPT}"; do
  [[ -f "$s" ]] || { echo "[char] ERROR: missing script — $s" >&2; exit 2; }
done

echo "================================================================"
echo "  4e / characterize + classify"
echo "  Registry dir:      ${REGISTRY_DIR_RUN}"
echo "  Out dir:           ${OUTDIR_RUN}"
echo "  RESULTS_REGISTRY:  ${RESULTS_REGISTRY_DIR:-(bridge missing)}"
echo "  SAMPLE_GROUP:      ${SAMPLE_GROUP:-all_226}"
echo "  Chroms:            ${CHROMS_CSV:-all}"
echo "  JobID:             ${SLURM_JOB_ID}"
echo "  Started:           $(date)"
echo "================================================================"

echo "--- 4e.1 characterize ---"
"${RSCRIPT_BIN}" "${CHAR_SCRIPT}" "${REGISTRY_DIR_RUN}" "${OUTDIR_RUN}"

echo ""
echo "--- 4e.2 compute candidate status (Q1–Q7 classification) ---"
"${RSCRIPT_BIN}" "${STATUS_SCRIPT}" "${REGISTRY_DIR_RUN}" "${OUTDIR_RUN}"

echo "================================================================"
echo "  4e complete — $(date)"
echo "  Characterize:   ${OUTDIR_RUN}/characterize_candidates.tsv"
echo "  Status table:   ${OUTDIR_RUN}/candidate_status.tsv"
echo ""
echo "  Integrity check (optional but recommended):"
echo "    Rscript -e 'source(file.path(Sys.getenv(\"BASE\"),"
echo "                \"utils/load_bridge.R\"));"
echo "                ic <- reg\$results\$integrity_check();"
echo "                print(ic); cat(\"all_pass=\", attr(ic, \"all_pass\"), \"\\n\")'"
echo "================================================================"
