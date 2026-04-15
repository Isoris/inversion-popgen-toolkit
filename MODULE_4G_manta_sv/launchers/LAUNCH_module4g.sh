#!/usr/bin/env bash
# =============================================================================
# LAUNCH_module4g.sh — Submit Manta SV pipeline to SLURM with dependencies
# =============================================================================
# Usage:
#   cd /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/MODULE_4G_ALL_Manta
#   bash run_manta_pipeline.sh
# =============================================================================
set -euo pipefail

# ── Config loading (same pattern as SLURM scripts) ─────────────────────────
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${SCRIPT_DIR}/../00_module4g_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }

set -a
source "${CONFIG}"
set +a

# ── Create log directory ───────────────────────────────────────────────────
LOG_DIR="${OUTDIR}/logs"
mkdir -p "${LOG_DIR}"

echo "=== Manta SV Pipeline Launcher ==="
echo "Script directory: ${SCRIPT_DIR}"
echo "Manta install:   ${MANTA_INSTALL}"
echo "BAM manifest:    ${BAM_MANIFEST}"
echo "Custom ini:      minCandidateVariantSize = 50"
echo ""

# Verify Manta exists
[[ -f "${MANTA_CONFIG}" ]] || { echo "ERROR: configManta.py not found: ${MANTA_CONFIG}" >&2; exit 1; }

# Verify manifest exists
[[ -f "${BAM_MANIFEST}" ]] || { echo "ERROR: Manifest not found: ${BAM_MANIFEST}" >&2; exit 1; }

echo "Submitting 5-job dependency chain..."
echo ""

JID1=$(sbatch --parsable \
  -o "${LOG_DIR}/01_prep.%j.out" -e "${LOG_DIR}/01_prep.%j.err" \
  "${SCRIPT_DIR}/../slurm/SLURM_A01_prep_inputs.sh")
echo "  [1/5] Prep inputs:       Job ${JID1}"

JID2=$(sbatch --parsable --dependency=afterok:${JID1} \
  -o "${LOG_DIR}/02_discovery.%j.out" -e "${LOG_DIR}/02_discovery.%j.err" \
  "${SCRIPT_DIR}/../slurm/SLURM_A02_manta_discovery.sh")
echo "  [2/5] Manta discovery:   Job ${JID2} (after ${JID1})"

JID3=$(sbatch --parsable --dependency=afterok:${JID2} \
  -o "${LOG_DIR}/03_merge_split.%j.out" -e "${LOG_DIR}/03_merge_split.%j.err" \
  "${SCRIPT_DIR}/../slurm/SLURM_A03_merge_and_split.sh")
echo "  [3/5] Merge+Split:       Job ${JID3} (after ${JID2})"

JID4=$(sbatch --parsable --dependency=afterok:${JID3} \
  -o "${LOG_DIR}/04_annotation.%j.out" -e "${LOG_DIR}/04_annotation.%j.err" \
  "${SCRIPT_DIR}/../slurm/SLURM_A04_annotation_filter.sh")
echo "  [4/5] Annotation:        Job ${JID4} (after ${JID3})"

JID5=$(sbatch --parsable --dependency=afterok:${JID4} \
  -o "${LOG_DIR}/05_summary.%j.out" -e "${LOG_DIR}/05_summary.%j.err" \
  "${SCRIPT_DIR}/../slurm/SLURM_A05_summary_report.sh")
echo "  [5/5] Summary:           Job ${JID5} (after ${JID4})"

# Optional: plots (uncomment when R script is ready)
# JID6=$(sbatch --parsable --dependency=afterok:${JID5} \
#   -o "${LOG_DIR}/06_plots.%j.out" -e "${LOG_DIR}/06_plots.%j.err" \
#   "${SCRIPT_DIR}/../slurm/SLURM_B01_plot_results.sh")
# echo "  [6/6] Plots:             Job ${JID6} (after ${JID5})"

echo ""
echo "=== Pipeline submitted ==="
echo "Chain: ${JID1} → ${JID2} → ${JID3} → ${JID4} → ${JID5}"
echo ""
echo "Monitor:  squeue -u \$(whoami)"
echo "Logs:     tail -f ${LOG_DIR}/*.out"
