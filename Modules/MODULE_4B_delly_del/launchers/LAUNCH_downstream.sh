#!/usr/bin/env bash
# =============================================================================
# run_downstream_analysis.sh — Submit downstream DEL analysis to SLURM
# =============================================================================
# Run this AFTER the main pipeline (run_delly_pipeline.sh) has completed.
#
# Usage:
#   cd /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/delly_sv
#   bash run_downstream_analysis.sh
# =============================================================================
set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${SCRIPT_DIR}/../00_module4b_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
set -a
source "${CONFIG}"
set +a

LOG_DIR="${OUTDIR}/logs"
mkdir -p "${LOG_DIR}"

echo "=== DELLY DEL Downstream Analysis Launcher ==="
echo "Script directory: ${SCRIPT_DIR}"
echo ""

# Verify upstream outputs exist
for f in "${DIR_FINAL}/catalog_226.DEL.vcf.gz" \
         "${DIR_SUMMARY}/DEL_binary_genotype_matrix.tsv" \
         "${DIR_ANNOT}/catalog_226.functional_class.tsv"; do
  [[ -f "$f" ]] || { echo "ERROR: Missing upstream output: $f" >&2; echo "Run the main pipeline first." >&2; exit 1; }
done

echo "Submitting 2-job dependency chain..."
echo ""

JID7=$(sbatch --parsable \
  -o "${LOG_DIR}/07_downstream.%j.out" -e "${LOG_DIR}/07_downstream.%j.err" \
  "${SCRIPT_DIR}/../slurm/SLURM_B02_downstream_analysis.sh")
echo "  [1/2] Downstream tables:  Job ${JID7}"

JID8=$(sbatch --parsable --dependency=afterok:${JID7} \
  -o "${LOG_DIR}/08_ext_plots.%j.out" -e "${LOG_DIR}/08_ext_plots.%j.err" \
  "${SCRIPT_DIR}/../slurm/SLURM_B03_extended_plots.sh")
echo "  [2/2] Extended plots:     Job ${JID8} (after ${JID7})"

echo ""
echo "=== Downstream analysis submitted ==="
echo "Chain: ${JID7} → ${JID8}"
echo ""
echo "Monitor:  squeue -u \$(whoami)"
echo "Logs:     tail -f ${LOG_DIR}/*.out"
echo ""
echo "After completion, results will be in:"
echo "  Tables: ${DIR_SUMMARY}/"
echo "  Plots:  ${DIR_PLOTS}/"
