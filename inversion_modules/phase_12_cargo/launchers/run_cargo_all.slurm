#!/usr/bin/env bash
#SBATCH --job-name=cargo_array
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/cargo_%A_%a.out
#SBATCH --error=logs/cargo_%A_%a.err
# SBATCH --array=1-N    # set N via --array on submission, derived from diagnostic_table.tsv
set -euo pipefail

# =============================================================================
# run_cargo_all.slurm — SLURM array driver for MODULE_6_Cargo
#
# Each array task processes one candidate:
#   1. STEP_C61 per-arrangement burden + NS/SS + LoF + FST + SFS  (Level 1)
#   2. STEP_C62 configuration spectrum                              (Level 2, if eligible)
#
# Submission pattern (after STEP_C60 has produced diagnostic_table.tsv):
#   N=$(($(wc -l < ${CARGO_INVENTORY_DIR}/diagnostic_table.tsv) - 1))
#   sbatch --array=1-${N} run_cargo_all.slurm
#
# Or to process only Level 1 candidates:
#   N=$(awk -F'\t' 'NR>1 && $13=="1"' ${CARGO_INVENTORY_DIR}/diagnostic_table.tsv | wc -l)
#   sbatch --array=1-${N} run_cargo_all.slurm
# =============================================================================

mkdir -p logs

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../00_cargo_config.sh"

DIAG="${CARGO_INVENTORY_DIR}/diagnostic_table.tsv"
if [[ ! -f "${DIAG}" ]]; then
  echo "[FATAL] Diagnostic table not found: ${DIAG}"
  echo "  Run STEP_C60_cargo_inventory.py first"
  exit 1
fi

# Build candidate list — all rows in diagnostic table, regardless of level.
# Level gate is applied INSIDE the python scripts, so failing-gate candidates
# exit cleanly.
CANDS=( $(awk -F'\t' 'NR>1 {print $1}' "${DIAG}") )

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
IDX=$((TASK_ID - 1))

if [[ ${IDX} -ge ${#CANDS[@]} ]]; then
  echo "[ERROR] task ${TASK_ID} > ${#CANDS[@]} candidates"
  exit 2
fi

CID="${CANDS[$IDX]}"
echo "==================================================================="
echo "[$(date)] Cargo array task ${TASK_ID}/${#CANDS[@]} → candidate ${CID}"
echo "==================================================================="

# Step C61: per-arrangement burden, NS/SS, LoF, FST, SFS
echo ""
echo "[$(date)] === STEP_C61 (Level 1) ==="
python3 "${SCRIPT_DIR}/../compute/STEP_C61_per_arrangement_burden.py" \
    --candidate-id "${CID}" \
    || { echo "[ERROR] STEP_C61 failed for ${CID}"; exit 11; }

# Step C62: configuration spectrum (skipped internally if level_2_ok != 1)
echo ""
echo "[$(date)] === STEP_C62 (Level 2) ==="
python3 "${SCRIPT_DIR}/../compute/STEP_C62_config_spectrum.py" \
    --candidate-id "${CID}" \
    || { echo "[WARN] STEP_C62 failed or skipped for ${CID}"; }

echo ""
echo "[$(date)] === Cargo task ${TASK_ID} complete: ${CID} ==="
