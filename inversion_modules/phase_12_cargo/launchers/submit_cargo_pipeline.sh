#!/usr/bin/env bash
# =============================================================================
# submit_cargo_pipeline.sh — submit the full MODULE_6_Cargo pipeline as a
# chain of SLURM jobs with afterok dependencies.
#
# Order:
#   1. STEP_60 build_gene_intervals     (one-time, idempotent)
#   2. STEP_61 eggnog_annotate          (one-time, optional)
#   3. STEP_C60 inventory + diagnostic gate
#   4. run_cargo_all.slurm  (array, one task per candidate; runs C61 + C62)
#   5. STEP_C63 functional enrichment
#   6. STEP_C64 cross-inversion synthesis
#   7. STEP_C65 per-inversion figures
#
# Usage:
#   bash submit_cargo_pipeline.sh
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../00_cargo_config.sh"

cd "${SCRIPT_DIR}/.."
mkdir -p logs

submit_dep() {
  local name="$1"; shift
  local dep="$1"; shift
  local script="$1"; shift
  if [[ -n "${dep}" ]]; then
    sbatch --parsable --dependency=afterok:"${dep}" "${script}" "$@"
  else
    sbatch --parsable "${script}" "$@"
  fi
}

# ── Step 60 ──
JID_60=""
if [[ ! -f "${GENE_BED}" || ! -f "${TRANSCRIPT_GENE_MAP}" ]]; then
  JID_60=$(sbatch --parsable STEP_60_build_gene_intervals.sh)
  echo "Submitted STEP_60: ${JID_60}"
else
  echo "STEP_60 outputs exist — skipping"
fi

# ── Step 61 (optional) ──
JID_61="${JID_60}"
if [[ ! -f "${GENE_FUNCTION_TSV}" ]]; then
  JID_61=$(submit_dep "step61" "${JID_60}" STEP_61_eggnog_annotate.sh)
  echo "Submitted STEP_61: ${JID_61}"
else
  echo "STEP_61 output exists — skipping"
fi

# ── Step C60 ──
JID_C60=$(submit_dep "C60" "${JID_61}" compute/STEP_C60_cargo_inventory.slurm)
echo "Submitted STEP_C60: ${JID_C60}"

# ── Cargo array (C61 + C62) ──
# We can't compute N until the diagnostic table exists, so submit the array
# with a placeholder and resize on the fly via a small shim.
#
# Approach: submit a single-task waker that, once C60 finishes, queues the
# real array job using sbatch from inside the waker.

WAKER=$(mktemp -p logs cargo_waker.XXXXXX.sh)
cat > "${WAKER}" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=cargo_waker
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --time=00:05:00
#SBATCH --mem=512M
#SBATCH --output=logs/waker_%j.out
#SBATCH --error=logs/waker_%j.err
#SBATCH --dependency=afterok:${JID_C60}
set -euo pipefail
source "${SCRIPT_DIR}/../00_cargo_config.sh"
N=\$(awk -F'\\t' 'NR>1 {n++} END{print n+0}' "\${CARGO_INVENTORY_DIR}/diagnostic_table.tsv")
echo "[waker] queueing cargo array with N=\$N tasks"
sbatch --array=1-\${N} "${SCRIPT_DIR}/../launchers/run_cargo_all.slurm"
EOF
chmod +x "${WAKER}"
JID_WAKER=$(sbatch --parsable "${WAKER}")
echo "Submitted cargo waker: ${JID_WAKER} (will queue array after diagnostic table is built)"

# ── Note ──
# The downstream analysis steps (C63/C64/C65) cannot use afterok on the array
# from this script because the array is queued by the waker, not directly by
# us. Submit them manually after the array completes, or use a second waker:
echo ""
echo "After array completes, submit the analysis layer:"
echo "  sbatch ${SCRIPT_DIR}/../analysis/STEP_C63_functional_enrichment.slurm"
echo "  sbatch ${SCRIPT_DIR}/../analysis/STEP_C64_cross_inversion_synthesis.slurm"
echo "  sbatch ${SCRIPT_DIR}/../plot/STEP_C65_per_inversion_profile.slurm"
