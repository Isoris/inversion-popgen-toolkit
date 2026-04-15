#!/usr/bin/env bash
# =============================================================================
# LAUNCH_module3.sh — Master runner for het/ROH/FROH module (steps 1-5)
# =============================================================================
# This module produces: het, local theta, ROH, FROH, summary tables, plots, stats.
# Breeding/scoring is a separate downstream module.
#
# Usage:
#   bash launchers/LAUNCH_module3.sh                    # Run everything sequentially
#   bash launchers/LAUNCH_module3.sh --step 2           # Run only step 2
#   bash launchers/LAUNCH_module3.sh --from 3           # Run from step 3 onward
#   bash launchers/LAUNCH_module3.sh --step 2 --slurm   # Submit step 2 as SLURM array
# =============================================================================
set -euo pipefail
STEPS="$(cd "$(dirname "$0")/../steps" && pwd)"
SLURM_DIR="$(cd "$(dirname "$0")/../slurm" && pwd)"
UTILS="$(cd "$(dirname "$0")/../utils" && pwd)"

STEP=""
FROM=1
USE_SLURM=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --step)  STEP="$2"; shift 2 ;;
    --from)  FROM="$2"; shift 2 ;;
    --slurm) USE_SLURM=true; shift ;;
    *) echo "Unknown arg: $1"; exit 1 ;;
  esac
done

run_step() {
  local step_num="$1" desc="$2"; shift 2
  [[ -n "${STEP}" && "${step_num}" != "${STEP}" ]] && return 0
  [[ "${step_num}" -lt "${FROM}" ]] && return 0
  echo ""
  echo "================================================================"
  echo "  STEP ${step_num}: ${desc}"
  echo "================================================================"
  "$@"
}

# ── Step 1: Prepare inputs ──────────────────────────────────────────────
run_step 1 "Prepare and validate inputs" \
  bash "${STEPS}/STEP_A01_prep_inputs.sh"

# ── Step 2: Per-sample heterozygosity ───────────────────────────────────
if [[ "${USE_SLURM}" == true && ("${STEP}" == "2" || -z "${STEP}") ]]; then
  echo ""
  echo "================================================================"
  echo "  STEP 2: Per-sample heterozygosity (SLURM array)"
  echo "================================================================"
  source "$(dirname "$0")/../00_module3_config.sh"
  N=$(wc -l < "${SAMPLE_LIST}")
  echo "Submitting SLURM array for ${N} samples..."
  sbatch --array=1-${N} "${SLURM_DIR}/SLURM_A02_heterozygosity_worker.sh"
  echo "Submitted. Wait for completion before running step 3."
  if [[ -z "${STEP}" ]]; then
    echo "Stopping sequential run. Re-run with --from 3 after SLURM completes."
    exit 0
  fi
else
  run_step 2 "Per-sample heterozygosity (sequential)" \
    bash "${STEPS}/STEP_A02_run_heterozygosity.sh"
fi

# ── Step 3: ngsF-HMM ──────────────────────────────────────────────────
run_step 3 "ngsF-HMM (multi-replicate)" \
  bash "${STEPS}/STEP_A03_run_ngsF_HMM.sh"

# ── Step 4: Parse ROH + het in/out ROH ────────────────────────────────
run_step 4 "Parse ROH, compute FROH, het in/out ROH" \
  bash "${STEPS}/STEP_A04_parse_roh_and_het.sh"

# ── Step 5: Plots + stats + report ───────────────────────────────────
run_step 5 "Generate all plots, statistics, and report" \
  bash "${STEPS}/STEP_B01_run_all_plots.sh"

echo ""
echo "================================================================"
echo "  ROH MODULE COMPLETE"
echo "  Breeding/scoring module can be run separately if needed."
echo "================================================================"
