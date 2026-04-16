#!/usr/bin/env bash
# =============================================================================
# LAUNCH_module2b_compute.sh — Submit all compute jobs with dependency chain
#
# Order:
#   1. STEP_A01: generate task list
#   2. STEP_A02: NGSadmix array (sbatch)
#   3. STEP_A03: evalAdmix array (sbatch --dependency=afterok:A02)
#
# Usage:
#   bash launchers/LAUNCH_module2b_compute.sh [--per-chr] [--dry-run]
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../00_module2b_config.sh"

SCRIPT_DIR="$(cd "$(dirname "$0")/../scripts" && pwd)"
SLURM_SCRIPTS="$(cd "$(dirname "$0")/../slurm" && pwd)"
MODULE_ROOT="$(cd "$(dirname "$0")/.." && pwd)"

PER_CHR=""
DRY_RUN=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --per-chr) PER_CHR="--per-chr"; shift;;
    --dry-run) DRY_RUN=1; shift;;
    *) echo "Usage: $0 [--per-chr] [--dry-run]"; exit 1;;
  esac
done

m2b_init_dirs
mkdir -p "${MODULE_ROOT}/logs"

# ── Step 1: Generate task list ───────────────────────────────────────────────
TASK_FILE="${MODULE_ROOT}/tasks_ngsadmix.tsv"
m2b_log "Generating task list..."
bash "${SCRIPT_DIR}/STEP_A01_generate_tasks.sh" ${PER_CHR} > "$TASK_FILE"

N_TASKS=$(tail -n +2 "$TASK_FILE" | wc -l)
m2b_log "Tasks generated: ${N_TASKS} (in ${TASK_FILE})"

if [[ "$N_TASKS" -eq 0 ]]; then
  m2b_die "No tasks generated — check BEAGLE files exist"
fi

# ── Step 2: Submit NGSadmix array ────────────────────────────────────────────
m2b_log "Submitting NGSadmix array (${N_TASKS} tasks, max 32 concurrent)..."

if [[ "$DRY_RUN" -eq 1 ]]; then
  echo "[DRY-RUN] sbatch --array=1-${N_TASKS}%32 ${SLURM_SCRIPTS}/SLURM_A02_ngsadmix_worker.sh ${TASK_FILE}"
  JOB_A02="DRYRUN"
else
  JOB_A02=$(sbatch --parsable \
    --array=1-${N_TASKS}%32 \
    "${SLURM_SCRIPTS}/SLURM_A02_ngsadmix_worker.sh" "${TASK_FILE}")
  m2b_log "NGSadmix job: ${JOB_A02}"
fi

# ── Step 3: Submit evalAdmix array (depends on NGSadmix) ────────────────────
m2b_log "Submitting evalAdmix array (depends on ${JOB_A02})..."

if [[ "$DRY_RUN" -eq 1 ]]; then
  echo "[DRY-RUN] sbatch --dependency=afterany:${JOB_A02} --array=1-${N_TASKS}%16 ${SLURM_SCRIPTS}/SLURM_A03_evaladmix_worker.sh ${TASK_FILE}"
else
  JOB_A03=$(sbatch --parsable \
    --dependency=afterany:${JOB_A02} \
    --array=1-${N_TASKS}%16 \
    "${SLURM_SCRIPTS}/SLURM_A03_evaladmix_worker.sh" "${TASK_FILE}")
  m2b_log "evalAdmix job: ${JOB_A03} (depends on ${JOB_A02})"
fi

# ── Summary ──────────────────────────────────────────────────────────────────
m2b_log "=== Compute submission complete ==="
m2b_log "Task file: ${TASK_FILE}"
m2b_log "Total tasks: ${N_TASKS}"
m2b_log ""
m2b_log "After jobs complete, run:"
m2b_log "  bash launchers/LAUNCH_module2b_figures.sh"
