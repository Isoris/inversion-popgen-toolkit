#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=80
#SBATCH --mem=237GB
#SBATCH -t 0-12:00:00
#SBATCH -A lt200308
#SBATCH -J m2b_evaladmix
#SBATCH -o logs/evaladmix.%A_%a.out
#SBATCH -e logs/evaladmix.%A_%a.err
# =============================================================================
# SLURM_A03_evaladmix_worker.sh — SLURM array: one evalAdmix run per task
#
# Same task file as STEP_A02. Requires NGSadmix to have completed.
#
# Usage:
#   N=$(tail -n +2 tasks_ngsadmix.tsv | wc -l)
#   sbatch --array=1-${N}%16 slurm/SLURM_A03_evaladmix_worker.sh tasks_ngsadmix.tsv
# =============================================================================
set -euo pipefail
source ~/.bashrc; mamba activate assembly
source "$(dirname "$0")/../00_module2b_config.sh"

TASK_FILE="${1:?Usage: $0 <tasks.tsv>}"
LINE_NUM=$((SLURM_ARRAY_TASK_ID + 1))
LINE=$(sed -n "${LINE_NUM}p" "$TASK_FILE")
[[ -n "$LINE" ]] || m2b_die "No line ${LINE_NUM}"

SCOPE=$(echo "$LINE"   | cut -f1)
THIN=$(echo "$LINE"    | cut -f2)
SSET=$(echo "$LINE"    | cut -f3)
NSAMP=$(echo "$LINE"   | cut -f4)
K=$(echo "$LINE"       | cut -f6)
SEED=$(echo "$LINE"    | cut -f7)
BEAGLE=$(echo "$LINE"  | cut -f8)

P=${SLURM_CPUS_PER_TASK}

RUN_TAG="$(build_run_tag "$SCOPE" "$THIN" "$SSET" "$NSAMP")"
RUN_DIR="${MODULE2B_RESULTS}/${RUN_TAG}"
mkdir -p "${MODULE2B_EVALADMIX}"

STEM="${RUN_TAG}_K$(printf "%02d" "$K")_seed${SEED}"
FOPT="${RUN_DIR}/${STEM}.fopt.gz"
QOPT="${RUN_DIR}/${STEM}.qopt"
COROUT="${MODULE2B_EVALADMIX}/${STEM}.corres.txt"
EVALLOG="${MODULE2B_EVALADMIX}/${STEM}.evaladmix.log"

# Need both NGSadmix outputs
[[ -s "$FOPT" ]] || { m2b_log "SKIP ${STEM} (no fopt)"; exit 0; }
[[ -s "$QOPT" ]] || { m2b_log "SKIP ${STEM} (no qopt)"; exit 0; }

# Skip if already done
[[ -s "$COROUT" ]] && { m2b_log "SKIP ${STEM} (done)"; exit 0; }

# Check binary
[[ -x "$EVALADMIX_BIN" ]] || m2b_die "evalAdmix not found: ${EVALADMIX_BIN}"

m2b_log "RUN evalAdmix ${STEM}"

"$EVALADMIX_BIN" \
  -beagle "$BEAGLE" -fname "$FOPT" -qname "$QOPT" \
  -P "$P" -nIts ${EVALADMIX_NITS} -misTol ${EVALADMIX_MISTOL} \
  -minMaf ${EVALADMIX_MINMAF} \
  -o "$COROUT" \
  > "$EVALLOG" 2>&1

m2b_log "DONE ${STEM}"
