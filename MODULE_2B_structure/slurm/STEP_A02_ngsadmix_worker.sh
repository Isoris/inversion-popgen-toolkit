#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=80
#SBATCH --mem=237GB
#SBATCH -t 1-00:00:00
#SBATCH -A lt200308
#SBATCH -J m2b_ngsadmix
#SBATCH -o logs/ngsadmix.%A_%a.out
#SBATCH -e logs/ngsadmix.%A_%a.err
# =============================================================================
# STEP_A02_ngsadmix_worker.sh — SLURM array: one NGSadmix run per task
#
# Reads task line from STEP_A01 output (skip header).
# Columns: scope, thin, sample_set, n_samples, samples_file, K, seed, beagle_path
#
# Usage:
#   N=$(tail -n +2 tasks_ngsadmix.tsv | wc -l)
#   sbatch --array=1-${N}%32 slurm/STEP_A02_ngsadmix_worker.sh tasks_ngsadmix.tsv
# =============================================================================
set -euo pipefail
source ~/.bashrc; mamba activate assembly
source "$(dirname "$0")/../00_module2b_config.sh"

TASK_FILE="${1:?Usage: $0 <tasks.tsv>}"
m2b_check_file "$TASK_FILE"

# +1 to skip header
LINE_NUM=$((SLURM_ARRAY_TASK_ID + 1))
LINE=$(sed -n "${LINE_NUM}p" "$TASK_FILE")
[[ -n "$LINE" ]] || m2b_die "No line ${LINE_NUM} in ${TASK_FILE}"

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
mkdir -p "${RUN_DIR}"

m2b_check_file "$BEAGLE"

STEM="${RUN_TAG}_K$(printf "%02d" "$K")_seed${SEED}"
PREFIX="${RUN_DIR}/${STEM}"

if [[ -s "${PREFIX}.qopt" ]]; then
  m2b_log "SKIP ${STEM} (exists)"
  exit 0
fi

m2b_log "RUN ${STEM}  K=${K} seed=${SEED}"

NGSadmix -likes "$BEAGLE" -K "$K" -minMaf ${NGSADMIX_MINMAF} \
  -P "$P" -seed "$SEED" -o "$PREFIX" \
  > "${PREFIX}.stdout.log" 2>&1

m2b_log "DONE ${STEM}"
