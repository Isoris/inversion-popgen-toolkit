#!/usr/bin/env bash
# =============================================================================
# STEP_A01_generate_tasks.sh — Generate NGSadmix task list for array submission
#
# Output: TSV to stdout with columns:
#   scope  thin  sample_set  n_samples  samples_file  K  seed  beagle_path
#
# Produces tasks for:
#   - Whole-genome × thin500 × all 226
#   - Whole-genome × thin1000 × all 226
#   - Whole-genome × thin500 × pruned 81
#   - Per-chromosome (LG01..LG28) × thin500 × all 226  (optional: --per-chr)
#
# Usage:
#   bash scripts/STEP_A01_generate_tasks.sh > tasks_ngsadmix.tsv
#   bash scripts/STEP_A01_generate_tasks.sh --per-chr > tasks_ngsadmix_with_chr.tsv
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../00_module2b_config.sh"

PER_CHR=0
[[ "${1:-}" == "--per-chr" ]] && PER_CHR=1

# Header
echo -e "scope\tthin\tsample_set\tn_samples\tsamples_file\tK\tseed\tbeagle_path"

for THIN in "${THIN_PANELS[@]}"; do
  # Whole-genome BEAGLE
  WG_BEAGLE="${BEAGLE_MERGED_DIR}/catfish.wholegenome.byRF.thin_${THIN}.beagle.gz"

  for K in $(seq "$K_MIN" "$K_MAX"); do
    for SEED in "${SEEDS[@]}"; do

      # Batch 1/2: all 226 × thin500/thin1000
      if [[ -s "$WG_BEAGLE" ]]; then
        echo -e "wholegenome\t${THIN}\tall\t${N_SAMPLES}\t${SAMPLE_LIST}\t${K}\t${SEED}\t${WG_BEAGLE}"
      fi

      # Batch 3: pruned 81 × thin500 only
      if [[ "$THIN" == "500" && -s "$PRUNED_LIST" ]]; then
        echo -e "wholegenome\t${THIN}\tpruned\t${N_PRUNED}\t${PRUNED_LIST}\t${K}\t${SEED}\t${WG_BEAGLE}"
      fi

    done
  done
done

# Per-chromosome tasks (thin500 only, all 226 only)
if [[ "$PER_CHR" -eq 1 ]]; then
  THIN=500
  for CHR_NUM in $(seq 1 "$N_CHROMOSOMES"); do
    LG=$(printf "LG%02d" "$CHR_NUM")
    CHR_BEAGLE="${BEAGLE_BYCHR_DIR}/thin_${THIN}/catfish.${LG}.byRF.thin_${THIN}.beagle.gz"
    [[ -s "$CHR_BEAGLE" ]] || continue

    for K in $(seq "$K_MIN" "$K_MAX"); do
      for SEED in "${SEEDS[@]}"; do
        echo -e "${LG}\t${THIN}\tall\t${N_SAMPLES}\t${SAMPLE_LIST}\t${K}\t${SEED}\t${CHR_BEAGLE}"
      done
    done
  done
fi
