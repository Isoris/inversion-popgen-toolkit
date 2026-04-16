#!/usr/bin/env bash
# =============================================================================
# 05_run_all_plots.sh — Run all ROH module plots + statistics + report
# =============================================================================
set -euo pipefail
STEPS="$(cd "$(dirname "$0")" && pwd)"
UTILS="$(cd "$(dirname "$0")/../utils" && pwd)"
source "$(cd "$(dirname "$0")/.." && pwd)/00_module3_config.sh"
hr_init_dirs

hr_log "=== Step 05: Run all plots and statistics ==="

Rscript -e 'for(p in c("data.table","ggplot2")) if(!require(p, character.only=TRUE)) stop(paste("Missing R package:", p))' || {
  hr_die "Missing required R packages: data.table, ggplot2"
}

MASTER="${DIR_TABLES}/master_summary.tsv"
PER_CHR="${DIR_TABLES}/per_chr_roh_summary.tsv"
BINS="${DIR_TABLES}/catfish_roh.per_sample_roh_bins_long.tsv"
HET_TSV="${DIR_HET}/04_summary/genomewide_heterozygosity.tsv"
CHROM_SIZES="${DIR_INPUTS}/chrom_sizes.tsv"

hr_check_file "${MASTER}" "master_summary.tsv"

# ── 1. Core heterozygosity plots ────────────────────────────────────────
hr_log "Running plot_heterozygosity_core.R..."
if [[ -f "${HET_TSV}" ]]; then
  Rscript "${UTILS}/plot_heterozygosity_core.R" \
    "${HET_TSV}" "${DIR_PLOTS_CORE}" \
    ${ANCESTRY_LABELS:+"${ANCESTRY_LABELS}"} \
    2>&1 | tee "${DIR_LOGS}/plot_het_core.log"
fi

# ── 2. Core ROH plots ──────────────────────────────────────────────────
hr_log "Running plot_roh_core.R..."
if [[ -f "${BINS}" ]]; then
  Rscript "${UTILS}/plot_roh_core.R" \
    "${MASTER}" "${BINS}" "${DIR_PLOTS_CORE}" \
    2>&1 | tee "${DIR_LOGS}/plot_roh_core.log"
else
  Rscript "${UTILS}/plot_roh_core.R" \
    "${MASTER}" "/dev/null" "${DIR_PLOTS_CORE}" \
    2>&1 | tee "${DIR_LOGS}/plot_roh_core.log" || true
fi

# ── 3. Scatter plots ───────────────────────────────────────────────────
hr_log "Running plot_scatter_stats.R..."
Rscript "${UTILS}/plot_scatter_stats.R" \
  "${MASTER}" "${DIR_PLOTS_CORE}" \
  "" \
  ${ANCESTRY_LABELS:+"${ANCESTRY_LABELS}"} \
  2>&1 | tee "${DIR_LOGS}/plot_scatter.log"

# ── 4. Chromosome-level plots ─────────────────────────────────────────
if [[ -f "${PER_CHR}" ]]; then
  hr_log "Running plot_roh_by_chromosome.R..."
  Rscript "${UTILS}/plot_roh_by_chromosome.R" \
    "${PER_CHR}" "${DIR_PLOTS_CORE}" \
    ${CHROM_SIZES:+"${CHROM_SIZES}"} \
    ${ANCESTRY_LABELS:+"${ANCESTRY_LABELS}"} \
    2>&1 | tee "${DIR_LOGS}/plot_chr.log"
fi

# ── 5. Theta ideogram plots (main scale) ────────────────────────────
if [[ -d "${DIR_HET}/03_theta" && -f "${CHROM_SIZES}" ]]; then
  hr_log "Running plot_theta_ideogram.R (main 500kb)..."
  Rscript "${UTILS}/plot_theta_ideogram.R" \
    "${DIR_HET}/03_theta" "${CHROM_SIZES}" "${DIR_PLOTS_CORE}" \
    "${SAMPLE_LIST}" 10 \
    2>&1 | tee "${DIR_LOGS}/plot_ideogram_main.log"

  # Multiscale theta ideograms
  if [[ "${RUN_EXTRA_THETA_SCALES:-0}" -eq 1 && -d "${DIR_HET}/03_theta/multiscale" ]]; then
    for i in "${!THETA_SCALE_LABELS[@]}"; do
      LABEL="${THETA_SCALE_LABELS[$i]}"
      hr_log "  Theta ideogram: ${LABEL}..."
      MS_OUT="${DIR_PLOTS_CORE}/theta_${LABEL}"
      mkdir -p "${MS_OUT}"
      Rscript "${UTILS}/plot_theta_ideogram.R" \
        "${DIR_HET}/03_theta/multiscale" "${CHROM_SIZES}" "${MS_OUT}" \
        "${SAMPLE_LIST}" 10 \
        2>&1 | tee "${DIR_LOGS}/plot_ideogram_${LABEL}.log" || true
    done
  fi
fi

# ── 6. Metadata overlay plots ────────────────────────────────────────
hr_log "Running plot_roh_metadata_overlays.R..."
META_ARGS=()
[[ -n "${ANCESTRY_LABELS}" && -f "${ANCESTRY_LABELS}" ]] && META_ARGS+=(--ancestry "${ANCESTRY_LABELS}")
[[ -n "${COVTREE_ORDER}" && -f "${COVTREE_ORDER}" ]] && META_ARGS+=(--order "${COVTREE_ORDER}")
[[ -n "${PRUNED81_SAMPLES}" && -f "${PRUNED81_SAMPLES}" ]] && META_ARGS+=(--pruned81 "${PRUNED81_SAMPLES}")
[[ -n "${NGSRELATE_PAIRS}" && -f "${NGSRELATE_PAIRS}" ]] && META_ARGS+=(--kinship "${NGSRELATE_PAIRS}")

if [[ -f "${PER_CHR}" ]]; then
  Rscript "${UTILS}/plot_roh_metadata_overlays.R" \
    "${MASTER}" "${PER_CHR}" "${BINS:-/dev/null}" "${DIR_PLOTS_META}" \
    "${META_ARGS[@]}" \
    2>&1 | tee "${DIR_LOGS}/plot_metadata.log"
fi

# ── 7. Statistics (sample-level + chromosome-level) ──────────────────
hr_log "Running run_stats.R..."
STATS_ARGS=("${MASTER}" "${DIR_STATS}")
[[ -n "${ANCESTRY_LABELS}" && -f "${ANCESTRY_LABELS}" ]] && STATS_ARGS+=("${ANCESTRY_LABELS}") || STATS_ARGS+=("")
[[ -f "${PER_CHR}" ]] && STATS_ARGS+=("${PER_CHR}") || STATS_ARGS+=("")

Rscript "${UTILS}/run_stats.R" "${STATS_ARGS[@]}" \
  2>&1 | tee "${DIR_LOGS}/run_stats.log"

# ── 8. Report ────────────────────────────────────────────────────────
hr_log "Running 06_write_report.py..."
python3 "${UTILS}/write_report.py" \
  --tables-dir "${DIR_TABLES}" \
  --stats-dir "${DIR_STATS}" \
  --out-dir "${DIR_REPORT}" \
  --ngsf-reps "${NGSFHMM_REPS}" \
  2>&1 | tee "${DIR_LOGS}/write_report.log"

hr_log "=== Step 05 complete ==="
hr_log "Core plots:     ${DIR_PLOTS_CORE}"
hr_log "Metadata plots: ${DIR_PLOTS_META}"
hr_log "Statistics:     ${DIR_STATS}"
hr_log "Report:         ${DIR_REPORT}"
