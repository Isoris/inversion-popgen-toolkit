#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=64G
#SBATCH -t 0-04:00:00
#SBATCH -J manta_plot
#SBATCH -o logs/06_plots.%j.out
#SBATCH -e logs/06_plots.%j.err
# =============================================================================
# 06_plot_results.sh — v3: tier composition + evidence panels
# =============================================================================
#
# Panels (generated for both 226 and 81 cohorts where applicable):
#   A. Tier composition waterfall (PASS → excluded → tier1 → tier2 → tier3)
#   B. Genome-wide DEL density heatmap (1-Mb windows)
#   C. Per-sample SV burden — PASS level (stacked by type)
#   D. Per-sample SV burden — tier2 level (stacked by type)
#   E. Size distribution per SV type (violin + boxplot)
#   F. Sharing spectrum per type
#   G. Per-chromosome SV counts
#   H. Evidence distribution (sum PR+SR alt) per type with tier thresholds
#   I. Precision composition (PRECISE vs IMPRECISE) per type
#   J. Pairwise shared DEL heatmap (226)
#
# Prerequisites: run filter_all_sv_types.sh AND 05_summary_report.sh first.
# =============================================================================
set -euo pipefail
source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${SCRIPT_DIR}/../00_module4g_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
set -a
source "${CONFIG}"
set +a

FILT_DIR="${OUTDIR}/10_filtered_catalogs"
RSCRIPT_BIN="/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript"

mv_log "=== STEP 6: Publication-style SV plots (v3) ==="

# Check R deps
mv_log "Checking R dependencies..."
"${RSCRIPT_BIN}" -e '
  pkgs <- c("data.table","ggplot2","cowplot","viridis","scales","optparse")
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]
  if(length(missing) > 0) {
    cat("Missing R packages:", paste(missing, collapse=", "), "\n")
    quit(status=1)
  }
  cat("All R dependencies OK\n")
'

mv_init_dirs

# Verify key inputs
for f in \
  "${DIR_SUMMARY}/tier_composition_226.tsv" \
  "${DIR_SUMMARY}/tier_composition_81.tsv" \
  "${DIR_SUMMARY}/master_counts.tsv"; do
  if [[ ! -f "$f" ]]; then
    mv_die "Missing: $(basename "$f") — run 05_summary_report.sh first"
  fi
done

if [[ ! -d "${FILT_DIR}" ]]; then
  mv_die "Missing: ${FILT_DIR} — run filter_all_sv_types.sh first"
fi

R_SCRIPT="${SCRIPT_DIR}/../utils/plot_manta_results.R"
mv_check_file "${R_SCRIPT}" "R plotting script"

mv_log "Running R plots..."
"${RSCRIPT_BIN}" "${R_SCRIPT}" \
  --summary_dir "${DIR_SUMMARY}" \
  --final_dir   "${DIR_FINAL}" \
  --filter_dir  "${FILT_DIR}" \
  --plot_dir    "${DIR_PLOTS}" \
  --ref_fai     "${REF_FAI}"

mv_log "=== PLOTS COMPLETE ==="
mv_log "Output: ${DIR_PLOTS}/"
ls -lh "${DIR_PLOTS}"/*.pdf "${DIR_PLOTS}"/*.png 2>/dev/null | while read line; do mv_log "  ${line}"; done
