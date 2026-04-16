#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=64G
#SBATCH -t 0-04:00:00
#SBATCH -J bnd_plot
#SBATCH -o logs/06_plots.%j.out
#SBATCH -e logs/06_plots.%j.err
set -euo pipefail
source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${SCRIPT_DIR}/../00_module4e_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
set -a
source "${CONFIG}"
set +a

RSCRIPT_BIN="/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript"
"${RSCRIPT_BIN}" -e 'suppressPackageStartupMessages({library(cowplot); library(ComplexHeatmap)}); cat("R deps OK\n")'

dv_init_dirs
dv_log "=== STEP 6: Publication-style BND plots ==="

for f in \
  "${DIR_SUMMARY}/per_sample_BND_counts.tsv" \
  "${DIR_SUMMARY}/per_chromosome_BND_counts.tsv" \
  "${DIR_SUMMARY}/pairwise_shared_BND.tsv" \
  "${DIR_SUMMARY}/BND_binary_genotype_matrix.tsv" \
  "${DIR_FINAL}/catalog_226.BND.GT_matrix.tsv" \
  "${DIR_SUMMARY}/BND_window_counts_1Mb.tsv"; do
  dv_check_file "$f" "$(basename "$f")"
done

"${RSCRIPT_BIN}" "${SCRIPT_DIR}/../utils/plot_BND_results.R" \
  --summary_dir "${DIR_SUMMARY}" \
  --final_dir   "${DIR_FINAL}" \
  --plot_dir    "${DIR_PLOTS}" \
  --ref_fai     "${REF_FAI}" \
  --samples_unrelated "${SAMPLES_UNRELATED}"

dv_log "=== PLOTS COMPLETE ==="
dv_log "Output: ${DIR_PLOTS}/"
ls -lh "${DIR_PLOTS}"/*.pdf "${DIR_PLOTS}"/*.png 2>/dev/null | while read -r line; do
  dv_log "  ${line}"
done
