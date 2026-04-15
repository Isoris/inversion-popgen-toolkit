#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=64G
#SBATCH -t 0-04:00:00
#SBATCH -J delly_ext_plots
#SBATCH -o logs/08_extended_plots.%j.out
#SBATCH -e logs/08_extended_plots.%j.err
# =============================================================================
# 08_extended_plots.sh — Run extended DEL plotting suite
# Requires: 07_downstream_analysis.sh completed (master annotation + summaries)
# =============================================================================
set -euo pipefail
source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${SCRIPT_DIR}/../00_module4b_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
set -a
source "${CONFIG}"
set +a

RSCRIPT_BIN="/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript"
"${RSCRIPT_BIN}" -e 'library(cowplot); cat("cowplot OK\n")'

dv_init_dirs
dv_log "=== EXTENDED DEL PLOTS ==="

# Check that downstream analysis has run
for f in "${DIR_SUMMARY}/master_DEL_annotation_226.tsv" \
         "${DIR_SUMMARY}/per_sample_DEL_summary.tsv"; do
  dv_check_file "$f" "$(basename "$f") — run 07_downstream_analysis.sh first"
done

# Resolve ancestry
ANCESTRY_ARGS=""
if [[ -f "${PA_BEST_SEED_TABLE}" ]]; then
  BEST_K=$(awk -F'\t' 'NR>1 { if ($5 > max || NR==2) { max=$5; k=$1 } } END { print k }' \
           "${PA_BEST_SEED_TABLE}")
  BEST_PREFIX=$(awk -F'\t' -v k="${BEST_K}" 'NR>1 && $1==k { print $4 }' "${PA_BEST_SEED_TABLE}")
  CANDIDATE_QOPT="${PA_NGSADMIX_DIR}/${BEST_PREFIX}.qopt"
  if [[ -f "${CANDIDATE_QOPT}" ]]; then
    ANCESTRY_ARGS="--qopt ${CANDIDATE_QOPT} --qopt_samples ${SAMPLES_ALL}"
    dv_log "  Using ancestry: K=${BEST_K}"
  fi
fi

"${RSCRIPT_BIN}" "${SCRIPT_DIR}/11_plot_DEL_extended.R" \
  --summary_dir "${DIR_SUMMARY}" \
  --final_dir   "${DIR_FINAL}" \
  --plot_dir    "${DIR_PLOTS}" \
  --ref_fai     "${REF_FAI}" \
  --samples_unrelated "${SAMPLES_UNRELATED}" \
  ${ANCESTRY_ARGS}

dv_log "=== EXTENDED PLOTS COMPLETE ==="
dv_log "Output: ${DIR_PLOTS}/"
ls -lh "${DIR_PLOTS}"/*.pdf "${DIR_PLOTS}"/*.png 2>/dev/null | while read line; do dv_log "  ${line}"; done
