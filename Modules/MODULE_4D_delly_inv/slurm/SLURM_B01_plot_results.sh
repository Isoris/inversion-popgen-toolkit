#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=64G
#SBATCH -t 0-04:00:00
#SBATCH -J inv_plot
#SBATCH -o logs/06_plots.%j.out
#SBATCH -e logs/06_plots.%j.err
set -euo pipefail

source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${SCRIPT_DIR}/../00_module4d_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }

set -a
source "${CONFIG}"
set +a

RSCRIPT_BIN="/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript"
"${RSCRIPT_BIN}" -e 'library(cowplot); library(ComplexHeatmap); cat("R deps OK\n")'

dv_init_dirs
dv_log "=== STEP 6: Publication-style INV plots ==="

# ------------------------------------------------------------------
# Core required files
# ------------------------------------------------------------------
GT_MATRIX="${DIR_FINAL}/catalog_226.INV.GT_matrix.tsv"
PER_SAMPLE="${DIR_SUMMARY}/per_sample_INV_counts.tsv"
PER_CHR="${DIR_SUMMARY}/per_chromosome_INV_counts.tsv"
PAIRWISE="${DIR_SUMMARY}/pairwise_shared_INV.tsv"
BINARY="${DIR_SUMMARY}/INV_binary_genotype_matrix.tsv"
SVLEN="${DIR_SUMMARY}/INV_svlen_distribution.tsv"
WINDOWS="${DIR_SUMMARY}/INV_window_counts_1Mb.tsv"

dv_check_file "${GT_MATRIX}" "catalog_226.INV.GT_matrix.tsv"
dv_check_file "${PER_CHR}"   "per_chromosome_INV_counts.tsv"
dv_check_file "${PAIRWISE}"  "pairwise_shared_INV.tsv"
dv_check_file "${BINARY}"    "INV_binary_genotype_matrix.tsv"
dv_check_file "${SVLEN}"     "INV_svlen_distribution.tsv"
dv_check_file "${WINDOWS}"   "INV_window_counts_1Mb.tsv"

# ------------------------------------------------------------------
# Rebuild per-sample counts if missing
# ------------------------------------------------------------------
if [[ ! -s "${PER_SAMPLE}" ]]; then
  dv_log "[WARN] Missing per_sample_INV_counts.tsv — rebuilding from GT matrix"

  python3 - << PYEOF
import csv

gt_matrix = r"${GT_MATRIX}"
per_out   = r"${PER_SAMPLE}"

with open(gt_matrix, newline="") as f:
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)
    samples = header[5:]
    counts = [0] * len(samples)

    alt_like = {"0/1","1/1","0|1","1|0","1|1","1/2","2/1","1|2","2|1","2/2"}

    for row in reader:
        gts = row[5:]
        for i, gt in enumerate(gts):
            if gt in alt_like:
                counts[i] += 1

with open(per_out, "w", newline="") as o:
    w = csv.writer(o, delimiter="\t")
    w.writerow(["sample", "n_INVs"])
    for s, c in zip(samples, counts):
        w.writerow([s, c])

print(f"[OK] rebuilt {per_out}")
PYEOF
fi

dv_check_file "${PER_SAMPLE}" "per_sample_INV_counts.tsv"

# ------------------------------------------------------------------
# Run plotting
# ------------------------------------------------------------------
"${RSCRIPT_BIN}" "${SCRIPT_DIR}/../utils/plot_INV_results.R" \
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
