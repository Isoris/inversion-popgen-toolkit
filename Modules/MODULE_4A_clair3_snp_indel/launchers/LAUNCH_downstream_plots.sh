#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=64G
#SBATCH -t 0-06:00:00
#SBATCH -A lt200308
#SBATCH -J c3_plots
#SBATCH -o logs/c3_plots.%j.out
#SBATCH -e logs/c3_plots.%j.err
# ============================================================
# run_plots.sh — Run all Clair3 downstream plotting scripts
#
# Usage:
#   bash run_plots.sh --chrom C_gar_LG01
#   sbatch run_plots.sh --chrom C_gar_LG01
# ============================================================
set -euo pipefail
source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${SCRIPT_DIR}/00_config.sh"

CHROM=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --chrom) CHROM="$2"; shift 2 ;;
        *) echo "Unknown: $1" >&2; exit 1 ;;
    esac
done
[[ -n "$CHROM" ]] || { echo "[ERROR] Specify --chrom" >&2; exit 1; }

ds_init_dirs "$CHROM"
ds_resolve_ancestry

RSCRIPT_BIN="/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript"
"${RSCRIPT_BIN}" -e 'library(cowplot); library(ComplexHeatmap); cat("R packages OK\n")'

DS_DIR="${DS_ROOT}/${CHROM}"

# Check that data scripts have run
for f in "${DS_DIR}/catalogs/COMBINED_catalog.tsv" \
         "${DS_DIR}/phase_prep/phase_block_catalog.tsv" \
         "${DS_DIR}/per_sample/per_sample_COMBINED_summary.tsv"; do
    [[ -f "$f" ]] || { echo "[ERROR] Missing: $f — run run_downstream.sh first" >&2; exit 1; }
done

ANCESTRY_ARGS=""
[[ -n "${QOPT:-}" ]] && ANCESTRY_ARGS="--qopt ${QOPT} --qopt_samples ${QOPT_SAMPLES}"

ds_log "=== 10: Main publication plots ==="
"${RSCRIPT_BIN}" "${SCRIPT_DIR}/10_plot_main.R" \
    --ds_chrom_dir "${DS_DIR}" \
    --chrom "${CHROM}" \
    --ref_fai "${REF_FAI}" \
    --samples_unrelated "${SAMPLES_UNRELATED}"

ds_log "=== 11: Extended plots ==="
"${RSCRIPT_BIN}" "${SCRIPT_DIR}/11_plot_extended.R" \
    --ds_chrom_dir "${DS_DIR}" \
    --chrom "${CHROM}" \
    --ref_fai "${REF_FAI}" \
    --samples_unrelated "${SAMPLES_UNRELATED}" \
    ${ANCESTRY_ARGS}

ds_log "=== 12: Phase signature plots ==="
"${RSCRIPT_BIN}" "${SCRIPT_DIR}/12_plot_phase_signatures.R" \
    --ds_chrom_dir "${DS_DIR}" \
    --chrom "${CHROM}" \
    --ref_fai "${REF_FAI}" \
    --samples_unrelated "${SAMPLES_UNRELATED}" \
    ${ANCESTRY_ARGS}

ds_log "=== All plots complete ==="
ds_log "Output: ${DS_DIR}/figures/"
ls -lh "${DS_DIR}/figures/"*.pdf "${DS_DIR}/figures/"*.png 2>/dev/null | wc -l | xargs -I{} echo "  Total figure files: {}"
