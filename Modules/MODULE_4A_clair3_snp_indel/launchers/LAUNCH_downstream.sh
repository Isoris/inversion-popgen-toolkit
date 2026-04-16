#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=64G
#SBATCH -t 1-00:00:00
#SBATCH -A lt200308
#SBATCH -J c3_downstream
#SBATCH -o logs/c3_downstream.%j.out
#SBATCH -e logs/c3_downstream.%j.err
# ============================================================
# run_downstream.sh — Run full Clair3 downstream analysis suite
#
# Usage:
#   bash run_downstream.sh --chrom C_gar_LG01
#   sbatch run_downstream.sh --chrom C_gar_LG01
# ============================================================
set -euo pipefail
source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${SCRIPT_DIR}/00_config.sh"

CHROM=""
SKIP_PLOTS=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --chrom)      CHROM="$2"; shift 2 ;;
        --skip_plots) SKIP_PLOTS=1; shift ;;
        *) echo "Unknown: $1" >&2; exit 1 ;;
    esac
done

[[ -n "$CHROM" ]] || { echo "[ERROR] Specify --chrom" >&2; exit 1; }

ds_init_dirs "$CHROM"
ds_resolve_ancestry

PP_CHROM_DIR="${PP_RESULTS}/${CHROM}"
[[ -d "$PP_CHROM_DIR" ]] || { echo "[ERROR] No postprocess results: $PP_CHROM_DIR" >&2; exit 1; }

ds_log "════════════════════════════════════════════════════════"
ds_log " Clair3 Downstream Analysis: ${CHROM}"
ds_log "════════════════════════════════════════════════════════"

# ── 01: Build variant catalogs + GT matrices ──
ds_log "=== 01: Build variant catalogs ==="
python3 "${SCRIPT_DIR}/01_build_variant_catalog.py" \
    --pp_results "${PP_CHROM_DIR}" \
    --chrom "${CHROM}" \
    --outdir "${DS_CHROM_DIR}" \
    --samples_unrelated "${SAMPLES_UNRELATED}"

# ── 02: Phase block catalog + inversion prep ──
ds_log "=== 02: Phase block catalog ==="
python3 "${SCRIPT_DIR}/02_phase_block_catalog.py" \
    --pp_results "${PP_CHROM_DIR}" \
    --chrom "${CHROM}" \
    --outdir "${DS_CHROM_DIR}"

# ── 03: Master annotation ──
ds_log "=== 03: Master annotation ==="
python3 "${SCRIPT_DIR}/03_master_annotation.py" \
    --catalogs_dir "${DS_CHROM_DIR}/catalogs" \
    --phase_prep_dir "${DS_CHROM_DIR}/phase_prep" \
    --chrom "${CHROM}" \
    --outdir "${DS_CHROM_DIR}" \
    --ref_fai "${REF_FAI}"

# ── 04: Per-sample summary ──
ds_log "=== 04: Per-sample summary ==="
ANCESTRY_ARGS=""
[[ -n "${QOPT}" ]] && ANCESTRY_ARGS="--qopt ${QOPT} --qopt_samples ${QOPT_SAMPLES}"
python3 "${SCRIPT_DIR}/04_per_sample_summary.py" \
    --matrices_dir "${DS_CHROM_DIR}/matrices" \
    --annotation_dir "${DS_CHROM_DIR}/annotation" \
    --phase_prep_dir "${DS_CHROM_DIR}/phase_prep" \
    --outdir "${DS_CHROM_DIR}" \
    --samples_unrelated "${SAMPLES_UNRELATED}" \
    ${ANCESTRY_ARGS}

# ── 05: Distance matrices + rare network ──
ds_log "=== 05: Distance matrices ==="
python3 "${SCRIPT_DIR}/05_sample_distance_matrices.py" \
    --matrices_dir "${DS_CHROM_DIR}/matrices" \
    --catalogs_dir "${DS_CHROM_DIR}/catalogs" \
    --phase_prep_dir "${DS_CHROM_DIR}/phase_prep" \
    --outdir "${DS_CHROM_DIR}" \
    --chrom "${CHROM}" \
    --ref_fai "${REF_FAI}"

# ── 06: Marker selection ──
ds_log "=== 06: Marker selection ==="
python3 "${SCRIPT_DIR}/06_marker_selection.py" \
    --annotation_dir "${DS_CHROM_DIR}/annotation" \
    --matrices_dir "${DS_CHROM_DIR}/matrices" \
    --outdir "${DS_CHROM_DIR}" \
    --per_sample_dir "${DS_CHROM_DIR}/per_sample"

# ── 07: Gene/chromosome tables ──
ds_log "=== 07: Gene/chromosome tables ==="
python3 "${SCRIPT_DIR}/07_gene_chr_tables.py" \
    --annotation_dir "${DS_CHROM_DIR}/annotation" \
    --chrom "${CHROM}" \
    --ref_fai "${REF_FAI}" \
    --outdir "${DS_CHROM_DIR}"

# ── Summary ──
ds_log ""
ds_log "════════════════════════════════════════════════════════"
ds_log " Downstream analysis complete for ${CHROM}"
ds_log " Output: ${DS_CHROM_DIR}/"
ds_log "════════════════════════════════════════════════════════"
ds_log ""
ds_log "Directory structure:"
find "${DS_CHROM_DIR}" -name "*.tsv" -type f | wc -l | xargs -I{} echo "  Total TSV files: {}"
for subdir in catalogs matrices annotation per_sample distances markers gene_tables phase_prep; do
    n=$(find "${DS_CHROM_DIR}/${subdir}" -name "*.tsv" -type f 2>/dev/null | wc -l)
    ds_log "  ${subdir}/: ${n} files"
done

ds_log ""
ds_log "Next: run plotting scripts (R) or use --skip_plots to defer"
ds_log "  Figures will need: run_plots.sh --chrom ${CHROM}"
