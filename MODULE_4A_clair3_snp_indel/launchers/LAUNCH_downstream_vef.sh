#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=32G
#SBATCH -t 0-06:00:00
#SBATCH -A lt200308
#SBATCH -J c3_vef
#SBATCH -o logs/c3_vef.%j.out
#SBATCH -e logs/c3_vef.%j.err
# ============================================================
# run_vef.sh — Run Variant Evidence Framework scoring
#
# Prerequisite: run_downstream.sh must have completed for this chrom.
#
# Usage:
#   bash run_vef.sh --chrom C_gar_LG01
#   sbatch run_vef.sh --chrom C_gar_LG01
#
# Optional inputs:
#   --csq_tsv    PATH   Pre-computed bcftools csq consequences
#   --roh_dir    PATH   Directory with per-sample .roh.bed files
#   --esm_tsv    PATH   ESM missense scores
# ============================================================
set -euo pipefail
source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${SCRIPT_DIR}/00_config.sh"

CHROM=""
CSQ_TSV=""
ROH_DIR=""
ESM_TSV=""
VTYPE="COMBINED"  # default: score all variants together

while [[ $# -gt 0 ]]; do
    case "$1" in
        --chrom)    CHROM="$2"; shift 2 ;;
        --csq_tsv)  CSQ_TSV="$2"; shift 2 ;;
        --roh_dir)  ROH_DIR="$2"; shift 2 ;;
        --esm_tsv)  ESM_TSV="$2"; shift 2 ;;
        --vtype)    VTYPE="$2"; shift 2 ;;
        *) echo "Unknown: $1" >&2; exit 1 ;;
    esac
done
[[ -n "$CHROM" ]] || { echo "[ERROR] Specify --chrom" >&2; exit 1; }

ds_init_dirs "$CHROM"
DS_DIR="${DS_ROOT}/${CHROM}"

# Check prerequisites
CATALOG="${DS_DIR}/catalogs/${VTYPE}_catalog.tsv"
GT_MATRIX="${DS_DIR}/matrices/GT_matrix.${VTYPE}.tsv"
[[ -f "$CATALOG" ]]  || { echo "[ERROR] Missing catalog: $CATALOG — run run_downstream.sh first" >&2; exit 1; }
[[ -f "$GT_MATRIX" ]] || { echo "[ERROR] Missing GT matrix: $GT_MATRIX" >&2; exit 1; }

# Try to find ROH directory if not specified
if [[ -z "$ROH_DIR" ]]; then
    # Check common locations
    for candidate in \
        "${BASE}/MODULE_3_HetROH/roh_beds" \
        "${BASE}/MODULE_3_HetROH/ngsF-HMM_results" \
        "${BASE}/pa_roary_results/06_het_roh/roh_beds"; do
        if [[ -d "$candidate" ]]; then
            ROH_DIR="$candidate"
            ds_log "Auto-detected ROH dir: $ROH_DIR"
            break
        fi
    done
fi

# Try to find functional class from DELLY annotation
FUNC_TSV=""
DELLY_FUNC="${BASE}/MODULE_5_DELLY_DEL/10_annotation/catalog_226.functional_class.tsv"
[[ -f "$DELLY_FUNC" ]] && FUNC_TSV="$DELLY_FUNC"

ds_log "════════════════════════════════════════════════════════"
ds_log " Variant Evidence Framework: ${CHROM} (${VTYPE})"
ds_log "════════════════════════════════════════════════════════"
ds_log "  Catalog:     $CATALOG"
ds_log "  GT matrix:   $GT_MATRIX"
ds_log "  CSQ:         ${CSQ_TSV:-not provided}"
ds_log "  ROH:         ${ROH_DIR:-not provided}"
ds_log "  ESM:         ${ESM_TSV:-not provided}"
ds_log "  Functional:  ${FUNC_TSV:-not provided}"

# ── Step 20: Annotate variants ──
ds_log ""
ds_log "=== 20: Annotate variant consequences ==="

ANNOT_ARGS="--catalog $CATALOG --gt_matrix $GT_MATRIX --chrom $CHROM --outdir $DS_DIR"
[[ -n "$CSQ_TSV" ]]  && ANNOT_ARGS="$ANNOT_ARGS --csq_tsv $CSQ_TSV"
[[ -n "$ROH_DIR" ]]  && ANNOT_ARGS="$ANNOT_ARGS --roh_dir $ROH_DIR"
[[ -n "$ESM_TSV" ]]  && ANNOT_ARGS="$ANNOT_ARGS --esm_tsv $ESM_TSV"
[[ -n "$FUNC_TSV" ]] && ANNOT_ARGS="$ANNOT_ARGS --functional_class_tsv $FUNC_TSV"

python3 "${SCRIPT_DIR}/20_annotate_variant_consequences.py" $ANNOT_ARGS

# ── Step 21: Score breeding concern ──
ds_log ""
ds_log "=== 21: Score breeding concern ==="

python3 "${SCRIPT_DIR}/21_score_breeding_concern.py" \
    --annotated_variants "${DS_DIR}/vef/annotated_variants.tsv" \
    --per_sample_genotypes "${DS_DIR}/vef/per_sample_genotypes.tsv" \
    --outdir "$DS_DIR"

# ── Summary ──
ds_log ""
ds_log "════════════════════════════════════════════════════════"
ds_log " VEF scoring complete for ${CHROM}"
ds_log "════════════════════════════════════════════════════════"
ds_log ""
ds_log "Outputs in ${DS_DIR}/vef/:"
ls -lh "${DS_DIR}/vef/"*.tsv 2>/dev/null | while read line; do ds_log "  $line"; done
ds_log ""
ds_log "Key files:"
ds_log "  scored_variants.tsv       — full per-variant × per-sample scored table"
ds_log "  variant_master_table.tsv  — per-variant summary with worst-case BC"
ds_log "  fish_burden_summary.tsv   — per-fish aggregate deleterious burden"
ds_log "  high_concern_report.tsv   — BC3+BC4 variants for manual review"
ds_log ""
ds_log "To enrich results, provide:"
ds_log "  --csq_tsv   bcftools csq output for consequence annotation"
ds_log "  --roh_dir   ROH BEDs from ngsF-HMM for zygosity context"
ds_log "  --esm_tsv   ESM scores for missense functional prediction"
