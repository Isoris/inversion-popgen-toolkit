#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=16G
#SBATCH -t 0-06:00:00
#SBATCH -A lt200308
#SBATCH -J c3_s07a
#SBATCH -o logs/c3_s07a.%A_%a.out
#SBATCH -e logs/c3_s07a.%A_%a.err

set -euo pipefail

# ============================================================
# run_step07a_array.sh — Per-sample regenotyping (SLURM array)
#
# Each task scans ONE sample's BAM against the shared catalog.
#
# Usage:
#   sbatch --array=1-226%50 run_step07a_array.sh --chrom C_gar_LG01
#   bash run_step07a_array.sh --chrom C_gar_LG01 --line 25
#   bash run_step07a_array.sh --chrom C_gar_LG01 --sample CGA097
# ============================================================

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
ROOT="${BASE}/MODULE_4A_SNP_INDEL50_Clair3"
SCRIPTS="${ROOT}/postprocess_scripts"
BAM_MANIFEST="${BASE}/pa_roary_results/00_manifests/sample_bam_minimap2_vs_P99TLENMAPQ30.tsv"
OUTBASE="${ROOT}/postprocess_results"

CHROM=""
LINE_NO=""
SAMPLE_NAME=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --chrom)   CHROM="$2"; shift 2 ;;
        --line)    LINE_NO="$2"; shift 2 ;;
        --sample)  SAMPLE_NAME="$2"; shift 2 ;;
        *) echo "Unknown: $1" >&2; exit 1 ;;
    esac
done

[[ -n "$CHROM" ]] || { echo "[ERROR] Specify --chrom" >&2; exit 1; }

POP_DIR="${OUTBASE}/${CHROM}/_population"
REGEN_CAT="${POP_DIR}/rescued_indel_regenotype_catalog.tsv"
[[ -s "$REGEN_CAT" ]] || { echo "[ERROR] No catalog: $REGEN_CAT" >&2; exit 1; }

# Resolve sample
if [[ -n "$SAMPLE_NAME" ]]; then
    SAMPLE="$SAMPLE_NAME"
    BAM=$(awk -F'\t' -v s="$SAMPLE" 'NR>1 && $1==s{print $3}' "$BAM_MANIFEST")
elif [[ -n "$LINE_NO" ]]; then
    ROW=$(( LINE_NO + 1 ))
    SAMPLE=$(awk -v r="$ROW" 'BEGIN{FS="\t"} NR==r{print $1}' "$BAM_MANIFEST")
    BAM=$(awk -v r="$ROW" 'BEGIN{FS="\t"} NR==r{print $3}' "$BAM_MANIFEST")
elif [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    LINE_NO="${SLURM_ARRAY_TASK_ID}"
    ROW=$(( LINE_NO + 1 ))
    SAMPLE=$(awk -v r="$ROW" 'BEGIN{FS="\t"} NR==r{print $1}' "$BAM_MANIFEST")
    BAM=$(awk -v r="$ROW" 'BEGIN{FS="\t"} NR==r{print $3}' "$BAM_MANIFEST")
else
    echo "[ERROR] Need --sample, --line, or SLURM_ARRAY_TASK_ID" >&2
    exit 1
fi

[[ -n "$SAMPLE" ]] || { echo "[ERROR] Could not resolve sample" >&2; exit 1; }
[[ -s "$BAM" ]]    || { echo "[ERROR] BAM not found: $BAM" >&2; exit 1; }

echo "[STEP07A] SAMPLE=$SAMPLE  CHROM=$CHROM  BAM=$BAM"

python "${SCRIPTS}/STEP07A_regenotype_one_sample.py" \
    --catalog "$REGEN_CAT" \
    --sample_id "$SAMPLE" \
    --bam "$BAM" \
    --outdir "$POP_DIR" \
    --chrom "$CHROM"
