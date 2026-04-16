#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=16G
#SBATCH -t 0-04:00:00
#SBATCH -A lt200308
#SBATCH -J c3pop_ps
#SBATCH -o logs/c3pop_ps.%A_%a.out
#SBATCH -e logs/c3pop_ps.%A_%a.err

set -euo pipefail

# ============================================================
# run_population_per_sample.sh — STEP08-10 (per-sample, SLURM array)
#
# Runs final classification, marker package, and figures for
# one sample at a time.  Designed for SLURM array submission.
#
# Usage:
#   # Array mode (typical):
#   sbatch --array=1-226 run_population_per_sample.sh --chrom C_gar_LG01
#
#   # Single sample by line number:
#   bash run_population_per_sample.sh --chrom C_gar_LG01 --line 25
#
#   # Single sample by name:
#   bash run_population_per_sample.sh --chrom C_gar_LG01 --sample CGA097
# ============================================================

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
ROOT="${BASE}/MODULE_4A_SNP_INDEL50_Clair3"

SCRIPTS="${ROOT}/postprocess_scripts"
OUTBASE="${ROOT}/postprocess_results"
REF="${BASE}/00-samples/fClaHyb_Gar_LG.fa"
REF_FAI="${REF}.fai"
CHROM=""
LINE_NO=""
SAMPLE_NAME=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --scripts)  SCRIPTS="$2"; shift 2 ;;
        --outbase)  OUTBASE="$2"; shift 2 ;;
        --chrom)    CHROM="$2"; shift 2 ;;
        --line)     LINE_NO="$2"; shift 2 ;;
        --sample)   SAMPLE_NAME="$2"; shift 2 ;;
        *) echo "Unknown: $1" >&2; exit 1 ;;
    esac
done

[[ -n "$CHROM" ]] || { echo "[ERROR] Specify --chrom" >&2; exit 1; }

PER_SAMPLE_DIR="${OUTBASE}/${CHROM}"
POP_DIR="${OUTBASE}/${CHROM}/_population"
SAMPLE_LIST="${POP_DIR}/sample_list_for_array.txt"

# ── Resolve sample name ──
if [[ -n "$SAMPLE_NAME" ]]; then
    # Explicit sample name given
    SAMPLE="$SAMPLE_NAME"
elif [[ -n "$LINE_NO" ]]; then
    # Line number given
    SAMPLE=$(sed -n "${LINE_NO}p" "$SAMPLE_LIST")
elif [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    # SLURM array mode
    LINE_NO="${SLURM_ARRAY_TASK_ID}"
    SAMPLE=$(sed -n "${LINE_NO}p" "$SAMPLE_LIST")
else
    echo "[ERROR] Need --sample, --line, or SLURM_ARRAY_TASK_ID" >&2
    exit 1
fi

[[ -n "$SAMPLE" ]] || { echo "[ERROR] Could not resolve sample (line=$LINE_NO)" >&2; exit 1; }

SAMPLE_DIR="${PER_SAMPLE_DIR}/${SAMPLE}"
[[ -d "$SAMPLE_DIR" ]] || { echo "[ERROR] Sample dir not found: $SAMPLE_DIR" >&2; exit 1; }

CHR_LEN=$(awk -v c="$CHROM" '$1==c{print $2}' "$REF_FAI")
[[ -n "$CHR_LEN" ]] || CHR_LEN=33300000

echo "════════════════════════════════════════"
echo "[POP-PERSAMPLE] SAMPLE=$SAMPLE  CHROM=$CHROM  LINE=${LINE_NO:-direct}"
echo "════════════════════════════════════════"

# ── Check prerequisites ──
STEP03="${SAMPLE_DIR}/all_classified_step03.tsv"
[[ -s "$STEP03" ]] || { echo "[ERROR] Missing step03 for $SAMPLE: $STEP03" >&2; exit 1; }

# ── Regenotype data (from shared step) ──
REGEN_DETAIL="${POP_DIR}/regenotyped_rescued_indels.tsv"
REGEN_SUMMARY="${POP_DIR}/regenotype_cohort_summary.tsv"

REGEN_ARGS=""
[[ -s "$REGEN_DETAIL" ]]  && REGEN_ARGS="$REGEN_ARGS --regen $REGEN_DETAIL"
[[ -s "$REGEN_SUMMARY" ]] && REGEN_ARGS="$REGEN_ARGS --regen_summary $REGEN_SUMMARY"

PHASE_ARG=""
[[ -f "${SAMPLE_DIR}/all_variants_with_phase.tsv" ]] && \
    PHASE_ARG="--phase ${SAMPLE_DIR}/all_variants_with_phase.tsv"

# ── STEP08: Final classification ──
echo "[STEP08] Final classification for $SAMPLE …"
python "${SCRIPTS}/STEP08_final_classification.py" \
    --step03 "$STEP03" \
    --weak "${SAMPLE_DIR}/weak_indel_candidates.tsv" \
    $REGEN_ARGS \
    --blocks "${SAMPLE_DIR}/local_overlap_blocks.tsv" \
    --outdir "$SAMPLE_DIR" \
    --sample_id "$SAMPLE"

# ── STEP09: Marker handoff package ──
echo "[STEP09] Marker package for $SAMPLE …"
python "${SCRIPTS}/STEP09_marker_handoff_package.py" \
    --classified "${SAMPLE_DIR}/final_variant_classification.tsv" \
    $PHASE_ARG \
    --ref "$REF" \
    --outdir "${SAMPLE_DIR}/01_Use_for_markers_multiplex" \
    --sample_id "$SAMPLE" \
    --flank_bp 150

# ── STEP10: Publication figure ──
echo "[STEP10] Figure for $SAMPLE …"
SAMPLE_FIGDIR="${SAMPLE_DIR}/figures"
mkdir -p "$SAMPLE_FIGDIR"
Rscript "${SCRIPTS}/STEP10_publication_figure.R" \
    "${SAMPLE_DIR}/final_variant_classification.tsv" \
    "${SAMPLE_DIR}/weak_indel_candidates.tsv" \
    "$SAMPLE_FIGDIR" \
    "$CHR_LEN" \
    "$SAMPLE" 2>/dev/null || echo "[WARN] Figure failed for $SAMPLE"

echo "[DONE] $SAMPLE"
