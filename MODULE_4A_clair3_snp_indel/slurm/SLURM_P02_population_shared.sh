#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=32G
#SBATCH -t 0-04:00:00
#SBATCH -A lt200308
#SBATCH -J c3pop_shared
#SBATCH -o logs/c3pop_shared.%j.out
#SBATCH -e logs/c3pop_shared.%j.err

set -euo pipefail

# ============================================================
# run_population_shared.sh — STEP05-06 (population clustering)
#
# Runs the shared clustering steps. STEP07 is now a separate
# array job (run_step07a_array.sh + STEP07B_merge).
#
# Usage:
#   sbatch run_population_shared.sh --chrom C_gar_LG01
#   sbatch run_population_shared.sh --chrom C_gar_LG01 --from_step 6
# ============================================================

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
ROOT="${BASE}/MODULE_4A_SNP_INDEL50_Clair3"

SCRIPTS="${ROOT}/postprocess_scripts"
OUTBASE="${ROOT}/postprocess_results"
BAM_MANIFEST="${BASE}/pa_roary_results/00_manifests/sample_bam_minimap2_vs_P99TLENMAPQ30.tsv"
CHROM=""
MIN_SAMPLES=2
FROM_STEP=5

while [[ $# -gt 0 ]]; do
    case "$1" in
        --scripts)      SCRIPTS="$2"; shift 2 ;;
        --outbase)      OUTBASE="$2"; shift 2 ;;
        --bam_manifest) BAM_MANIFEST="$2"; shift 2 ;;
        --chrom)        CHROM="$2"; shift 2 ;;
        --min_samples)  MIN_SAMPLES="$2"; shift 2 ;;
        --from_step)    FROM_STEP="$2"; shift 2 ;;
        *) echo "Unknown: $1" >&2; exit 1 ;;
    esac
done

[[ -n "$CHROM" ]] || { echo "[ERROR] Specify --chrom" >&2; exit 1; }

PER_SAMPLE_DIR="${OUTBASE}/${CHROM}"
POP_DIR="${OUTBASE}/${CHROM}/_population"
mkdir -p "$POP_DIR" "${ROOT}/logs"

echo "════════════════════════════════════════"
echo "[POP-SHARED] CHROM=$CHROM  FROM_STEP=$FROM_STEP"
echo "════════════════════════════════════════"

# ── STEP05 ──
if [[ "$FROM_STEP" -le 5 ]]; then
    WEAK_MANIFEST="${POP_DIR}/weak_candidate_file_manifest.txt"
    find "$PER_SAMPLE_DIR" -maxdepth 2 -name "weak_indel_candidates.tsv" -size +0c \
        | sort > "$WEAK_MANIFEST"
    N_FILES=$(wc -l < "$WEAK_MANIFEST")
    echo "[POP] Found $N_FILES weak candidate files"

    if [[ "$N_FILES" -ge 2 ]]; then
        echo "[STEP05] Group weak candidates …"
        python "${SCRIPTS}/STEP05_group_weak_candidates_across_samples.py" \
            --manifest "$WEAK_MANIFEST" --outdir "$POP_DIR" --min_samples "$MIN_SAMPLES"
    else
        echo "[POP] Fewer than 2 weak candidate files — skipping STEP05"
    fi
else
    echo "[STEP05] Skipped (--from_step $FROM_STEP)"
fi

# ── STEP06 ──
if [[ "$FROM_STEP" -le 6 ]]; then
    POP_CATALOG="${POP_DIR}/weak_indel_population_catalog.tsv"
    if [[ -s "$POP_CATALOG" ]]; then
        echo "[STEP06] Prepare regenotype catalog …"
        python "${SCRIPTS}/STEP06_prepare_regenotype_catalog.py" \
            --catalog "$POP_CATALOG" \
            --outdir "$POP_DIR" --min_samples "$MIN_SAMPLES"
    else
        echo "[STEP06] No population catalog — skipping"
    fi
else
    echo "[STEP06] Skipped (--from_step $FROM_STEP)"
fi

# ── Build sample list for arrays ──
SAMPLE_LIST="${POP_DIR}/sample_list_for_array.txt"
ls -d "$PER_SAMPLE_DIR"/*/ 2>/dev/null \
    | while read d; do
        s=$(basename "$d")
        [[ "$s" == _* ]] && continue
        [[ -s "${d}/all_classified_step03.tsv" ]] || continue
        echo "$s"
    done | sort > "$SAMPLE_LIST"

N_SAMPLES=$(wc -l < "$SAMPLE_LIST")

# Check if catalog exists for STEP07
REGEN_CAT="${POP_DIR}/rescued_indel_regenotype_catalog.tsv"
N_REGEN=0
[[ -s "$REGEN_CAT" ]] && N_REGEN=$(tail -n +2 "$REGEN_CAT" | wc -l)

echo ""
echo "════════════════════════════════════════"
echo "[POP-SHARED] Complete."
echo "  Samples: $N_SAMPLES"
echo "  Regenotype candidates: $N_REGEN"
echo "  Sample list: $SAMPLE_LIST"
echo "════════════════════════════════════════"
