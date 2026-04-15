#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# 03_launch_postprocessing.sh
#
# Pipeline phases:
#   Phase 1: Per-sample array   → STEP01-04+02B
#   Phase 2: Shared single job  → STEP05-06 (cluster + catalog)
#   Phase 3: Per-sample array   → STEP07A (regenotype 1 BAM each)
#   Phase 4: Shared single job  → STEP07B (merge regenotype results)
#   Phase 5: Per-sample array   → STEP08-10 (classify + markers + figures)
#
# Usage:
#   bash 03_launch_postprocessing.sh --chrom C_gar_LG01
#   bash 03_launch_postprocessing.sh --all_chroms
#
#   # Resume from any step:
#   bash 03_launch_postprocessing.sh --chrom C_gar_LG01 --step 7
#   bash 03_launch_postprocessing.sh --chrom C_gar_LG01 --step 8
#
# Step map:
#   --step 1   Phase 1+2+3+4+5  (full run)
#   --step 5   Phase 2+3+4+5    (STEP05 onward)
#   --step 6   Phase 2+3+4+5    (from STEP06)
#   --step 7   Phase 3+4+5      (STEP07A array + merge + STEP08-10)
#   --step 7b  Phase 4+5        (merge only + STEP08-10)
#   --step 8   Phase 5           (STEP08-10 per-sample only)
# ============================================================

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
ROOT="${BASE}/MODULE_4A_SNP_INDEL50_Clair3"
SCRIPTS="${ROOT}/postprocess_scripts"
MANIFEST="${ROOT}/meta/clair3_sample_manifest.tsv"
BAM_MANIFEST="${BASE}/pa_roary_results/00_manifests/sample_bam_minimap2_vs_P99TLENMAPQ30.tsv"
CHROM_LIST="${ROOT}/meta/chromosome_list.txt"

MODE="single"
TARGET_CHROM=""
FROM_STEP="1"
MAX_CONCURRENT=50
MIN_SAMPLES=2

while [[ $# -gt 0 ]]; do
    case "$1" in
        --chrom)          TARGET_CHROM="$2"; shift 2 ;;
        --all_chroms)     MODE="all"; shift ;;
        --scripts)        SCRIPTS="$2"; shift 2 ;;
        --step)           FROM_STEP="$2"; shift 2 ;;
        --max_concurrent) MAX_CONCURRENT="$2"; shift 2 ;;
        --min_samples)    MIN_SAMPLES="$2"; shift 2 ;;
        *) echo "Unknown: $1" >&2; exit 1 ;;
    esac
done

[[ -s "$MANIFEST" ]] || { echo "[ERROR] Missing manifest: $MANIFEST" >&2; exit 1; }
N_SAMPLES=$(tail -n +2 "$MANIFEST" | wc -l)

if [[ "$MODE" == "all" ]]; then
    [[ -s "$CHROM_LIST" ]] || { echo "[ERROR] Missing chromosome list" >&2; exit 1; }
    CHROMS=($(cat "$CHROM_LIST"))
elif [[ -n "$TARGET_CHROM" ]]; then
    CHROMS=("$TARGET_CHROM")
else
    echo "[ERROR] Specify --chrom or --all_chroms" >&2
    exit 1
fi

# Determine which phases to run
RUN_P1=1; RUN_P2=1; RUN_P3=1; RUN_P4=1; RUN_P5=1
P2_FROM=5

case "$FROM_STEP" in
    1|2|3|4)  ;;                                              # run all
    5)        RUN_P1=0; P2_FROM=5 ;;                          # Phase 2-5
    6)        RUN_P1=0; P2_FROM=6 ;;                          # Phase 2-5, skip STEP05
    7)        RUN_P1=0; RUN_P2=0 ;;                           # Phase 3-5
    7b|7B)    RUN_P1=0; RUN_P2=0; RUN_P3=0 ;;                # Phase 4-5 (merge+classify)
    8|9|10)   RUN_P1=0; RUN_P2=0; RUN_P3=0; RUN_P4=0 ;;     # Phase 5 only
    *)        echo "[ERROR] Invalid --step: $FROM_STEP" >&2; exit 1 ;;
esac

mkdir -p "${ROOT}/logs"

echo "════════════════════════════════════════════════════════"
echo "[LAUNCH] Post-processing pipeline v4"
echo "  Samples:        $N_SAMPLES"
echo "  Chromosomes:    ${#CHROMS[@]}"
echo "  Resume from:    --step $FROM_STEP"
echo "  Phase 1 (01-04):  $([ $RUN_P1 -eq 1 ] && echo 'RUN' || echo 'SKIP')"
echo "  Phase 2 (05-06):  $([ $RUN_P2 -eq 1 ] && echo "RUN (from $P2_FROM)" || echo 'SKIP')"
echo "  Phase 3 (07A):    $([ $RUN_P3 -eq 1 ] && echo 'RUN (array)' || echo 'SKIP')"
echo "  Phase 4 (07B):    $([ $RUN_P4 -eq 1 ] && echo 'RUN (merge)' || echo 'SKIP')"
echo "  Phase 5 (08-10):  $([ $RUN_P5 -eq 1 ] && echo 'RUN (array)' || echo 'SKIP')"
echo "════════════════════════════════════════════════════════"

for CHROM in "${CHROMS[@]}"; do
    echo ""
    echo "── $CHROM ──"

    VCF_DIR="${ROOT}/vcf/${CHROM}"
    OUTBASE="${ROOT}/postprocess_results"
    POP_DIR="${OUTBASE}/${CHROM}/_population"
    mkdir -p "$POP_DIR"

    PREV_JID=""

    # ── Phase 1: Per-sample STEP01-04+02B ──
    if [[ "$RUN_P1" -eq 1 ]]; then
        N_VCF=$(find "$VCF_DIR" -name "*.vcf.gz" 2>/dev/null | wc -l)
        if [[ "$N_VCF" -eq 0 ]]; then
            echo "  [WARN] No VCFs in $VCF_DIR — skipping $CHROM"
            continue
        fi
        PREV_JID=$(sbatch --parsable \
            --array=1-${N_SAMPLES}%${MAX_CONCURRENT} \
            "${SCRIPTS}/run_per_sample_steps_v2.sh" \
            --chrom "$CHROM" --vcf_dir "$VCF_DIR" \
            --bed "${ROOT}/beds/${CHROM}.bed"
        )
        echo "  Phase 1: job $PREV_JID (${N_SAMPLES} tasks, STEP01-04)"
    else
        echo "  Phase 1: skipped"
    fi

    # ── Phase 2: Shared STEP05-06 ──
    if [[ "$RUN_P2" -eq 1 ]]; then
        DEP=""; [[ -n "$PREV_JID" ]] && DEP="--dependency=afterok:${PREV_JID}"
        PREV_JID=$(sbatch --parsable $DEP \
            "${SCRIPTS}/run_population_shared.sh" \
            --chrom "$CHROM" --bam_manifest "$BAM_MANIFEST" \
            --min_samples "$MIN_SAMPLES" --from_step "$P2_FROM"
        )
        echo "  Phase 2: job $PREV_JID (STEP${P2_FROM}-06)"
    else
        echo "  Phase 2: skipped"
        # Build sample list if needed
        SAMPLE_LIST="${POP_DIR}/sample_list_for_array.txt"
        if [[ ! -f "$SAMPLE_LIST" ]]; then
            ls -d "${OUTBASE}/${CHROM}"/*/ 2>/dev/null \
                | while read d; do
                    s=$(basename "$d")
                    [[ "$s" == _* ]] && continue
                    [[ -s "${d}/all_classified_step03.tsv" ]] || continue
                    echo "$s"
                done | sort > "$SAMPLE_LIST"
        fi
    fi

    # ── Phase 3: STEP07A per-sample array ──
    if [[ "$RUN_P3" -eq 1 ]]; then
        DEP=""; [[ -n "$PREV_JID" ]] && DEP="--dependency=afterok:${PREV_JID}"
        PREV_JID=$(sbatch --parsable $DEP \
            --array=1-${N_SAMPLES}%${MAX_CONCURRENT} \
            "${SCRIPTS}/run_step07a_array.sh" \
            --chrom "$CHROM"
        )
        echo "  Phase 3: job $PREV_JID (${N_SAMPLES} tasks, STEP07A regenotype)"
    else
        echo "  Phase 3: skipped"
    fi

    # ── Phase 4: STEP07B merge ──
    if [[ "$RUN_P4" -eq 1 ]]; then
        DEP=""; [[ -n "$PREV_JID" ]] && DEP="--dependency=afterok:${PREV_JID}"
        PREV_JID=$(sbatch --parsable $DEP \
            "${SCRIPTS}/run_step07b_merge.sh" \
            --chrom "$CHROM"
        )
        echo "  Phase 4: job $PREV_JID (STEP07B merge)"
    else
        echo "  Phase 4: skipped"
    fi

    # ── Phase 5: Per-sample STEP08-10 array ──
    if [[ "$RUN_P5" -eq 1 ]]; then
        SAMPLE_LIST="${POP_DIR}/sample_list_for_array.txt"
        if [[ -f "$SAMPLE_LIST" ]]; then
            N_PS=$(wc -l < "$SAMPLE_LIST")
        else
            N_PS="$N_SAMPLES"
        fi

        DEP=""; [[ -n "$PREV_JID" ]] && DEP="--dependency=afterok:${PREV_JID}"
        PREV_JID=$(sbatch --parsable $DEP \
            --array=1-${N_PS}%${MAX_CONCURRENT} \
            "${SCRIPTS}/run_population_per_sample.sh" \
            --chrom "$CHROM"
        )
        echo "  Phase 5: job $PREV_JID (${N_PS} tasks, STEP08-10)"
    else
        echo "  Phase 5: skipped"
    fi

done

echo ""
echo "════════════════════════════════════════════════════════"
echo "[LAUNCH] All jobs submitted."
echo "  Monitor: squeue -u \$USER"
echo "  Logs:    ${ROOT}/logs/"
echo "════════════════════════════════════════════════════════"
