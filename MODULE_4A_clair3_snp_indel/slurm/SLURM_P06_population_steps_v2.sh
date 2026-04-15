#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# run_population_steps_v2.sh — Population pipeline launcher
#
# Submits two SLURM jobs with dependency:
#   1. Shared population steps (STEP05-07) — single job
#   2. Per-sample array (STEP08-10) — one task per sample
#
# The array job waits for the shared job to complete.
#
# Usage:
#   bash run_population_steps_v2.sh --chrom C_gar_LG01
#   bash run_population_steps_v2.sh --chrom C_gar_LG01 --min_samples 5
#
# For a single sample (no SLURM, runs locally):
#   bash run_population_steps_v2.sh --chrom C_gar_LG01 --sample CGA097
#
# To run only the per-sample array (shared steps already done):
#   bash run_population_steps_v2.sh --chrom C_gar_LG01 --skip_shared
# ============================================================

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
ROOT="${BASE}/MODULE_4A_SNP_INDEL50_Clair3"
SCRIPTS="${ROOT}/postprocess_scripts"
OUTBASE="${ROOT}/postprocess_results"
BAM_MANIFEST="${BASE}/pa_roary_results/00_manifests/sample_bam_minimap2_vs_P99TLENMAPQ30.tsv"

CHROM=""
MIN_SAMPLES=2
SINGLE_SAMPLE=""
SKIP_SHARED=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --chrom)        CHROM="$2"; shift 2 ;;
        --min_samples)  MIN_SAMPLES="$2"; shift 2 ;;
        --sample)       SINGLE_SAMPLE="$2"; shift 2 ;;
        --skip_shared)  SKIP_SHARED=1; shift ;;
        --scripts)      SCRIPTS="$2"; shift 2 ;;
        *) echo "Unknown: $1" >&2; exit 1 ;;
    esac
done

[[ -n "$CHROM" ]] || { echo "[ERROR] Specify --chrom" >&2; exit 1; }

POP_DIR="${OUTBASE}/${CHROM}/_population"
mkdir -p "$POP_DIR" "${ROOT}/logs"

echo "════════════════════════════════════════════════════════"
echo "[LAUNCH] Population pipeline v2 — $CHROM"
echo "════════════════════════════════════════════════════════"

# ── Single sample mode (local, no SLURM) ──
if [[ -n "$SINGLE_SAMPLE" ]]; then
    echo "[MODE] Single sample: $SINGLE_SAMPLE"

    if [[ "$SKIP_SHARED" -eq 0 ]]; then
        echo "[SHARED] Running STEP05-07 …"
        bash "${SCRIPTS}/run_population_shared.sh" \
            --chrom "$CHROM" --min_samples "$MIN_SAMPLES"
    fi

    echo "[PER-SAMPLE] Running STEP08-10 for $SINGLE_SAMPLE …"
    bash "${SCRIPTS}/run_population_per_sample.sh" \
        --chrom "$CHROM" --sample "$SINGLE_SAMPLE"

    echo "[DONE] Single sample $SINGLE_SAMPLE"
    exit 0
fi

# ── Full SLURM array mode ──

# Step 1: Submit shared job (STEP05-07)
if [[ "$SKIP_SHARED" -eq 0 ]]; then
    SHARED_JID=$(sbatch --parsable \
        "${SCRIPTS}/run_population_shared.sh" \
        --chrom "$CHROM" --min_samples "$MIN_SAMPLES" \
        --bam_manifest "$BAM_MANIFEST"
    )
    echo "[SHARED] Submitted: job $SHARED_JID (STEP05-07)"
else
    SHARED_JID=""
    echo "[SHARED] Skipped (--skip_shared)"
fi

# Step 2: Count samples for array
SAMPLE_LIST="${POP_DIR}/sample_list_for_array.txt"

if [[ "$SKIP_SHARED" -eq 1 ]]; then
    # Build sample list now (shared step would normally do this)
    PER_SAMPLE_DIR="${OUTBASE}/${CHROM}"
    ls -d "$PER_SAMPLE_DIR"/*/ 2>/dev/null \
        | while read d; do
            s=$(basename "$d")
            [[ "$s" == _* ]] && continue
            [[ -s "${d}/all_classified_step03.tsv" ]] || continue
            echo "$s"
        done | sort > "$SAMPLE_LIST"
fi

# If shared job was submitted, we need to wait for it before we know N_SAMPLES.
# Use a wrapper approach: submit the array with a generous --array range,
# and the worker script will exit cleanly if its line doesn't exist.
if [[ -f "$SAMPLE_LIST" ]]; then
    N_SAMPLES=$(wc -l < "$SAMPLE_LIST")
else
    # Estimate from directory count
    N_SAMPLES=$(ls -d "${OUTBASE}/${CHROM}"/*/ 2>/dev/null | grep -v '/_' | wc -l)
fi

[[ "$N_SAMPLES" -gt 0 ]] || { echo "[ERROR] No samples found" >&2; exit 1; }

echo "[ARRAY] Samples: $N_SAMPLES"

# Step 3: Submit per-sample array with dependency
DEP_ARG=""
[[ -n "$SHARED_JID" ]] && DEP_ARG="--dependency=afterok:${SHARED_JID}"

ARRAY_JID=$(sbatch --parsable \
    --array=1-${N_SAMPLES}%50 \
    $DEP_ARG \
    "${SCRIPTS}/run_population_per_sample.sh" \
    --chrom "$CHROM"
)
echo "[ARRAY] Submitted: job $ARRAY_JID (${N_SAMPLES} tasks, STEP08-10)"
[[ -n "$SHARED_JID" ]] && echo "  Depends on: $SHARED_JID"

echo ""
echo "════════════════════════════════════════════════════════"
echo "[LAUNCH] All jobs submitted for $CHROM"
echo "  Shared (STEP05-07): ${SHARED_JID:-skipped}"
echo "  Array  (STEP08-10): $ARRAY_JID (${N_SAMPLES} tasks)"
echo ""
echo "  Monitor: squeue -u \$USER"
echo "  Logs:    ${ROOT}/logs/c3pop_shared.*.out"
echo "           ${ROOT}/logs/c3pop_ps.${ARRAY_JID}_*.out"
echo "════════════════════════════════════════════════════════"
