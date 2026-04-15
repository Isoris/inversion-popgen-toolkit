#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# run_full_pipeline.sh
#
# Master controller for the complete Clair3 discovery →
# post-processing pipeline.
#
# Phases:
#   PHASE 1: Prepare inputs (BEDs, manifests, dispatch table)
#   PHASE 2: Run Clair3 discovery (all samples × all chromosomes)
#   PHASE 3: Validate discovery outputs
#   PHASE 4: Launch post-processing (per-sample + population)
#
# Usage:
#   # Full pipeline, all chromosomes:
#   bash run_full_pipeline.sh --all_chroms
#
#   # Single chromosome (e.g. pilot):
#   bash run_full_pipeline.sh --chrom C_gar_LG28
#
#   # Just one phase:
#   bash run_full_pipeline.sh --phase prepare
#   bash run_full_pipeline.sh --phase discover --chrom C_gar_LG28
#   bash run_full_pipeline.sh --phase validate --chrom C_gar_LG28
#   bash run_full_pipeline.sh --phase postprocess --chrom C_gar_LG28
# ============================================================

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
ROOT="${BASE}/MODULE_4A_SNP_INDEL50_Clair3"

PHASE="all"
TARGET_CHROM=""
ALL_CHROMS=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --phase)      PHASE="$2"; shift 2 ;;
        --chrom)      TARGET_CHROM="$2"; shift 2 ;;
        --all_chroms) ALL_CHROMS=1; shift ;;
        *) echo "Unknown: $1" >&2; exit 1 ;;
    esac
done

echo ""
echo "╔══════════════════════════════════════════════════════════╗"
echo "║  Clair3 Discovery + Post-Processing Pipeline            ║"
echo "║  Catfish population genomics                            ║"
echo "╚══════════════════════════════════════════════════════════╝"
echo ""

# ── PHASE 1: Prepare ─────────────────────────────────────────────────────────
if [[ "$PHASE" == "all" || "$PHASE" == "prepare" ]]; then
    echo "━━━ PHASE 1: Prepare inputs ━━━"
    bash "${ROOT}/00_prepare_discovery_inputs.sh"
    echo ""
fi

# Determine chromosomes
if [[ "$ALL_CHROMS" -eq 1 ]]; then
    CHROMS=($(cat "${ROOT}/meta/chromosome_list.txt"))
elif [[ -n "$TARGET_CHROM" ]]; then
    CHROMS=("$TARGET_CHROM")
else
    if [[ -f "${ROOT}/meta/chromosome_list.txt" ]]; then
        CHROMS=($(cat "${ROOT}/meta/chromosome_list.txt"))
    else
        echo "[ERROR] Specify --chrom or --all_chroms" >&2
        exit 1
    fi
fi

N_SAMPLES=$(tail -n +2 "${ROOT}/meta/clair3_sample_manifest.tsv" 2>/dev/null | wc -l || echo 0)

# ── PHASE 2: Discovery ──────────────────────────────────────────────────────
if [[ "$PHASE" == "all" || "$PHASE" == "discover" ]]; then
    echo "━━━ PHASE 2: Clair3 discovery ━━━"
    echo "  Samples: $N_SAMPLES"
    echo "  Chromosomes: ${#CHROMS[@]}"
    echo ""

    DISC_JIDS=()

    for CHROM in "${CHROMS[@]}"; do
        JID=$(sbatch --parsable \
            --array=1-${N_SAMPLES} \
            "${ROOT}/01_run_clair3_discovery.sh" \
            --chrom "$CHROM" \
        )
        DISC_JIDS+=("$JID")
        echo "  $CHROM: submitted array job $JID (${N_SAMPLES} tasks)"
    done

    echo ""
    echo "  Discovery jobs submitted: ${#DISC_JIDS[@]} arrays"
    echo "  Total tasks: $(( N_SAMPLES * ${#CHROMS[@]} ))"
    echo ""
    echo "  Wait for completion before proceeding."
    echo "  Monitor: squeue -u \$USER | grep clair3"
    echo ""

    if [[ "$PHASE" == "discover" ]]; then
        echo "  Next: bash run_full_pipeline.sh --phase validate $(
            [[ "$ALL_CHROMS" -eq 1 ]] && echo '--all_chroms' || echo "--chrom $TARGET_CHROM")"
        exit 0
    fi
fi

# ── PHASE 3: Validate ───────────────────────────────────────────────────────
if [[ "$PHASE" == "all" || "$PHASE" == "validate" ]]; then
    echo "━━━ PHASE 3: Validate discovery outputs ━━━"

    if [[ "$ALL_CHROMS" -eq 1 ]]; then
        bash "${ROOT}/02_validate_discovery_outputs.sh" --all_chroms
    else
        for CHROM in "${CHROMS[@]}"; do
            bash "${ROOT}/02_validate_discovery_outputs.sh" --chrom "$CHROM"
        done
    fi

    echo ""
fi

# ── PHASE 4: Post-processing ────────────────────────────────────────────────
if [[ "$PHASE" == "all" || "$PHASE" == "postprocess" ]]; then
    echo "━━━ PHASE 4: Post-processing pipeline ━━━"

    if [[ "$ALL_CHROMS" -eq 1 ]]; then
        bash "${ROOT}/03_launch_postprocessing.sh" --all_chroms
    else
        for CHROM in "${CHROMS[@]}"; do
            bash "${ROOT}/03_launch_postprocessing.sh" --chrom "$CHROM"
        done
    fi
fi

echo ""
echo "╔══════════════════════════════════════════════════════════╗"
echo "║  Pipeline dispatched.                                   ║"
echo "║  Monitor: squeue -u \$USER                              ║"
echo "║  Logs:    ${ROOT}/logs/                                 ║"
echo "╚══════════════════════════════════════════════════════════╝"
