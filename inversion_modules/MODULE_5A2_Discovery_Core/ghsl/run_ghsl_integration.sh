#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# run_ghsl_integration.sh — Clair3 → Snake 3 GHSL integration
#
# Phase 1: Prepare Clair3 phased data (per-chromosome array)
#           Produces: TSV, VCF, BED, marker files
# Phase 2: Run Snake 3 v2 with adaptive thresholds
#
# Usage:
#   # Full run (all chromosomes):
#   bash run_ghsl_integration.sh
#
#   # Single chromosome:
#   bash run_ghsl_integration.sh --chrom C_gar_LG01
#
#   # Skip prep (already done), just run Snake 3:
#   bash run_ghsl_integration.sh --skip_prep
#
#   # Run prep only (no Snake 3):
#   bash run_ghsl_integration.sh --prep_only
# ============================================================

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
C3_ROOT="${BASE}/MODULE_4A_SNP_INDEL50_Clair3"
INV_ROOT="${BASE}/inversion_localpca_v7"
CODEBASE="${BASE}/inversion_codebase_v8.4"
SCRIPTS="$(cd "$(dirname "$0")" && pwd)"

REF_FAI="${BASE}/00-samples/fClaHyb_Gar_LG.fa.fai"
SAMPLES_IND="${INV_ROOT}/samples.ind"
CHROM_LIST="${C3_ROOT}/meta/chromosome_list.txt"

PHASED_OUTDIR="${INV_ROOT}/phased_summaries"
SNAKE3_OUTDIR="${INV_ROOT}/snake3_v2_results"
STEP10_PREFIX="${INV_ROOT}/step10_output"   # adjust to your actual prefix

RSCRIPT_BIN="/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript"

CHROM=""
SKIP_PREP=0
PREP_ONLY=0
TIER="all"
MAX_CONCURRENT=28

while [[ $# -gt 0 ]]; do
    case "$1" in
        --chrom)          CHROM="$2"; shift 2 ;;
        --skip_prep)      SKIP_PREP=1; shift ;;
        --prep_only)      PREP_ONLY=1; shift ;;
        --tier)           TIER="$2"; shift 2 ;;
        --step10_prefix)  STEP10_PREFIX="$2"; shift 2 ;;
        --max_concurrent) MAX_CONCURRENT="$2"; shift 2 ;;
        *) echo "Unknown: $1" >&2; exit 1 ;;
    esac
done

mkdir -p "$PHASED_OUTDIR" "$SNAKE3_OUTDIR" "${C3_ROOT}/logs"

echo "════════════════════════════════════════════════════════"
echo "[GHSL] Clair3 → Snake 3 integration"
echo "  Clair3 results: $C3_ROOT"
echo "  Phased output:  $PHASED_OUTDIR"
echo "  Snake 3 output: $SNAKE3_OUTDIR"
echo "  Tier filter:    $TIER"
echo "  Skip prep:      $SKIP_PREP"
echo "  Prep only:      $PREP_ONLY"
echo "════════════════════════════════════════════════════════"

# ── Phase 1: Prepare Clair3 data ──
PREP_JID=""
if [[ "$SKIP_PREP" -eq 0 ]]; then
    if [[ -n "$CHROM" ]]; then
        # Single chromosome (local)
        echo "[PHASE 1] Preparing $CHROM …"
        python "${SCRIPTS}/prepare_clair3_for_ghsl.py" \
            --pp_results "${C3_ROOT}/postprocess_results/${CHROM}" \
            --chrom "$CHROM" \
            --outdir "$PHASED_OUTDIR" \
            --ref_fai "$REF_FAI" \
            --tier "$TIER"
    else
        # All chromosomes via SLURM array
        N_CHROMS=$(wc -l < "$CHROM_LIST")
        PREP_JID=$(sbatch --parsable \
            --array=1-${N_CHROMS}%${MAX_CONCURRENT} \
            "${SCRIPTS}/run_prepare_ghsl_array.sh" \
            --outdir "$PHASED_OUTDIR" \
            --tier "$TIER"
        )
        echo "[PHASE 1] Prep array submitted: job $PREP_JID ($N_CHROMS chromosomes)"
    fi
fi

if [[ "$PREP_ONLY" -eq 1 ]]; then
    echo "[DONE] Prep only mode — Snake 3 not launched."
    exit 0
fi

# ── Phase 2: Run Snake 3 v2 ──
if [[ -n "$CHROM" && "$SKIP_PREP" -eq 0 ]]; then
    # Single chromosome, prep already ran locally
    echo "[PHASE 2] Running Snake 3 v2 for $CHROM …"
    "${RSCRIPT_BIN}" "${SCRIPTS}/STEP10h_snake3_ghsl_haplotype_contrast_v2.R" \
        "$STEP10_PREFIX" \
        "$PHASED_OUTDIR" \
        "$SAMPLES_IND" \
        "$SNAKE3_OUTDIR"
elif [[ -n "$PREP_JID" ]]; then
    # Submit Snake 3 as dependency
    cat > "/tmp/run_snake3_v2_$$.sh" << JOBEOF
#!/bin/bash
#SBATCH -p compute
#SBATCH -N 1 -n 8 --mem=32G
#SBATCH -t 0-12:00:00
#SBATCH -A lt200308
#SBATCH -J snake3_v2
#SBATCH -o ${C3_ROOT}/logs/snake3_v2.%j.out
#SBATCH -e ${C3_ROOT}/logs/snake3_v2.%j.err
set -euo pipefail
source ~/.bashrc
mamba activate assembly
${RSCRIPT_BIN} ${SCRIPTS}/STEP10h_snake3_ghsl_haplotype_contrast_v2.R \\
    ${STEP10_PREFIX} \\
    ${PHASED_OUTDIR} \\
    ${SAMPLES_IND} \\
    ${SNAKE3_OUTDIR}
JOBEOF

    SNAKE3_JID=$(sbatch --parsable \
        --dependency=afterok:${PREP_JID} \
        "/tmp/run_snake3_v2_$$.sh"
    )
    echo "[PHASE 2] Snake 3 v2 submitted: job $SNAKE3_JID (depends on $PREP_JID)"
    rm -f "/tmp/run_snake3_v2_$$.sh"
else
    # Skip prep mode, run directly
    echo "[PHASE 2] Running Snake 3 v2 …"
    "${RSCRIPT_BIN}" "${SCRIPTS}/STEP10h_snake3_ghsl_haplotype_contrast_v2.R" \
        "$STEP10_PREFIX" \
        "$PHASED_OUTDIR" \
        "$SAMPLES_IND" \
        "$SNAKE3_OUTDIR"
fi

echo ""
echo "════════════════════════════════════════════════════════"
echo "[GHSL] Integration complete."
echo "  Phased data:    $PHASED_OUTDIR/"
echo "  Snake 3 output: $SNAKE3_OUTDIR/"
echo "════════════════════════════════════════════════════════"
