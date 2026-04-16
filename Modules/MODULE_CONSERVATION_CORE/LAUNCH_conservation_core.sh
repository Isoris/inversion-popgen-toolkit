#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# LAUNCH_conservation_core.sh — Run the CORE conservation pipeline
# =============================================================================
# Dependency chain:
#   STEP_00 → STEP_03 → { STEP_04 → STEP_04B }
#                        { STEP_04 → STEP_05 }
#                        { STEP_04 → STEP_05B (GPU) }
#                      → STEP_14 → STEP_15
#
# Three scoring tools: bcftools csq + SIFT4G + VESM
# No cross-species alignment.
#
# Usage:
#   bash LAUNCH_conservation_core.sh [--from STEP] [--dry-run] [--no-vesm]
#
# Options:
#   --from N     Start from step N (skip earlier steps if outputs exist)
#   --dry-run    Print sbatch commands without submitting
#   --no-vesm    Skip VESM (no GPU needed; score with csq + SIFT only)
#   --chr LIST   Comma-separated chromosomes (default: all 28)
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config/pipeline.config.sh"

# Parse arguments
FROM_STEP=0
DRY_RUN=false
NO_VESM=false
CHR_LIST=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --from)    FROM_STEP="$2"; shift 2 ;;
        --dry-run) DRY_RUN=true; shift ;;
        --no-vesm) NO_VESM=true; shift ;;
        --chr)     CHR_LIST="$2"; shift 2 ;;
        *)         echo "Unknown option: $1"; exit 1 ;;
    esac
done

cons_init_dirs

submit() {
    local script="$1"
    local dep="${2:-}"
    local extra="${3:-}"

    local cmd="sbatch"
    [[ -n "$dep" ]] && cmd="$cmd --dependency=afterok:${dep}"
    [[ -n "$extra" ]] && cmd="$cmd ${extra}"
    cmd="$cmd ${SCRIPT_DIR}/${script}"

    if $DRY_RUN; then
        echo "[DRY-RUN] $cmd"
        echo "DRY_$(basename "$script" .sh)_$$"
    else
        local out
        out=$($cmd)
        local jid
        jid=$(echo "$out" | grep -oP '\d+')
        echo "  Submitted ${script} → job ${jid}"
        echo "$jid"
    fi
}

echo "=== LAUNCH: Conservation Pipeline (CORE-ONLY) ==="
echo "  Script dir: ${SCRIPT_DIR}"
echo "  From step:  ${FROM_STEP}"
echo "  VESM:       $(if $NO_VESM; then echo SKIP; else echo YES; fi)"
echo ""

# ── STEP 00: Setup ───
JID_00=""
if [[ $FROM_STEP -le 0 ]]; then
    echo "[STEP 00] Reference preparation..."
    JID_00=$(submit STEP_00_setup.sh)
fi

# ── STEP 03: Merge + Normalize ───
JID_03=""
if [[ $FROM_STEP -le 3 ]]; then
    echo "[STEP 03] Merge + normalize Clair3 VCFs..."
    JID_03=$(submit STEP_03_merge_normalize.sh "$JID_00")
fi

# ── STEP 04: SnpEff ───
JID_04=""
if [[ $FROM_STEP -le 4 ]]; then
    echo "[STEP 04] SnpEff annotation..."
    JID_04=$(submit STEP_04_snpeff_annotate.sh "$JID_03")
fi

# ── STEP 04B: bcftools csq (parallel with 05, 05B) ───
JID_04B=""
if [[ $FROM_STEP -le 4 ]]; then
    echo "[STEP 04B] bcftools csq..."
    JID_04B=$(submit STEP_04B_bcftools_csq.sh "$JID_04")
fi

# ── STEP 05: SIFT4G (parallel with 04B, 05B) ───
JID_05=""
if [[ $FROM_STEP -le 5 ]]; then
    echo "[STEP 05] SIFT4G database + annotation..."
    JID_05=$(submit STEP_05_sift4g.sh "$JID_04")
fi

# ── STEP 05B: VESM (GPU, parallel with 04B, 05) ───
JID_05B=""
if [[ $FROM_STEP -le 5 ]] && ! $NO_VESM; then
    echo "[STEP 05B] VESM protein language model (GPU)..."
    JID_05B=$(submit STEP_05B_vesm_predict.sh "$JID_04" "--gres=gpu:1")
fi

# ── Wait for annotation convergence ───
DEPS_14=""
for jid in $JID_04B $JID_05 $JID_05B; do
    [[ -n "$jid" ]] && DEPS_14="${DEPS_14:+${DEPS_14}:}${jid}"
done

# ── STEP 14: Score variants ───
JID_14=""
if [[ $FROM_STEP -le 14 ]]; then
    echo "[STEP 14] Score variants (csq + SIFT + VESM + splice)..."
    JID_14=$(submit STEP_14_score_variants.sh "$DEPS_14")
fi

# ── STEP 15: Burden tables ───
JID_15=""
if [[ $FROM_STEP -le 15 ]]; then
    echo "[STEP 15] Burden tables + genotype matrix..."
    JID_15=$(submit STEP_15_burden_tables.sh "$JID_14")
fi

echo ""
echo "=== Pipeline submitted ==="
echo ""
echo "  Dependency chain:"
echo "    00 → 03 → 04 → { 04B, 05, 05B } → 14 → 15"
echo ""
echo "  Three output files (after STEP 15):"
echo "    ${DIR_MERGED}/variant_master_final.tsv"
echo "    ${DIR_BURDEN}/deleterious_variant_genotype_matrix.tsv"
echo "    ${DIR_BURDEN}/deleterious_burden_per_individual.tsv"
echo ""
echo "  Feed into C01f_c:"
echo "    Rscript STEP_C01f_c_burden_regression.R \\"
echo "      --burden ${DIR_BURDEN}/deleterious_burden_per_individual.tsv \\"
echo "      --del-geno ${DIR_BURDEN}/deleterious_variant_genotype_matrix.tsv \\"
echo "      --del-master ${DIR_MERGED}/variant_master_final.tsv \\"
echo "      --decomp <your_decomposition.tsv.gz> \\"
echo "      --catalogue <your_classification.tsv>"
echo ""
if [[ -n "$JID_15" ]] && ! $DRY_RUN; then
    echo "  Monitor: squeue -u \$USER | grep CONS"
    echo "  Final job: $JID_15"
fi
