#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=16G
#SBATCH -t 0-02:00:00
#SBATCH -A lt200308
#SBATCH -J ghsl_prep
#SBATCH -o logs/ghsl_prep.%A_%a.out
#SBATCH -e logs/ghsl_prep.%A_%a.err

set -euo pipefail

# ============================================================
# run_prepare_ghsl_array.sh — Prepare Clair3 data for Snake 3
#
# SLURM array: one task per chromosome.
#
# Usage:
#   sbatch --array=1-28 run_prepare_ghsl_array.sh
#   bash run_prepare_ghsl_array.sh --chrom C_gar_LG01
# ============================================================

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
C3_ROOT="${BASE}/MODULE_4A_SNP_INDEL50_Clair3"
SCRIPTS="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="${BASE}/inversion_localpca_v7/phased_summaries"
REF_FAI="${BASE}/00-samples/fClaHyb_Gar_LG.fa.fai"
CHROM_LIST="${C3_ROOT}/meta/chromosome_list.txt"
TIER="all"   # "all" or "tier1"
CHROM=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --chrom)   CHROM="$2"; shift 2 ;;
        --outdir)  OUTDIR="$2"; shift 2 ;;
        --tier)    TIER="$2"; shift 2 ;;
        *) echo "Unknown: $1" >&2; exit 1 ;;
    esac
done

mkdir -p "$OUTDIR" "${C3_ROOT}/logs"

# Resolve chromosome
if [[ -z "$CHROM" ]]; then
    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        CHROM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CHROM_LIST")
    else
        echo "[ERROR] Specify --chrom or use SLURM array" >&2
        exit 1
    fi
fi

[[ -n "$CHROM" ]] || { echo "[ERROR] Could not resolve chromosome" >&2; exit 1; }

PP_DIR="${C3_ROOT}/postprocess_results/${CHROM}"
[[ -d "$PP_DIR" ]] || { echo "[ERROR] No postprocess results: $PP_DIR" >&2; exit 1; }

echo "[GHSL-PREP] Chromosome: $CHROM"
echo "[GHSL-PREP] Tier: $TIER"
echo "[GHSL-PREP] Input: $PP_DIR"
echo "[GHSL-PREP] Output: $OUTDIR"

python "${SCRIPTS}/prepare_clair3_for_ghsl.py" \
    --pp_results "$PP_DIR" \
    --chrom "$CHROM" \
    --outdir "$OUTDIR" \
    --ref_fai "$REF_FAI" \
    --tier "$TIER"

echo "[GHSL-PREP] Done: $CHROM"
