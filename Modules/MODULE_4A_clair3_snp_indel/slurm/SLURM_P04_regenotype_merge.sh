#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 80
#SBATCH --mem=237G
#SBATCH -t 0-02:00:00
#SBATCH -A lt200308
#SBATCH -J c3_s07b
#SBATCH -o logs/c3_s07b.%j.out
#SBATCH -e logs/c3_s07b.%j.err

set -euo pipefail

# ============================================================
# run_step07b_merge.sh — Merge STEP07A per-sample outputs
# Usage:
#   sbatch run_step07b_merge.sh --chrom C_gar_LG01
# ============================================================

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
ROOT="${BASE}/MODULE_4A_SNP_INDEL50_Clair3"
SCRIPTS="${ROOT}/postprocess_scripts"
OUTBASE="${ROOT}/postprocess_results"
CHROM=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --chrom) CHROM="$2"; shift 2 ;;
        *) echo "Unknown: $1" >&2; exit 1 ;;
    esac
done

[[ -n "$CHROM" ]] || { echo "[ERROR] Specify --chrom" >&2; exit 1; }

POP_DIR="${OUTBASE}/${CHROM}/_population"

python "${SCRIPTS}/STEP07B_merge_regenotype.py" \
    --indir "$POP_DIR" \
    --outdir "$POP_DIR"
