#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --mem=180G
#SBATCH -t 3-00:00:00
#SBATCH -A lt200308
#SBATCH -J clair3_all
#SBATCH -o logs/clair3_all.%A_%a.out
#SBATCH -e logs/clair3_all.%A_%a.err

set -euo pipefail

# ============================================================
# 01_run_clair3_discovery_allchr.sh
#
# Alternative to 01_run_clair3_discovery.sh.
# Uses the dispatch table (sample × chromosome) so a single
# --array covers ALL combinations.
#
# Usage:
#   # Total tasks = N_samples × N_chromosomes
#   NTASKS=$(tail -n +2 meta/dispatch_table.tsv | wc -l)
#   sbatch --array=1-${NTASKS} 01_run_clair3_discovery_allchr.sh
#
# The dispatch table (meta/dispatch_table.tsv) has columns:
#   TASK_ID  SAMPLE  BAM  CHROM  BED
# ============================================================

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
ROOT="${BASE}/MODULE_4A_SNP_INDEL50_Clair3"
DISPATCH="${ROOT}/meta/dispatch_table.tsv"

[[ -s "$DISPATCH" ]] || { echo "[ERROR] Missing dispatch table: $DISPATCH" >&2; exit 1; }

TASK_ID="${SLURM_ARRAY_TASK_ID:-}"
[[ -n "$TASK_ID" ]] || { echo "[ERROR] Need SLURM_ARRAY_TASK_ID" >&2; exit 1; }

# +1 for header
ROW=$(( TASK_ID + 1 ))

SAMPLE="$(awk -v r="$ROW" 'BEGIN{FS="\t"} NR==r{print $2}' "$DISPATCH")"
BAM="$(awk    -v r="$ROW" 'BEGIN{FS="\t"} NR==r{print $3}' "$DISPATCH")"
CHROM="$(awk  -v r="$ROW" 'BEGIN{FS="\t"} NR==r{print $4}' "$DISPATCH")"

[[ -n "$SAMPLE" && -n "$BAM" && -n "$CHROM" ]] || {
    echo "[ERROR] Could not parse dispatch row $TASK_ID" >&2
    exit 1
}

echo "[DISPATCH] TASK=$TASK_ID  SAMPLE=$SAMPLE  CHROM=$CHROM"

# Delegate to the main discovery script
bash "${ROOT}/01_run_clair3_discovery.sh" \
    --chrom "$CHROM" \
    --bam "$BAM" \
    --sample "$SAMPLE"
