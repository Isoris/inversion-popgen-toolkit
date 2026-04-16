#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# 00_prepare_discovery_inputs.sh
#
# Prepares everything needed for cohort-wide Clair3 discovery:
#   1. Per-chromosome BED files from the reference FAI
#   2. A Clair3-specific manifest linking samples to filtered BAMs
#   3. A chromosome list for array dispatching
#   4. Validates all BAMs exist and are indexed
#
# Usage:
#   bash 00_prepare_discovery_inputs.sh [--dry_run]
# ============================================================

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
ROOT="${BASE}/MODULE_4A_SNP_INDEL50_Clair3"
REF="${BASE}/00-samples/fClaHyb_Gar_LG.fa"
REF_FAI="${REF}.fai"
BAM_MANIFEST="${BASE}/pa_roary_results/00_manifests/sample_bam_minimap2_vs_P99TLENMAPQ30.tsv"

DRY_RUN=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry_run) DRY_RUN=1; shift ;;
        *) echo "Unknown: $1" >&2; exit 1 ;;
    esac
done

# ── directories ──────────────────────────────────────────────────────────────
mkdir -p "${ROOT}/beds" "${ROOT}/meta" "${ROOT}/logs" "${ROOT}/vcf" "${ROOT}/gvcf"

echo "════════════════════════════════════════════════════════"
echo "[PREP] Clair3 discovery input preparation"
echo "════════════════════════════════════════════════════════"

# ── 1. Validate reference ────────────────────────────────────────────────────
[[ -s "$REF" ]]     || { echo "[ERROR] Missing reference: $REF" >&2; exit 1; }
[[ -s "$REF_FAI" ]] || { echo "[ERROR] Missing FAI: $REF_FAI" >&2; exit 1; }

N_CHROMS=$(wc -l < "$REF_FAI")
echo "[PREP] Reference: $REF"
echo "[PREP] Chromosomes in FAI: $N_CHROMS"

# ── 2. Generate per-chromosome BED files ─────────────────────────────────────
echo "[PREP] Generating per-chromosome BED files …"

CHROM_LIST="${ROOT}/meta/chromosome_list.txt"
> "$CHROM_LIST"

while IFS=$'\t' read -r CHROM LENGTH OFFSET LINEBASES LINEWIDTH; do
    BED_FILE="${ROOT}/beds/${CHROM}.bed"
    echo -e "${CHROM}\t0\t${LENGTH}" > "$BED_FILE"
    echo "$CHROM" >> "$CHROM_LIST"
done < "$REF_FAI"

echo "[PREP] BED files written to: ${ROOT}/beds/"
echo "[PREP] Chromosome list: $CHROM_LIST  ($N_CHROMS chromosomes)"

# ── 3. Build Clair3 sample manifest ──────────────────────────────────────────
echo "[PREP] Building Clair3 sample manifest …"

[[ -s "$BAM_MANIFEST" ]] || { echo "[ERROR] Missing BAM manifest: $BAM_MANIFEST" >&2; exit 1; }

CLAIR3_MANIFEST="${ROOT}/meta/clair3_sample_manifest.tsv"
echo -e "Sample\tBAM" > "$CLAIR3_MANIFEST"

N_SAMPLES=0
N_MISSING=0
N_NO_INDEX=0

# Skip header, extract sample (col1) and filtered BAM (col3)
tail -n +2 "$BAM_MANIFEST" | while IFS=$'\t' read -r SAMPLE BAM_MM2 BAM_FILT REST; do
    if [[ -z "$BAM_FILT" || "$BAM_FILT" == "NA" ]]; then
        echo "[WARN] No filtered BAM for $SAMPLE" >&2
        continue
    fi

    if [[ ! -s "$BAM_FILT" ]]; then
        echo "[WARN] BAM not found: $BAM_FILT ($SAMPLE)" >&2
        N_MISSING=$(( N_MISSING + 1 ))
        continue
    fi

    # Check index
    if [[ ! -f "${BAM_FILT}.bai" && ! -f "${BAM_FILT%.bam}.bai" ]]; then
        if [[ "$DRY_RUN" -eq 0 ]]; then
            echo "[INFO] Indexing: $BAM_FILT"
            samtools index "$BAM_FILT"
        else
            echo "[DRY] Would index: $BAM_FILT"
            N_NO_INDEX=$(( N_NO_INDEX + 1 ))
        fi
    fi

    echo -e "${SAMPLE}\t${BAM_FILT}" >> "$CLAIR3_MANIFEST"
    N_SAMPLES=$(( N_SAMPLES + 1 ))
done

# Count actual samples written (excluding header)
ACTUAL_SAMPLES=$(tail -n +2 "$CLAIR3_MANIFEST" | wc -l)

echo "[PREP] Clair3 manifest: $CLAIR3_MANIFEST"
echo "[PREP] Samples with valid BAMs: $ACTUAL_SAMPLES"

# ── 4. Generate array dispatch metadata ──────────────────────────────────────
echo "[PREP] Generating dispatch metadata …"

# Sample count for SLURM array
echo "$ACTUAL_SAMPLES" > "${ROOT}/meta/n_samples.txt"
echo "$N_CHROMS" > "${ROOT}/meta/n_chromosomes.txt"

# Combined dispatch table: sample × chromosome
DISPATCH="${ROOT}/meta/dispatch_table.tsv"
echo -e "TASK_ID\tSAMPLE\tBAM\tCHROM\tBED" > "$DISPATCH"

TASK=0
tail -n +2 "$CLAIR3_MANIFEST" | while IFS=$'\t' read -r SAMPLE BAM; do
    while read -r CHROM; do
        TASK=$(( TASK + 1 ))
        echo -e "${TASK}\t${SAMPLE}\t${BAM}\t${CHROM}\t${ROOT}/beds/${CHROM}.bed"
    done < "$CHROM_LIST"
done >> "$DISPATCH"

TOTAL_TASKS=$(tail -n +2 "$DISPATCH" | wc -l)
echo "[PREP] Dispatch table: $DISPATCH"
echo "[PREP] Total tasks (sample × chromosome): $TOTAL_TASKS"

# ── 5. Summary ───────────────────────────────────────────────────────────────
echo ""
echo "════════════════════════════════════════════════════════"
echo "[PREP] SUMMARY"
echo "  Reference:      $REF"
echo "  Chromosomes:    $N_CHROMS"
echo "  Samples:        $ACTUAL_SAMPLES"
echo "  Total tasks:    $TOTAL_TASKS"
echo ""
echo "  BED dir:        ${ROOT}/beds/"
echo "  Manifest:       $CLAIR3_MANIFEST"
echo "  Chromosome list: $CHROM_LIST"
echo "  Dispatch table: $DISPATCH"
echo ""
echo "  Next step:"
echo "    # Per-chromosome mode (recommended for large cohorts):"
echo "    for CHR in \$(cat ${CHROM_LIST}); do"
echo "      sbatch --array=1-${ACTUAL_SAMPLES} 01_run_clair3_discovery.sh --chrom \$CHR"
echo "    done"
echo ""
echo "    # Or single-dispatch mode (all at once):"
echo "    sbatch --array=1-${TOTAL_TASKS} 01_run_clair3_discovery_allchr.sh"
echo "════════════════════════════════════════════════════════"
