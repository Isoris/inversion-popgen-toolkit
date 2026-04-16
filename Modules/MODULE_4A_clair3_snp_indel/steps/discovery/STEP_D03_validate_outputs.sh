#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# 02_validate_discovery_outputs.sh
#
# Run AFTER all Clair3 discovery jobs complete.
# Validates that every expected VCF/GVCF exists, is indexed,
# and has a non-zero number of records.
#
# Also audits phase output:
#   - PS header present?
#   - number of phased GT records?
#
# Usage:
#   bash 02_validate_discovery_outputs.sh [--chrom C_gar_LG28]
#   bash 02_validate_discovery_outputs.sh --all_chroms
# ============================================================

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
ROOT="${BASE}/MODULE_4A_SNP_INDEL50_Clair3"
MANIFEST="${ROOT}/meta/clair3_sample_manifest.tsv"
CHROM_LIST="${ROOT}/meta/chromosome_list.txt"

MODE="single"
TARGET_CHROM=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --chrom)      TARGET_CHROM="$2"; shift 2 ;;
        --all_chroms) MODE="all"; shift ;;
        *) echo "Unknown: $1" >&2; exit 1 ;;
    esac
done

[[ -s "$MANIFEST" ]] || { echo "[ERROR] Missing manifest: $MANIFEST" >&2; exit 1; }

# Build chromosome list
if [[ "$MODE" == "all" ]]; then
    CHROMS=($(cat "$CHROM_LIST"))
elif [[ -n "$TARGET_CHROM" ]]; then
    CHROMS=("$TARGET_CHROM")
else
    echo "[ERROR] Specify --chrom or --all_chroms" >&2
    exit 1
fi

REPORT="${ROOT}/meta/discovery_validation_report.tsv"
echo -e "SAMPLE\tCHROM\tVCF_EXISTS\tVCF_INDEXED\tVCF_RECORDS\tGVCF_EXISTS\tHAS_PS_HEADER\tN_PHASED_GT\tSTATUS" > "$REPORT"

for CHROM in "${CHROMS[@]}"; do
    tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r SAMPLE BAM; do
        VCF="${ROOT}/vcf/${CHROM}/${SAMPLE}.${CHROM}.clair3.discovery.vcf.gz"
        GVCF="${ROOT}/gvcf/${CHROM}/${SAMPLE}.${CHROM}.clair3.discovery.gvcf.gz"

        VCF_EXISTS="no"
        VCF_INDEXED="no"
        VCF_RECORDS=0
        GVCF_EXISTS="no"
        HAS_PS_HEADER="no"
        N_PHASED_GT=0
        STATUS="MISSING"

        if [[ -s "$VCF" ]]; then
            VCF_EXISTS="yes"

            if [[ -f "${VCF}.csi" || -f "${VCF}.tbi" ]]; then
                VCF_INDEXED="yes"
            fi

            VCF_RECORDS=$(bcftools view -H "$VCF" 2>/dev/null | wc -l || echo 0)

            if bcftools view -h "$VCF" | grep -q '^##FORMAT=<ID=PS,'; then
                HAS_PS_HEADER="yes"
            fi

            N_PHASED_GT=$(bcftools view -H "$VCF" | awk '$10 ~ /\|/' | wc -l)

            if [[ "$VCF_RECORDS" -gt 0 ]]; then
                STATUS="OK"
            else
                STATUS="EMPTY"
            fi
        fi

        if [[ -s "$GVCF" ]]; then
            GVCF_EXISTS="yes"
        fi

        echo -e "${SAMPLE}\t${CHROM}\t${VCF_EXISTS}\t${VCF_INDEXED}\t${VCF_RECORDS}\t${GVCF_EXISTS}\t${HAS_PS_HEADER}\t${N_PHASED_GT}\t${STATUS}"
    done >> "$REPORT"
done

echo ""
echo "════════════════════════════════════════════════════════"
echo "[VALIDATE] Discovery output validation"
echo "  Chromosomes checked: ${#CHROMS[@]}"
echo "  Report: $REPORT"
echo ""

N_OK=$(awk -F'\t' 'NR>1 && $9=="OK"' "$REPORT" | wc -l)
N_MISSING=$(awk -F'\t' 'NR>1 && $9=="MISSING"' "$REPORT" | wc -l)
N_EMPTY=$(awk -F'\t' 'NR>1 && $9=="EMPTY"' "$REPORT" | wc -l)
N_WITH_PS=$(awk -F'\t' 'NR>1 && $7=="yes"' "$REPORT" | wc -l)
N_WITH_PHASED=$(awk -F'\t' 'NR>1 && $8>0' "$REPORT" | wc -l)
N_TOTAL=$(tail -n +2 "$REPORT" | wc -l)

echo "  Total expected:      $N_TOTAL"
echo "  OK:                  $N_OK"
echo "  Missing:             $N_MISSING"
echo "  Empty:               $N_EMPTY"
echo "  With PS header:      $N_WITH_PS"
echo "  With phased GTs:     $N_WITH_PHASED"

if [[ "$N_MISSING" -gt 0 ]]; then
    echo ""
    echo "  MISSING samples:"
    awk -F'\t' 'NR>1 && $9=="MISSING" {print "    " $1 " / " $2}' "$REPORT" | head -20
    REMAINING=$(( N_MISSING - 20 ))
    [[ "$REMAINING" -gt 0 ]] && echo "    ... and $REMAINING more"
fi

if [[ "$N_MISSING" -eq 0 && "$N_EMPTY" -eq 0 ]]; then
    echo ""
    echo "  ✓ All discovery outputs valid. Ready for post-processing."
else
    echo ""
    echo "  ✗ Some outputs missing or empty. Check failed jobs in logs/."
fi

echo "════════════════════════════════════════════════════════"
