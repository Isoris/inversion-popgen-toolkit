#!/usr/bin/env bash
set -euo pipefail

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
ROOT="${BASE}/MODULE_4A_SNP_INDEL50_Clair3"
META="${ROOT}/meta"
OUT="${META}/phase_audit.discovery.tsv"

TMP_FILES=$(find "$META" -maxdepth 1 -type f -name 'phase_audit.discovery.*.tmp.tsv' | sort)

if [[ -z "$TMP_FILES" ]]; then
    echo "[ERROR] No phase audit temp files found in $META" >&2
    exit 1
fi

{
    echo -e "SAMPLE\tCHROM\tFOUND_PHASED_MERGED\tVCF_SOURCE_TYPE\tHAS_PS_HEADER\tN_PHASED_GT\tN_TOTAL_RECORDS\tVCF_PATH"
    first=1
    while IFS= read -r f; do
        if [[ $first -eq 1 ]]; then
            awk 'FNR>1' "$f"
            first=0
        else
            awk 'FNR>1' "$f"
        fi
    done <<< "$TMP_FILES"
} > "$OUT"

echo "[DONE] Combined phase audit written to:"
echo "  $OUT"
echo ""

awk -F'\t' '
BEGIN{
    total=0; found_phased=0; has_ps=0; with_gt=0
}
NR>1{
    total++
    if($3==1) found_phased++
    if($5==1) has_ps++
    if($6>0) with_gt++
}
END{
    print "Summary:"
    print "  total rows:             " total
    print "  found phased_merge:     " found_phased
    print "  copied VCF has PS:      " has_ps
    print "  copied VCF has phased:  " with_gt
}' "$OUT"
