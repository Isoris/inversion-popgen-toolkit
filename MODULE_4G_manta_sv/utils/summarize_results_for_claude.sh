#!/usr/bin/env bash
# =============================================================================
# summarize_manta_results_for_claude.sh
# Quick summary of Manta pipeline outputs for sharing with Claude
# =============================================================================
set -euo pipefail
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${SCRIPT_DIR}/00_manta_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config" >&2; exit 1; }
set -a; source "${CONFIG}"; set +a

echo "=== MODULE_4G_ALL_Manta results summary ==="
echo "Date: $(date '+%F %T')"
echo ""

echo "--- SV counts (PASS) ---"
for cohort in catalog_226 catalog_81; do
  echo "  ${cohort}:"
  for svtype in DEL DUP INV BND INS_small INS_large; do
    vcf="${DIR_FINAL}/${cohort}.${svtype}.PASS.vcf.gz"
    if [[ -f "${vcf}" ]]; then
      n=$(bcftools view -H "${vcf}" 2>/dev/null | wc -l)
      printf "    %-12s %s\n" "${svtype}" "${n}"
    fi
  done
done

echo ""
echo "--- File tree ---"
find "${OUTDIR}" -maxdepth 2 -type f \( -name '*.vcf.gz' -o -name '*.tsv' -o -name '*.bed' -o -name '*.txt' -o -name '*.pdf' -o -name '*.png' \) \
  | sort | head -80

echo ""
echo "--- Summary report ---"
[[ -f "${DIR_SUMMARY}/manta_SV_summary_report.txt" ]] && cat "${DIR_SUMMARY}/manta_SV_summary_report.txt"
