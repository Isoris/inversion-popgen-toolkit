#!/usr/bin/env bash
# =============================================================================
# harvest_large_het_dels.sh
#
# Extracts large (≥5kb) het-DELs from the DELLY DEL catalog for use as
# annotation in Snake 3 GHSL and inversion breakpoint analysis.
#
# Outputs:
#   large_het_dels.bed       — BED with chrom, start, end, length, carrier counts
#   large_het_dels_detail.tsv — same + per-sample genotypes
#
# Usage:
#   bash harvest_large_het_dels.sh [min_size_bp]
#
# Sources 00_inversion_config.sh for paths. Default min_size = 5000 bp.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source config
CONFIG="${SCRIPT_DIR}/../00_inversion_config.sh"
if [[ -f "${CONFIG}" ]]; then
  set -a; source "${CONFIG}"; set +a
else
  echo "[ERROR] Missing config: ${CONFIG}" >&2
  exit 1
fi

MIN_SIZE="${1:-5000}"

# ── Locate DELLY DEL catalog ──
DELLY_DEL_DIR="${BASE}/MODULE_4B_DEL_Delly"
CATALOG_226="${DELLY_DEL_DIR}/07_final_catalogs/catalog_226.DEL.vcf.gz"
CATALOG_81="${DELLY_DEL_DIR}/07_final_catalogs/catalog_81.DEL.germline.PASS.vcf.gz"

# Use 226 if available, fall back to 81
if [[ -f "${CATALOG_226}" ]]; then
  VCF="${CATALOG_226}"
  COHORT="226"
elif [[ -f "${CATALOG_81}" ]]; then
  VCF="${CATALOG_81}"
  COHORT="81"
else
  echo "[ERROR] No DELLY DEL catalog found in ${DELLY_DEL_DIR}/07_final_catalogs/" >&2
  echo "[ERROR] Expected: catalog_226.DEL.vcf.gz or catalog_81.DEL.germline.PASS.vcf.gz" >&2
  exit 1
fi

OUTDIR="${MDS_DIR}/snake_regions_multiscale"
mkdir -p "${OUTDIR}"

OUTBED="${OUTDIR}/large_het_dels_${MIN_SIZE}bp.bed"
OUTDETAIL="${OUTDIR}/large_het_dels_${MIN_SIZE}bp_detail.tsv"
OUTSUMMARY="${OUTDIR}/large_het_dels_${MIN_SIZE}bp_summary.tsv"

echo "[harvest-dels] DELLY catalog: ${VCF} (cohort=${COHORT})"
echo "[harvest-dels] Min size: ${MIN_SIZE} bp"
echo "[harvest-dels] Output: ${OUTBED}"

# ── Get sample names from VCF header ──
SAMPLE_NAMES=$(bcftools query -l "${VCF}" | tr '\n' '\t')

# ── Extract large DELs with per-sample genotypes ──
# Header for detail file
echo -e "chrom\tstart\tend\tsvlen\tn_het\tn_hom_alt\tn_carriers\tcarrier_freq\tcarrier_samples" > "${OUTDETAIL}"

# Detect whether SVLEN is populated or needs END-POS fallback
HAS_SVLEN=$(bcftools query -f '%INFO/SVLEN\n' "${VCF}" | head -20 | grep -cv '^\.$')
if [[ "${HAS_SVLEN}" -gt 5 ]]; then
  echo "[harvest-dels] SVLEN field populated — using SVLEN for length"
  LEN_MODE="svlen"
else
  echo "[harvest-dels] SVLEN is missing/empty — using END-POS for length"
  LEN_MODE="endpos"
fi

bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVLEN\t[%SAMPLE=%GT\t]\n' "${VCF}" | \
awk -F'\t' -v min_size="${MIN_SIZE}" -v len_mode="${LEN_MODE}" '
{
  # Compute length: prefer SVLEN if available, else END-POS
  if (len_mode == "svlen" && $4 != "." && $4 != "") {
    svlen = $4 + 0
    if (svlen < 0) svlen = -svlen
  } else {
    # Fallback: END - POS
    svlen = ($3 + 0) - ($2 + 0)
    if (svlen < 0) svlen = -svlen
  }
  if (svlen < min_size) next

  n_het = 0; n_hom = 0; n_total = 0; carriers = ""
  for (i = 5; i <= NF; i++) {
    if ($i == "" || $i == ".") continue
    # Format: SAMPLE=GT
    split($i, parts, "=")
    samp = parts[1]; gt = parts[2]
    # Skip missing genotypes
    if (gt == "./." || gt == ".|." || gt == ".") continue
    n_total++
    if (gt == "0/1" || gt == "0|1" || gt == "1|0") {
      n_het++
      carriers = carriers samp ","
    } else if (gt == "1/1" || gt == "1|1") {
      n_hom++
      carriers = carriers samp ","
    }
  }

  n_carriers = n_het + n_hom
  if (n_carriers == 0) next

  freq = (n_total > 0) ? n_carriers / n_total : 0
  # Remove trailing comma
  gsub(/,$/, "", carriers)

  # BED line (to stdout)
  printf "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.4f\n", $1, $2, $3, svlen, n_het, n_hom, n_carriers, freq

  # Detail line (to stderr for capture)
  printf "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%s\n", $1, $2, $3, svlen, n_het, n_hom, n_carriers, freq, carriers > "/dev/stderr"
}' > "${OUTBED}" 2>> "${OUTDETAIL}"

N_TOTAL=$(wc -l < "${OUTBED}")

# ── Per-chromosome summary ──
echo -e "chrom\tn_large_dels\tmin_size\tmax_size\tmedian_size\tmean_carriers" > "${OUTSUMMARY}"
awk -F'\t' '{
  chr[$1]++
  size[$1] = size[$1] " " $4
  carriers[$1] += $7
}
END {
  for (c in chr) {
    n = chr[c]
    mean_c = carriers[c] / n
    # Get sizes for this chr
    split(size[c], s, " ")
    min_s = 999999999; max_s = 0
    for (i in s) {
      if (s[i]+0 > 0) {
        if (s[i]+0 < min_s) min_s = s[i]+0
        if (s[i]+0 > max_s) max_s = s[i]+0
      }
    }
    printf "%s\t%d\t%d\t%d\t.\t%.1f\n", c, n, min_s, max_s, mean_c
  }
}' "${OUTBED}" | sort -k1,1V >> "${OUTSUMMARY}"

echo ""
echo "[harvest-dels] === RESULTS ==="
echo "  Total large DELs (>=${MIN_SIZE} bp): ${N_TOTAL}"
echo "  BED:     ${OUTBED}"
echo "  Detail:  ${OUTDETAIL}"
echo "  Summary: ${OUTSUMMARY}"
echo ""

if [[ "${N_TOTAL}" -gt 0 ]]; then
  echo "[harvest-dels] Size distribution:"
  awk -F'\t' '{
    if ($4 >= 100000) bin="100kb+"
    else if ($4 >= 50000) bin="50-100kb"
    else if ($4 >= 20000) bin="20-50kb"
    else if ($4 >= 10000) bin="10-20kb"
    else bin="5-10kb"
    counts[bin]++
  }
  END {
    for (b in counts) printf "  %s: %d\n", b, counts[b]
  }' "${OUTBED}" | sort

  echo ""
  echo "[harvest-dels] Per-chromosome counts:"
  cut -f1 "${OUTBED}" | sort | uniq -c | sort -k2,2V | awk '{printf "  %s: %d\n", $2, $1}'
fi

echo ""
echo "[DONE]"
