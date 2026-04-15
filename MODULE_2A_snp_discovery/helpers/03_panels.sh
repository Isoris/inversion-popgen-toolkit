#!/usr/bin/env bash
###############################################################################
# helpers/03_panels.sh — Merge biSNP sites + distance-thin + ANGSD-index
# Combines old: S07 (200/500/1000), S08 (5k/10k/25k), S09 (summary)
# Called by: run_step1.sh build_panels
###############################################################################
set -euo pipefail
source "$(dirname "$0")/../config.sh"

MAF_DIR="${GLOBAL_DIR}/02_snps"
SITESDIR="${THIN_DIR}/03_sites"
mkdir -p "${SITESDIR}"

MAF_GLOB="${MAF_DIR}/catfish.*.mafs.gz"
NFILES=$(ls -1 ${MAF_GLOB} 2>/dev/null | wc -l)
(( NFILES > 0 )) || { echo "[ERROR] No .mafs.gz files in ${MAF_DIR}" >&2; exit 1; }
echo "[$(timestamp)] Found ${NFILES} mafs.gz files"

# ---- 1) Merge: 4-col (chr, pos, major, minor) ----
RAW_MAJMIN="${SITESDIR}/sites.raw.ALL.majmin.tsv"
if [[ ! -s "$RAW_MAJMIN" ]]; then
  echo "[$(timestamp)] Merging to 4-col majmin..."
  zcat ${MAF_GLOB} \
    | awk 'BEGIN{OFS="\t"} $1!="chromo" {print $1,$2,$3,$4}' \
    | awk 'BEGIN{OFS="\t"} length($3)==1 && length($4)==1 && $3!="N" && $3!="n" && $4!="N" && $4!="n" && $3!=$4 {print}' \
    | sort -k1,1 -k2,2n > "${RAW_MAJMIN}"
fi
echo "[$(timestamp)] Raw biSNPs: $(wc -l < "${RAW_MAJMIN}")"

# ---- 2) Thin function ----
thin_majmin() {
  local W="$1"
  local OUT="${SITESDIR}/sites.thin_${W}.majmin.tsv"
  if [[ -s "$OUT" ]]; then echo "[SKIP] thin_${W} exists"; return; fi
  awk -v W="${W}" 'BEGIN{OFS="\t"} {c=$1;p=$2;if(c!=pc){pc=c;last=-1e18} if(p-last>=W){print;last=p}}' \
    "${RAW_MAJMIN}" > "${OUT}"
  angsd sites index "${OUT}"
  echo "[$(timestamp)] thin_${W}: $(wc -l < "${OUT}") sites"
}

thin_2col() {
  local W="$1"
  local OUT="${SITESDIR}/sites.thin_${W}.tsv"
  if [[ -s "$OUT" ]]; then echo "[SKIP] thin_${W} exists"; return; fi
  awk -v W="${W}" 'BEGIN{OFS="\t"} {c=$1;p=$2;if(c!=pc){pc=c;last=-1e18} if(p-last>=W){print $1,$2;last=p}}' \
    "${RAW_MAJMIN}" > "${OUT}"
  angsd sites index "${OUT}"
  echo "[$(timestamp)] thin_${W}: $(wc -l < "${OUT}") sites"
}

for W in "${THIN_FINE[@]}"; do thin_majmin "$W"; done
for W in "${THIN_BROAD[@]}"; do thin_2col "$W"; done

# ---- 3) Summary for manuscript ----
echo ""
echo "=== Site panel summary ==="
echo -e "panel\tn_sites"
for W in "${THIN_ALL[@]}"; do
  for suf in majmin.tsv tsv; do
    F="${SITESDIR}/sites.thin_${W}.${suf}"
    [[ -s "$F" ]] && { echo -e "thin_${W}\t$(wc -l < "$F")"; break; }
  done
done

echo "[$(timestamp)] [DONE] 03_panels"
