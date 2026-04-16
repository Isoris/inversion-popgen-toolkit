#!/usr/bin/env bash
set -euo pipefail

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
INVDIR="${BASE}/inversion_localpca_v7"
OUTDIR="${INVDIR}/supplementary_tables"
mkdir -p "${OUTDIR}"

CAND="${INVDIR}/06_mds_candidates/inversion_localpca.candidate_regions.tsv.gz"
OUT="${OUTDIR}/supp_candidate_counts_per_chr.tsv"

zcat -f "${CAND}" \
| tail -n +2 \
| cut -f3 \
| sort \
| uniq -c \
| awk 'BEGIN{OFS="\t"; print "n_candidates","chrom"} {print $1,$2}' \
| sort -k1,1nr -k2,2 \
> "${OUT}"

echo "[DONE] ${OUT}"
cat "${OUT}"
