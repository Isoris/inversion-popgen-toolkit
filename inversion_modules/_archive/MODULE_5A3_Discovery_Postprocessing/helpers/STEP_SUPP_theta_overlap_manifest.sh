#!/usr/bin/env bash
set -euo pipefail

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
INVDIR="${BASE}/inversion_localpca_v7"
OUTDIR="${INVDIR}/supplementary_tables"
HETDIR="${INVDIR}/07_het_overlap"
CHRLIST="${INVDIR}/chrom.list"

mkdir -p "${OUTDIR}"

OUT="${OUTDIR}/supp_theta_overlap_manifest.tsv"

{
  echo -e "chrom\thet_bridge_tsv\tpopulation_mean_tsv\tcandidate_theta_summary\tcandidate_theta_detail\tstatus"
  while read -r chrom; do
    [[ -z "${chrom}" ]] && continue

    A="${HETDIR}/het_bridge_${chrom}.sample_tP_windows.tsv.gz"
    B="${HETDIR}/population_mean_tP_windows.${chrom}.tsv.gz"
    C="${HETDIR}/inversion_localpca.${chrom}.candidate_theta_summary.tsv.gz"
    D="${HETDIR}/inversion_localpca.${chrom}.candidate_theta_detail.tsv.gz"

    SA="no"; SB="no"; SC="no"; SD="no"
    [[ -s "${A}" ]] && SA="yes"
    [[ -s "${B}" ]] && SB="yes"
    [[ -s "${C}" ]] && SC="yes"
    [[ -s "${D}" ]] && SD="yes"

    STATUS="missing"
    if [[ "${SA}" == "yes" && "${SB}" == "yes" && "${SC}" == "yes" && "${SD}" == "yes" ]]; then
      STATUS="complete"
    elif [[ "${SA}" == "yes" || "${SB}" == "yes" || "${SC}" == "yes" || "${SD}" == "yes" ]]; then
      STATUS="partial"
    fi

    echo -e "${chrom}\t${SA}\t${SB}\t${SC}\t${SD}\t${STATUS}"
  done < "${CHRLIST}"
} > "${OUT}"

echo "[DONE] ${OUT}"
column -t -s $'\t' "${OUT}"
