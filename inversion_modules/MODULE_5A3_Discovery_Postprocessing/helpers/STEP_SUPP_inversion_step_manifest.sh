#!/usr/bin/env bash
set -euo pipefail

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
INVDIR="${BASE}/inversion_localpca_v7"
OUTDIR="${INVDIR}/supplementary_tables"
mkdir -p "${OUTDIR}"

OUT="${OUTDIR}/supp_inversion_step_manifest.tsv"

n_step08_sites=$(ls "${INVDIR}"/04_dosage_by_chr/C_gar_LG*.sites.tsv.gz 2>/dev/null | wc -l)
n_step08_dosage=$(ls "${INVDIR}"/04_dosage_by_chr/C_gar_LG*.dosage.tsv.gz 2>/dev/null | wc -l)

n_step09_meta=$(ls "${INVDIR}"/05_local_pca/STEP09_C_gar_LG*.window_meta.tsv.gz 2>/dev/null | wc -l)
n_step09_rds=$(ls "${INVDIR}"/05_local_pca/STEP09_C_gar_LG*.window_pca.rds 2>/dev/null | wc -l)
n_step09_tsv=$(ls "${INVDIR}"/05_local_pca/STEP09_C_gar_LG*.window_pca.tsv.gz 2>/dev/null | wc -l)

step10_mds="no"
step10_cand="no"
[[ -s "${INVDIR}/06_mds_candidates/inversion_localpca.window_mds.tsv.gz" ]] && step10_mds="yes"
[[ -s "${INVDIR}/06_mds_candidates/inversion_localpca.candidate_regions.tsv.gz" ]] && step10_cand="yes"

n_step11_sum=$(ls "${INVDIR}"/07_het_overlap/inversion_localpca.C_gar_LG*.candidate_theta_summary.tsv.gz 2>/dev/null | wc -l)
n_step11_det=$(ls "${INVDIR}"/07_het_overlap/inversion_localpca.C_gar_LG*.candidate_theta_detail.tsv.gz 2>/dev/null | wc -l)

{
  echo -e "step\tmetric\tvalue"
  echo -e "STEP08\tchr_sites_files\t${n_step08_sites}"
  echo -e "STEP08\tchr_dosage_files\t${n_step08_dosage}"
  echo -e "STEP09\twindow_meta_files\t${n_step09_meta}"
  echo -e "STEP09\twindow_pca_rds_files\t${n_step09_rds}"
  echo -e "STEP09\twindow_pca_tsv_files\t${n_step09_tsv}"
  echo -e "STEP10\twindow_mds_exists\t${step10_mds}"
  echo -e "STEP10\tcandidate_regions_exists\t${step10_cand}"
  echo -e "STEP11\tcandidate_theta_summary_files\t${n_step11_sum}"
  echo -e "STEP11\tcandidate_theta_detail_files\t${n_step11_det}"
} > "${OUT}"

echo "[DONE] ${OUT}"
column -t -s $'\t' "${OUT}"
