#!/usr/bin/env bash
set -euo pipefail

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
INVDIR="${BASE}/inversion_localpca_v7"
OUTDIR="${INVDIR}/supplementary_tables"
mkdir -p "${OUTDIR}"

COUNTS="${OUTDIR}/supp_candidate_counts_per_chr.tsv"
MANIFEST="${OUTDIR}/supp_inversion_step_manifest.tsv"
OUTTXT="${OUTDIR}/supp_inversion_workflow_notes.txt"

{
  echo "Inversion local-PCA supplementary notes"
  echo "======================================="
  echo
  echo "Candidate counts per chromosome:"
  cat "${COUNTS}"
  echo
  echo "Step manifest:"
  cat "${MANIFEST}"
  echo
  echo "Key files:"
  echo "  candidate regions: ${INVDIR}/06_mds_candidates/inversion_localpca.candidate_regions.tsv.gz"
  echo "  window MDS:        ${INVDIR}/06_mds_candidates/inversion_localpca.window_mds.tsv.gz"
  echo "  local PCA:         ${INVDIR}/05_local_pca/"
  echo "  theta overlap:     ${INVDIR}/07_het_overlap/"
} > "${OUTTXT}"

echo "[DONE] ${OUTTXT}"
