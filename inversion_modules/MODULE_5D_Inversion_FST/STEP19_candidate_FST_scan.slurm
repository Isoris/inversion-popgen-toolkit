#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128GB
#SBATCH -t 1-00:00:00
#SBATCH -J inv_step19
#SBATCH -o logs/inv_step19.%j.out
#SBATCH -e logs/inv_step19.%j.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

RSCRIPT_BIN="/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript"

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
INVDIR="${BASE}/inversion_localpca_v7"
SCRIPTDIR="${BASE}/inversion_codebase_v8.4/MODULE_5D_Inversion_FST"

CAND_TABLE="${INVDIR}/06_mds_candidates/inversion_localpca.candidate_regions.tsv.gz"
OUTDIR="${INVDIR}/14_fst_scans"

mkdir -p "${OUTDIR}" logs

"${RSCRIPT_BIN}" "${SCRIPTDIR}/STEP19_candidate_FST_scan.R" \
  "${INVDIR}/13_scoring_table/inversion_candidate_scoring_table.tsv.gz" \
  "${CAND_TABLE}" \
  "${INVDIR}/04_dosage_by_chr" \
  "${OUTDIR}" \
  all \
  0.10 \
  50

echo "[DONE] STEP19"
