#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH -t 6:00:00
#SBATCH -J inv_step15
#SBATCH -o logs/inv_step15.%j.out
#SBATCH -e logs/inv_step15.%j.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

RSCRIPT_BIN="/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript"

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
INVDIR="${BASE}/inversion_localpca_v7"
SCRIPTDIR="${BASE}/inversion_codebase_v8.4/MODULE_5B_Inversion_Followup/legacy_followup"

CAND_TABLE="${INVDIR}/06_mds_candidates/inversion_localpca.candidate_regions.tsv.gz"
OUTDIR="${INVDIR}/11_candidate_group_tables"

mkdir -p "${OUTDIR}" logs

"${RSCRIPT_BIN}" "${SCRIPTDIR}/STEP15_candidate_group_table_builder.R" \
  "${INVDIR}/08_regional_pca" \
  "${CAND_TABLE}" \
  "${OUTDIR}"

echo "[DONE] STEP15"
