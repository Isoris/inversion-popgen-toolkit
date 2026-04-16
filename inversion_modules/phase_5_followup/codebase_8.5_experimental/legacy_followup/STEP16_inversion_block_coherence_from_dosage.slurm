#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128GB
#SBATCH -t 1-00:00:00
#SBATCH -J inv_step16
#SBATCH -o logs/inv_step16.%j.out
#SBATCH -e logs/inv_step16.%j.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

RSCRIPT_BIN="/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript"

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
INVDIR="${BASE}/inversion_localpca_v7"
SCRIPTDIR="${BASE}/inversion_codebase_v8.4/MODULE_5B_Inversion_Followup/legacy_followup"

CAND_TABLE="${INVDIR}/06_mds_candidates/inversion_localpca.candidate_regions.tsv.gz"
OUTDIR="${INVDIR}/12_block_coherence"

mkdir -p "${OUTDIR}" logs

"${RSCRIPT_BIN}" "${SCRIPTDIR}/STEP16_inversion_block_coherence_from_dosage.R" \
  "${CAND_TABLE}" \
  "${INVDIR}/08_regional_pca" \
  "${INVDIR}/04_dosage_by_chr" \
  "${OUTDIR}" \
  0.5 \
  20

echo "[DONE] STEP16"
