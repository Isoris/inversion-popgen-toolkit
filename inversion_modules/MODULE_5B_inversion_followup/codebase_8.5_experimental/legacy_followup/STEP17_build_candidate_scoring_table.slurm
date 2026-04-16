#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH -t 6:00:00
#SBATCH -J inv_step17
#SBATCH -o logs/inv_step17.%j.out
#SBATCH -e logs/inv_step17.%j.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

RSCRIPT_BIN="/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript"

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
INVDIR="${BASE}/inversion_localpca_v7"
SCRIPTDIR="${BASE}/inversion_codebase_v8.4/MODULE_5B_Inversion_Followup/legacy_followup"

CAND_TABLE="${INVDIR}/06_mds_candidates/inversion_localpca.candidate_regions.tsv.gz"
OUTDIR="${INVDIR}/13_scoring_table"

mkdir -p "${OUTDIR}" logs

# Prepare K4 ancestry table if not already present
ANCESTRY_FULL="${BASE}/popstruct_thin/05_ngsadmix_global/sample_main_ancestry_by_K.tsv"
ANCESTRY_K4="${INVDIR}/sample_main_ancestry_K4.tsv"

if [[ -s "${ANCESTRY_FULL}" && ! -s "${ANCESTRY_K4}" ]]; then
  awk -F'\t' 'NR==1 || $1=="4"' "${ANCESTRY_FULL}" > "${ANCESTRY_K4}"
  echo "[INFO] Created K4 ancestry table: ${ANCESTRY_K4}"
fi

ANCESTRY_ARG="NONE"
[[ -s "${ANCESTRY_K4}" ]] && ANCESTRY_ARG="${ANCESTRY_K4}"

"${RSCRIPT_BIN}" "${SCRIPTDIR}/STEP17_build_candidate_scoring_table.R" \
  "${INVDIR}/11_candidate_group_tables/inversion_candidate_sample_groups.tsv.gz" \
  "${INVDIR}/12_block_coherence/inversion_block_coherence_per_sample.tsv.gz" \
  "${CAND_TABLE}" \
  "${OUTDIR}" \
  "${ANCESTRY_ARG}" \
  NONE

echo "[DONE] STEP17"
