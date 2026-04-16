#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH -t 12:00:00
#SBATCH -J inv_step13
#SBATCH -o logs/inv_step13.%A_%a.out
#SBATCH -e logs/inv_step13.%A_%a.err

set -euo pipefail

source ~/.bashrc
mamba activate assembly

RSCRIPT_BIN="/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript"

CANDIDATE_LIST="${1:?ERROR: provide candidate_ids.list as argv1}"

IDX="${SLURM_ARRAY_TASK_ID}"
LINE=$((IDX + 1))
CANDIDATE_ID="$(sed -n "${LINE}p" "${CANDIDATE_LIST}")"

if [[ -z "${CANDIDATE_ID}" ]]; then
  echo "[ERROR] Could not read candidate ID from ${CANDIDATE_LIST} line ${LINE}" >&2
  exit 1
fi

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
INVDIR="${BASE}/inversion_localpca_v7"
SCRIPTDIR="${BASE}/inversion_codebase_v8.4/MODULE_5A3_Discovery_Postprocessing"
STEP12DIR="${INVDIR}/08_regional_pca"
OUTDIR="${INVDIR}/09_combined_plots"
mkdir -p "${OUTDIR}" logs

# get chromosome for this candidate
CAND_TABLE="${INVDIR}/06_mds_candidates/inversion_localpca.candidate_regions.tsv.gz"
CHROM="$(zcat -f "${CAND_TABLE}" | awk -F'\t' -v cid="${CANDIDATE_ID}" '
NR==1{
  for(i=1;i<=NF;i++){
    if($i=="candidate_id") cidcol=i
    if($i=="chrom") chrcol=i
  }
  next
}
$cidcol==cid {print $chrcol; exit}
')"

[[ -n "${CHROM}" ]] || { echo "[ERROR] Could not determine chromosome for candidate ${CANDIDATE_ID}" >&2; exit 1; }

STEP12_PREFIX="${STEP12DIR}/STEP12_${CHROM}"
OUTPREFIX="${OUTDIR}/STEP13_${CHROM}"

[[ -s "${STEP12_PREFIX}.candidate_${CANDIDATE_ID}.regional_pca_samples.tsv.gz" ]] || { echo "[ERROR] Missing STEP12 output for candidate ${CANDIDATE_ID}" >&2; exit 1; }

if [[ -s "${OUTPREFIX}.candidate_${CANDIDATE_ID}.combined_ABC.pdf" \
   && -s "${OUTPREFIX}.candidate_${CANDIDATE_ID}.combined_ABC.png" ]]; then
  echo "[INFO] STEP13 outputs already exist for candidate ${CANDIDATE_ID}, skipping"
  exit 0
fi

echo "[INFO] CANDIDATE_ID=${CANDIDATE_ID}"
echo "[INFO] CHROM=${CHROM}"
echo "[INFO] STEP12_PREFIX=${STEP12_PREFIX}"
echo "[INFO] OUTPREFIX=${OUTPREFIX}"

"${RSCRIPT_BIN}" "${SCRIPTDIR}/steps/STEP13_make_combined_panels_ABC.R" \
  "${STEP12_PREFIX}" \
  "${CANDIDATE_ID}" \
  "${OUTPREFIX}" \
  8 \
  13

echo "[DONE] candidate ${CANDIDATE_ID}"
