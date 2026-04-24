#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH -t 12:00:00
#SBATCH -J inv_step18
#SBATCH -o logs/inv_step18.%A_%a.out
#SBATCH -e logs/inv_step18.%A_%a.err

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
SCRIPTDIR="${BASE}/inversion_codebase_v8.4/MODULE_5B_Inversion_Followup/legacy_followup"
OUTDIR="${INVDIR}/15_composite_figures"

mkdir -p "${OUTDIR}" logs

# Get chromosome for this candidate
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

# Skip if output already exists
if [[ -s "${OUTDIR}/composite_candidate_${CANDIDATE_ID}_${CHROM}.png" ]]; then
  echo "[INFO] Composite figure already exists for candidate ${CANDIDATE_ID}, skipping"
  exit 0
fi

# LD file (optional — will render placeholder if missing)
LD_FILE="${INVDIR}/08_ld_pruned81/${CHROM}.pruned81.pairs_r2.tsv"
LD_ARG="NONE"
[[ -s "${LD_FILE}" ]] && LD_ARG="${LD_FILE}"

echo "[INFO] CANDIDATE_ID=${CANDIDATE_ID} CHROM=${CHROM}"

"${RSCRIPT_BIN}" "${SCRIPTDIR}/STEP18_composite_candidate_figure.R" \
  "${CANDIDATE_ID}" \
  "${CAND_TABLE}" \
  "${INVDIR}/08_regional_pca" \
  "${INVDIR}/06_mds_candidates/inversion_localpca.mds.rds" \
  "${INVDIR}/04_dosage_by_chr" \
  "${OUTDIR}" \
  "${LD_ARG}" \
  NONE \
  250000

echo "[DONE] candidate ${CANDIDATE_ID}"
