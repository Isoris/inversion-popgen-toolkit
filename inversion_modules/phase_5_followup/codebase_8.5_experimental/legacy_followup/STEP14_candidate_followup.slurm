#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH -t 12:00:00
#SBATCH -J inv_step14
#SBATCH -o logs/inv_step14.%A_%a.out
#SBATCH -e logs/inv_step14.%A_%a.err

set -euo pipefail

source ~/.bashrc
mamba activate assembly

RSCRIPT_BIN="/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript"

CANDIDATE_LIST="${1:?ERROR: provide candidate_ids.list as argv1}"
FAMILY_TABLE="${2:-NONE}"

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

CAND_TABLE="${INVDIR}/06_mds_candidates/inversion_localpca.candidate_regions.tsv.gz"
CHROM_SIZES="${BASE}/het_roh/01_inputs_check/chrom_sizes.tsv"

CHROM="$(python3 - <<PY
import gzip, csv
cand_id = str(${CANDIDATE_ID})
with gzip.open(r"${CAND_TABLE}", "rt") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        if str(row["candidate_id"]) == cand_id:
            print(row["chrom"])
            break
PY
)"

[[ -n "${CHROM}" ]] || { echo "[ERROR] Could not determine chromosome for candidate ${CANDIDATE_ID}" >&2; exit 1; }

STEP12_PREFIX="${INVDIR}/08_regional_pca/STEP12_${CHROM}"
PCS_FILE="${STEP12_PREFIX}.candidate_${CANDIDATE_ID}.regional_pca_samples.tsv.gz"
TP_FILE="${INVDIR}/07_het_overlap/het_bridge_${CHROM}.sample_tP_windows.tsv.gz"

OUTDIR="${INVDIR}/10_candidate_followup"
mkdir -p "${OUTDIR}" logs
OUTPREFIX="${OUTDIR}/STEP14_${CHROM}"

[[ -x "${RSCRIPT_BIN}" ]] || { echo "[ERROR] Rscript not executable: ${RSCRIPT_BIN}" >&2; exit 1; }
[[ -s "${CAND_TABLE}" ]] || { echo "[ERROR] Missing ${CAND_TABLE}" >&2; exit 1; }
[[ -s "${PCS_FILE}" ]] || { echo "[ERROR] Missing ${PCS_FILE}" >&2; exit 1; }
[[ -s "${TP_FILE}" ]] || { echo "[ERROR] Missing ${TP_FILE}" >&2; exit 1; }
[[ -s "${CHROM_SIZES}" ]] || { echo "[ERROR] Missing ${CHROM_SIZES}" >&2; exit 1; }

echo "[INFO] CANDIDATE_ID=${CANDIDATE_ID}"
echo "[INFO] CHROM=${CHROM}"
echo "[INFO] PCS_FILE=${PCS_FILE}"
echo "[INFO] TP_FILE=${TP_FILE}"
echo "[INFO] FAMILY_TABLE=${FAMILY_TABLE}"
echo "[INFO] OUTPREFIX=${OUTPREFIX}"

"${RSCRIPT_BIN}" "${SCRIPTDIR}/STEP14_candidate_followup.R" \
  "${CAND_TABLE}" \
  "${CANDIDATE_ID}" \
  "${PCS_FILE}" \
  "${TP_FILE}" \
  "${OUTPREFIX}" \
  "${CHROM_SIZES}" \
  "${FAMILY_TABLE}" \
  0.10 \
  20000

echo "[DONE] candidate ${CANDIDATE_ID}"
