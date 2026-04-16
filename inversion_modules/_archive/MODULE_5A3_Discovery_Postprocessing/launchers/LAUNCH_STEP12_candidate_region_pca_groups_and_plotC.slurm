#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH -t 1-00:00:00
#SBATCH -J inv_step12
#SBATCH -o logs/inv_step12.%A_%a.out
#SBATCH -e logs/inv_step12.%A_%a.err

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

CAND_TABLE="${INVDIR}/06_mds_candidates/inversion_localpca.candidate_regions.tsv.gz"
SAMPLE_ORDER="${BASE}/het_roh/01_inputs_check/samples.ind"

# safer candidate->chrom lookup (no pipefail/SIGPIPE issue)
CHROM="$(python3 - <<PY
import gzip, csv
cand_id = str(${CANDIDATE_ID})
path = r"${CAND_TABLE}"
with gzip.open(path, "rt") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        if str(row["candidate_id"]) == cand_id:
            print(row["chrom"])
            break
PY
)"

if [[ -z "${CHROM}" ]]; then
  echo "[ERROR] Could not determine chromosome for candidate_id=${CANDIDATE_ID}" >&2
  exit 1
fi

SITES="${INVDIR}/04_dosage_by_chr/${CHROM}.sites.tsv.gz"
DOSAGE="${INVDIR}/04_dosage_by_chr/${CHROM}.dosage.tsv.gz"
TP_WINDOWS="${INVDIR}/07_het_overlap/het_bridge_${CHROM}.sample_tP_windows.tsv.gz"

OUTDIR="${INVDIR}/08_regional_pca"
mkdir -p "${OUTDIR}" logs

OUTPREFIX="${OUTDIR}/STEP12_${CHROM}"

[[ -x "${RSCRIPT_BIN}" ]] || { echo "[ERROR] Rscript not executable: ${RSCRIPT_BIN}" >&2; exit 1; }
[[ -s "${CAND_TABLE}" ]] || { echo "[ERROR] Missing ${CAND_TABLE}" >&2; exit 1; }
[[ -s "${SAMPLE_ORDER}" ]] || { echo "[ERROR] Missing ${SAMPLE_ORDER}" >&2; exit 1; }
[[ -s "${SITES}" ]] || { echo "[ERROR] Missing ${SITES}" >&2; exit 1; }
[[ -s "${DOSAGE}" ]] || { echo "[ERROR] Missing ${DOSAGE}" >&2; exit 1; }
[[ -s "${TP_WINDOWS}" ]] || { echo "[ERROR] Missing ${TP_WINDOWS}" >&2; exit 1; }

# skip if already done
if [[ -s "${OUTPREFIX}.candidate_${CANDIDATE_ID}.regional_pca_samples.tsv.gz" \
   && -s "${OUTPREFIX}.candidate_${CANDIDATE_ID}.plotC_summary.tsv.gz" \
   && -s "${OUTPREFIX}.candidate_${CANDIDATE_ID}.panelA_regional_PCA.pdf" \
   && -s "${OUTPREFIX}.candidate_${CANDIDATE_ID}.panelB_regional_PCA_het.pdf" \
   && -s "${OUTPREFIX}.candidate_${CANDIDATE_ID}.panelC_chrom_tP_by_group.pdf" ]]; then
  echo "[INFO] STEP12 outputs already exist for candidate ${CANDIDATE_ID}, skipping"
  exit 0
fi

echo "[INFO] CANDIDATE_ID=${CANDIDATE_ID}"
echo "[INFO] CHROM=${CHROM}"
echo "[INFO] CAND_TABLE=${CAND_TABLE}"
echo "[INFO] SITES=${SITES}"
echo "[INFO] DOSAGE=${DOSAGE}"
echo "[INFO] TP_WINDOWS=${TP_WINDOWS}"
echo "[INFO] SAMPLE_ORDER=${SAMPLE_ORDER}"
echo "[INFO] OUTPREFIX=${OUTPREFIX}"

"${RSCRIPT_BIN}" "${SCRIPTDIR}/steps/STEP12_candidate_region_pca_groups_and_plotC.R" \
  "${CAND_TABLE}" \
  "${CANDIDATE_ID}" \
  "${SITES}" \
  "${DOSAGE}" \
  "${TP_WINDOWS}" \
  "${SAMPLE_ORDER}" \
  "${OUTPREFIX}" \
  3

echo "[DONE] candidate ${CANDIDATE_ID}"
