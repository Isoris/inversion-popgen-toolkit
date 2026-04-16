#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB
#SBATCH -t 1-00:00:00
#SBATCH -J inv_step10b11
#SBATCH -o logs/inv_step10b11.%A_%a.out
#SBATCH -e logs/inv_step10b11.%A_%a.err

set -euo pipefail

source ~/.bashrc
mamba activate assembly

RSCRIPT_BIN="/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript"

CHROM_LIST="${1:?ERROR: provide chrom list as argv1}"

IDX="${SLURM_ARRAY_TASK_ID}"
LINE=$((IDX + 1))
CHROM="$(sed -n "${LINE}p" "${CHROM_LIST}")"

if [[ -z "${CHROM}" ]]; then
  echo "[ERROR] Could not read chromosome from ${CHROM_LIST} line ${LINE}" >&2
  exit 1
fi

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
INVDIR="${BASE}/inversion_localpca_v7"
HETDIR="${BASE}/het_roh"
SCRIPTDIR="${BASE}/inversion_codebase_clean_v7_3/clean/MODULE_5A_Inversion_Discovery"

mkdir -p "${INVDIR}/07_het_overlap" "${SCRIPTDIR}/logs"

THETA_PREFIX="${INVDIR}/07_het_overlap/het_bridge_${CHROM}"
POP_TP="${INVDIR}/07_het_overlap/population_mean_tP_windows.${CHROM}.tsv.gz"
OUTPREFIX="${INVDIR}/07_het_overlap/inversion_localpca.${CHROM}"

echo "[INFO] CHROM=${CHROM}"
echo "[INFO] Using Rscript: ${RSCRIPT_BIN}"

if [[ -s "${OUTPREFIX}.candidate_theta_summary.tsv.gz" && -s "${OUTPREFIX}.candidate_theta_detail.tsv.gz" ]]; then
  echo "[INFO] STEP10b/11 outputs already exist for ${CHROM}, skipping"
  exit 0
fi

"${RSCRIPT_BIN}" "${SCRIPTDIR}/core_code/STEP10b_parse_pestPG_to_sample_theta_windows_faster.R" \
  "${HETDIR}/02_heterozygosity/03_theta/multiscale" \
  "${HETDIR}/01_inputs_check/samples.ind" \
  "${THETA_PREFIX}" \
  "win50000.step10000" \
  "${CHROM}"

[[ -s "${THETA_PREFIX}.sample_tP_windows.tsv.gz" ]] || { echo "[ERROR] Missing parsed theta output for ${CHROM}" >&2; exit 1; }

"${RSCRIPT_BIN}" -e '
library(data.table)
dt <- fread("'"${THETA_PREFIX}"'.sample_tP_windows.tsv.gz")
pop_mean <- dt[, .(tP = mean(tP, na.rm=TRUE), nSamples = .N),
               by = .(chrom, WinStart, WinStop, WinCenter)]
fwrite(pop_mean, "'"${POP_TP}"'", sep="\t")
'

[[ -s "${POP_TP}" ]] || { echo "[ERROR] Missing population mean theta output for ${CHROM}" >&2; exit 1; }

"${RSCRIPT_BIN}" "${SCRIPTDIR}/core_code/STEP11_overlap_candidate_regions_with_theta_and_het.R" \
  "${INVDIR}/06_mds_candidates/inversion_localpca.candidate_regions.tsv.gz" \
  "${POP_TP}" \
  "${OUTPREFIX}" \
  chrom WinStart WinStop tP

echo "[DONE] ${CHROM}"
echo "  ${THETA_PREFIX}.sample_tP_windows.tsv.gz"
echo "  ${POP_TP}"
echo "  ${OUTPREFIX}.candidate_theta_summary.tsv.gz"
echo "  ${OUTPREFIX}.candidate_theta_detail.tsv.gz"
