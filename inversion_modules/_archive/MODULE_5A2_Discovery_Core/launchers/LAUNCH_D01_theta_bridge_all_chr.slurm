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

# =============================================================================
# LAUNCH_STEP10b_step11_all_chr.slurm  (v8.5)
#
# Per-chromosome theta bridge + candidate overlap.
# Runs across ALL chromosomes (not just candidate chromosomes).
#
# Usage:
#   sbatch --array=0-27%8 LAUNCH_STEP10b_step11_all_chr.slurm chrom.list
# =============================================================================

source ~/.bashrc
mamba activate assembly

CONFIG="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_codebase_v8.5/00_inversion_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }

set -a
source "${CONFIG}"
set +a

inv_init_dirs

CHROM_LIST="${1:?ERROR: provide chrom list as argv1}"

IDX="${SLURM_ARRAY_TASK_ID}"
LINE=$((IDX + 1))
CHROM="$(sed -n "${LINE}p" "${CHROM_LIST}")"

if [[ -z "${CHROM}" ]]; then
  echo "[ERROR] Could not read chromosome from ${CHROM_LIST} line ${LINE}" >&2
  exit 1
fi

A2_STEPS="${DISCOVERY2DIR}/steps"
A3_STEPS="${DISCOVERY3DIR}/steps"

THETA_PREFIX="${BRIDGE_DIR}/het_bridge_${CHROM}"
POP_TP="${BRIDGE_DIR}/population_mean_tP_windows.${CHROM}.tsv.gz"
OUTPREFIX="${BRIDGE_DIR}/inversion_localpca.${CHROM}"

echo "[INFO] CHROM=${CHROM}"
echo "[INFO] Using Rscript: ${RSCRIPT_BIN}"

if [[ -s "${OUTPREFIX}.candidate_theta_summary.tsv.gz" && -s "${OUTPREFIX}.candidate_theta_detail.tsv.gz" ]]; then
  echo "[INFO] STEP10b/11 outputs already exist for ${CHROM}, skipping"
  exit 0
fi

"${RSCRIPT_BIN}" "${A2_STEPS}/STEP_D01_theta_bridge.R" \
  "${HETDIR}/02_heterozygosity/03_theta/multiscale" \
  "${SAMPLES_IND}" \
  "${THETA_PREFIX}" \
  "${THETA_SCALE_TAG_DEFAULT}" \
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

"${RSCRIPT_BIN}" "${A3_STEPS}/STEP11_overlap_candidate_regions_with_theta_and_het.R" \
  "${CAND_FILE}" \
  "${POP_TP}" \
  "${OUTPREFIX}" \
  chrom WinStart WinStop tP

echo "[DONE] ${CHROM}"
echo "  ${THETA_PREFIX}.sample_tP_windows.tsv.gz"
echo "  ${POP_TP}"
echo "  ${OUTPREFIX}.candidate_theta_summary.tsv.gz"
echo "  ${OUTPREFIX}.candidate_theta_detail.tsv.gz"
