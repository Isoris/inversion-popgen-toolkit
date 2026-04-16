#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${SCRIPT_DIR}/../00_module4c_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }

set -a
source "${CONFIG}"
set +a

LOG_DIR="${OUTDIR}/logs"
mkdir -p "${LOG_DIR}"

for f in "${DIR_FINAL}/catalog_226.DUP.vcf.gz" \
         "${DIR_FINAL}/catalog_226.DUP.GT_matrix.tsv" \
         "${DIR_FINAL}/catalog_226.DUP.bed"; do
  [[ -f "$f" ]] || { echo "Missing required upstream output: $f" >&2; exit 1; }
done

JID1=$(sbatch --parsable "${SCRIPT_DIR}/04_annotation_layers.sh")
echo "[1/3] DUP annotation:   ${JID1}"

JID2=$(sbatch --parsable --dependency=afterok:${JID1} "${SCRIPT_DIR}/../slurm/SLURM_B01_downstream_analysis.sh")
echo "[2/3] DUP downstream:   ${JID2}"

JID3=$(sbatch --parsable --dependency=afterok:${JID2} "${SCRIPT_DIR}/../slurm/SLURM_B02_extended_plots.sh")
echo "[3/3] DUP plots:        ${JID3}"

echo "Chain: ${JID1} -> ${JID2} -> ${JID3}"
