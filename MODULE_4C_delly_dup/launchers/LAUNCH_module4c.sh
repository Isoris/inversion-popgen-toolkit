#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${SCRIPT_DIR}/../00_module4c_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }

set -a
source "${CONFIG}"
set +a

mkdir -p "${DIR_LOGS}"

echo "=== DELLY ${SVTYPE} Pipeline Launcher ==="
echo "Script directory: ${SCRIPT_DIR}"
echo "Reused markdup BAM dir: ${DIR_MARKDUP}"
echo "Reused exclude BED:     ${EXCL_BED}"
echo

JID1=$(sbatch --parsable "${SCRIPT_DIR}/../slurm/SLURM_A01_prep_inputs.sh")
echo "[1/4] Prep/check reuse:   ${JID1}"

JID2=$(sbatch --parsable --dependency=afterok:${JID1} "${SCRIPT_DIR}/../slurm/SLURM_A02_delly_discovery.sh")
echo "[2/4] Discovery:          ${JID2}"

JID3=$(sbatch --parsable --dependency=afterok:${JID2} "${SCRIPT_DIR}/../slurm/SLURM_A03_merge_genotype.sh")
echo "[3/4] Merge+genotype:     ${JID3}"

JID4=$(sbatch --parsable --dependency=afterok:${JID3} "${SCRIPT_DIR}/../slurm/SLURM_A05_summary_report.sh")
echo "[4/4] Summary:            ${JID4}"

echo
echo "Chain: ${JID1} -> ${JID2} -> ${JID3} -> ${JID4}"
