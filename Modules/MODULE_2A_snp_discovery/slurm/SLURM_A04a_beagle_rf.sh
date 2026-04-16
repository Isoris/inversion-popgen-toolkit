#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=80
#SBATCH --mem=237GB
#SBATCH -t 2-00:00:00
#SBATCH -J beagle_rf
#SBATCH -o beagle_rf.%A_%a.out
#SBATCH -e beagle_rf.%A_%a.err
set -euo pipefail
source ~/.bashrc; mamba activate assembly
source "$(dirname "$0")/../config.sh"

RFLIST="${1:?provide chunk_rf.list}"; W="${2:?provide thinning bp}"
SITES="${THIN_DIR}/03_sites/sites.thin_${W}.majmin.tsv"
BEAGLEDIR="${THIN_DIR}/04_beagle_byRF_majmin/thin_${W}"
mkdir -p "${BEAGLEDIR}" "${THIN_DIR}/logs"
P="${SLURM_CPUS_PER_TASK}"

[[ -s "${SITES}" ]] || { echo "[ERROR] Missing sites"; exit 1; }
[[ -s "${SITES}.idx" ]] || angsd sites index "${SITES}"

RF_FILE="$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$RFLIST")"
TAG="$(basename "$RF_FILE" .rf.txt)"
PREFIX="${BEAGLEDIR}/catfish.${TAG}.thin_${W}"

angsd -b "${BAMLIST}" -ref "${REF}" \
  -GL ${ANGSD_GL} -doMajorMinor 3 -doGlf 2 \
  -doCounts 1 -doDepth 1 \
  -minQ ${ANGSD_MINQ} -minMapQ ${ANGSD_MINMAPQ} -baq ${ANGSD_BAQ} -C ${ANGSD_C} \
  -setMinDepthInd ${ANGSD_MINDEPTHIND} -setMaxDepthInd ${ANGSD_MAXDEPTHIND} \
  -minInd ${ANGSD_MININD} \
  -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 \
  -rf "${RF_FILE}" -sites "${SITES}" \
  -P "${P}" -out "${PREFIX}" \
  > "${THIN_DIR}/logs/beagle_rf.${TAG}.thin_${W}.log" 2>&1
