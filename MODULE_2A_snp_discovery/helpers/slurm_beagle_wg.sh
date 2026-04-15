#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=80
#SBATCH --mem=237GB
#SBATCH -t 2-00:00:00
#SBATCH -J beagle_wg
#SBATCH -o beagle_wg.%j.out
#SBATCH -e beagle_wg.%j.err
set -euo pipefail
source ~/.bashrc; mamba activate assembly
source "$(dirname "$0")/../config.sh"

W="${1:?provide thinning bp}"
SITES="${THIN_DIR}/03_sites/sites.thin_${W}.majmin.tsv"
BEAGLEDIR="${THIN_DIR}/04_beagle_wholegenome_majmin"
mkdir -p "${BEAGLEDIR}" "${THIN_DIR}/logs"
P="${SLURM_CPUS_PER_TASK}"

[[ -s "${SITES}" ]] || { echo "[ERROR] Missing sites"; exit 1; }
PREFIX="${BEAGLEDIR}/catfish.wholegenome.thin_${W}"

angsd -b "${BAMLIST}" -ref "${REF}" \
  -GL ${ANGSD_GL} -doMajorMinor 3 -doGlf 2 \
  -doCounts 1 -doDepth 1 \
  -minQ ${ANGSD_MINQ} -minMapQ ${ANGSD_MINMAPQ} -baq ${ANGSD_BAQ} -C ${ANGSD_C} \
  -setMinDepthInd ${ANGSD_MINDEPTHIND} -setMaxDepthInd ${ANGSD_MAXDEPTHIND} \
  -minInd ${ANGSD_MININD} \
  -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 \
  -sites "${SITES}" -P "${P}" -out "${PREFIX}" \
  > "${THIN_DIR}/logs/beagle_wg.thin_${W}.log" 2>&1
