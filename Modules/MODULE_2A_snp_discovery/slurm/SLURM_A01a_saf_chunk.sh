#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=80
#SBATCH --mem=64GB
#SBATCH -t 1-00:00:00
#SBATCH -J angsd_saf
#SBATCH -o angsd_saf.%A_%a.out
#SBATCH -e angsd_saf.%A_%a.err
set -euo pipefail
source ~/.bashrc; mamba activate assembly
source "$(dirname "$0")/../config.sh"

CHUNKLIST="${1:?ERROR: provide chunk_rf.list}"
RF_FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$CHUNKLIST")
[[ -n "$RF_FILE" && -s "$RF_FILE" ]] || { echo "[ERROR] Bad RF" >&2; exit 1; }
TAG=$(basename "$RF_FILE" .rf.txt)

CHUNKDIR="${GLOBAL_DIR}/chunks/${TAG}"
mkdir -p "${CHUNKDIR}"/{01_saf,02_sfs,logs}
SAF_PREFIX="${CHUNKDIR}/01_saf/catfish.${TAG}"
SFS_OUT="${CHUNKDIR}/02_sfs/catfish.${TAG}.sfs"

P="${SLURM_CPUS_PER_TASK}"

angsd -b "${BAMLIST}" -ref "${REF}" \
  -GL ${ANGSD_GL} -doSaf 5 -doMajorMinor 1 -doMaf 1 \
  -doCounts 1 -doDepth 1 \
  -minQ ${ANGSD_MINQ} -minMapQ ${ANGSD_MINMAPQ} -baq ${ANGSD_BAQ} -C ${ANGSD_C} \
  -setMinDepthInd ${ANGSD_MINDEPTHIND} -setMaxDepthInd ${ANGSD_MAXDEPTHIND} \
  -minInd ${ANGSD_MININD} \
  -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 \
  -rf "${RF_FILE}" -sites "${CALLABLE_SITES}" \
  -P "${P}" -out "${SAF_PREFIX}" \
  2> "${CHUNKDIR}/logs/angsd_saf.${TAG}.log"

realSFS "${SAF_PREFIX}.saf.idx" -bootstrap ${SFS_BOOTSTRAPS} -P "${P}" \
  > "${SFS_OUT}" 2> "${CHUNKDIR}/logs/realSFS.${TAG}.log"

echo "[DONE] ${TAG}"
