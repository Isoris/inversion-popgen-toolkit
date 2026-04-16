#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=80
#SBATCH --mem=237GB
#SBATCH -t 2-00:00:00
#SBATCH -J angsd_bisnp
#SBATCH -o angsd_bisnp.%A_%a.out
#SBATCH -e angsd_bisnp.%A_%a.err
set -euo pipefail
source ~/.bashrc; mamba activate assembly
source "$(dirname "$0")/../config.sh"

RFLIST="${1:?ERROR: provide chunk_rf.list}"
PEST="${GLOBAL_DIR}/global_sfs/catfish.global.folded.mean.pest"
P="${SLURM_CPUS_PER_TASK}"

RF_FILE="$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$RFLIST")"
[[ -n "$RF_FILE" ]] || exit 1
TAG="$(basename "$RF_FILE" .rf.txt)"

mkdir -p "${GLOBAL_DIR}"/{02_snps,logs}
SNP_PREFIX="${GLOBAL_DIR}/02_snps/catfish.${TAG}"

angsd -b "$BAMLIST" -ref "$REF" \
  -GL ${ANGSD_GL} -doMajorMinor 1 -doMaf 1 -pest "$PEST" \
  -doCounts 1 -doDepth 1 \
  -minQ ${ANGSD_MINQ} -minMapQ ${ANGSD_MINMAPQ} -baq ${ANGSD_BAQ} -C ${ANGSD_C} \
  -setMinDepthInd ${ANGSD_MINDEPTHIND} -setMaxDepthInd ${ANGSD_MAXDEPTHIND} \
  -minInd ${ANGSD_MININD} \
  -SNP_pval ${SNP_PVAL} -minMaf ${MIN_MAF} \
  -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 \
  -rf "$RF_FILE" -sites "$CALLABLE_SITES" \
  -P "$P" -out "$SNP_PREFIX" \
  > "${GLOBAL_DIR}/logs/bisnp.${TAG}.log" 2>&1

echo "[DONE] ${TAG}"
