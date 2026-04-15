#!/usr/bin/env bash
#SBATCH -J ngsadmix_local
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH -t 2-00:00:00
#SBATCH -o ngsadmix_local.%A_%a.out
#SBATCH -e ngsadmix_local.%A_%a.err
set -euo pipefail
source ~/.bashrc; mamba activate assembly
source "$(dirname "$0")/../config.sh"

LIST="${1:?provide beagle list}"
WORKDIR="${THIN_DIR}/05_ngsadmix_local_byRF"
mkdir -p "${WORKDIR}/logs"; cd "${WORKDIR}"
P=${SLURM_CPUS_PER_TASK:-32}

task=${SLURM_ARRAY_TASK_ID}
IN=$(sed -n "${task}p" "$LIST")
[[ "$IN" == /* ]] || IN="${THIN_DIR}/${IN}"
[[ -s "$IN" ]] || { echo "[ERROR] Missing: $IN"; exit 1; }

bn=$(basename "$IN"); tag="${bn%.beagle.gz}"
LG=$(echo "$tag" | sed -n 's/.*\(LG[0-9]\+\).*/\1/p'); [[ -n "$LG" ]] || LG="LGNA"
thin=$(echo "$tag" | sed -n 's/.*thin_\([0-9]\+\).*/\1/p'); [[ -n "$thin" ]] || thin="NA"

OUTROOT="runs_${LG}_thin${thin}"; mkdir -p "$OUTROOT"

for SEED in "${SEEDS[@]}"; do
  for K in $(seq $K_MIN $K_MAX); do
    OUTPREFIX="${OUTROOT}/${tag}_K$(printf "%02d" "$K")_seed${SEED}"
    NGSadmix -likes "$IN" -K "$K" -minMaf ${NGSADMIX_MINMAF} -P "$P" -seed "$SEED" -o "$OUTPREFIX"
  done
done
