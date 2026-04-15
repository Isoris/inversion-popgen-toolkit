#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=80
#SBATCH --mem=237GB
#SBATCH -t 2-00:00:00
#SBATCH -J ngsadmix_global
#SBATCH -o ngsadmix_global.%j.out
#SBATCH -e ngsadmix_global.%j.err
set -euo pipefail
source ~/.bashrc; mamba activate assembly
source "$(dirname "$0")/../config.sh"

# Global NGSadmix: thin-500, whole-genome BEAGLE, all K × seeds
BEAGLE="${THIN_DIR}/03_merged_beagle/catfish.wholegenome.byRF.thin_${RELATE_THIN}.beagle.gz"
[[ -s "$BEAGLE" ]] || { echo "[ERROR] Missing global BEAGLE: $BEAGLE"; exit 1; }

RUN_DIR="${THIN_DIR}/05_ngsadmix_global/runs_thin${RELATE_THIN}"
mkdir -p "$RUN_DIR"
P=${SLURM_CPUS_PER_TASK}

for SEED in "${SEEDS[@]}"; do
  for K in $(seq $K_MIN $K_MAX); do
    PREFIX="${RUN_DIR}/thin${RELATE_THIN}_K$(printf "%02d" "$K")_seed${SEED}"
    [[ -s "${PREFIX}.qopt" ]] && { echo "[SKIP] ${PREFIX}"; continue; }
    echo "[RUN] K=$K seed=$SEED"
    NGSadmix -likes "$BEAGLE" -K "$K" -minMaf ${NGSADMIX_MINMAF} -P "$P" -seed "$SEED" -o "$PREFIX"
  done
done

echo "[DONE] Global NGSadmix thin-${RELATE_THIN}"
