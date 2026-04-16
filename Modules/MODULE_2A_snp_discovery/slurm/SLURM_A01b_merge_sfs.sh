#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 124
#SBATCH --mem=237GB
#SBATCH -t 2-00:00:00
#SBATCH -J merge_sfs
#SBATCH -o merge_sfs.%j.out
#SBATCH -e merge_sfs.%j.err
set -euo pipefail
source ~/.bashrc; mamba activate assembly
source "$(dirname "$0")/../config.sh"

# ---- Merge chunk SAFs ----
GLOBALDIR="${GLOBAL_DIR}/global_sfs"
mkdir -p "$GLOBALDIR"

SAF_LIST="${GLOBALDIR}/saf.idx.list"
GLOBAL_PREFIX="${GLOBALDIR}/catfish.global"
P="${SLURM_CPUS_PER_TASK:-124}"

find "${GLOBAL_DIR}/chunks" -type f -name "*.saf.idx" | sort > "$SAF_LIST"
N=$(wc -l < "$SAF_LIST")
(( N > 0 )) || { echo "[ERROR] No .saf.idx files" >&2; exit 1; }
echo "[INFO] Merging $N SAF files..."

rm -f "${GLOBAL_PREFIX}.saf.idx" "${GLOBAL_PREFIX}.saf.gz" "${GLOBAL_PREFIX}.saf.pos.gz"
realSFS cat -P "$P" -outnames "$GLOBAL_PREFIX" $(cat "$SAF_LIST")

# ---- Folded SFS + mean pest ----
OUT_SFS="${GLOBALDIR}/catfish.global.folded.sfs"
OUT_PEST="${GLOBALDIR}/catfish.global.folded.mean.pest"

echo "[INFO] Computing folded SFS (${SFS_BOOTSTRAPS} bootstraps)..."
realSFS "${GLOBAL_PREFIX}.saf.idx" -P "$P" -bootstrap ${SFS_BOOTSTRAPS} -fold 1 \
  > "$OUT_SFS" 2> "${GLOBALDIR}/realSFS.global.folded.log"

# Average bootstrap rows into mean pest prior
awk '{ for(i=1;i<=NF;i++) s[i]+=$i; n++ }
END { for(i=1;i<=NF;i++) printf "%.10f%s", s[i]/n, (i<NF?" ":"\n") }' "$OUT_SFS" > "$OUT_PEST"

echo "[DONE] SAF: ${GLOBAL_PREFIX}.saf.idx"
echo "[DONE] SFS: $OUT_SFS"
echo "[DONE] PEST: $OUT_PEST"
