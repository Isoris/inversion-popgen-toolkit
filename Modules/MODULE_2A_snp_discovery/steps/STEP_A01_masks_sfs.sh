#!/usr/bin/env bash
###############################################################################
# STEP_A01_masks_sfs.sh — Callable mask, RF chunks, SAF, merged SFS, pest
#
# Combines old: S00_mask_regions, S01_callable_mask, S02_chr_rf_chunks,
#               S03_angsd_saf_chunks, S04_merge_saf, S05_folded_sfs
#
# This is NOT a SLURM script — it prints sbatch commands for array jobs.
# For interactive/sequential use, set RUN_MODE=local.
#
# Called by: run_step1.sh masks_sfs
###############################################################################
set -euo pipefail
source "$(dirname "$0")/../config.sh"

MASK_SCRIPT="$(dirname "$0")/mask_regions_from_fasta.py"
OUTDIR="${GLOBAL_DIR}"
mkdir -p "${OUTDIR}"/{chunks,global_sfs,rf_files,logs}

ARGFILE="${OUTDIR}/01_masks_sfs.arg"
{
  echo -e "key\tvalue"
  echo -e "step\t01_masks_sfs"
  echo -e "datetime\t$(timestamp)"
  echo -e "host\t$(hostname)"
  echo -e "reference\t${REF}"
  echo -e "mask_prefix\t${MASK_PREFIX}"
  echo -e "callable_sites\t${CALLABLE_SITES}"
  echo -e "angsd_doSaf\t5"
  echo -e "sfs_bootstraps\t${SFS_BOOTSTRAPS}"
  for v in ANGSD_GL ANGSD_MINQ ANGSD_MINMAPQ ANGSD_BAQ ANGSD_C ANGSD_MINDEPTHIND ANGSD_MAXDEPTHIND ANGSD_MININD; do
    echo -e "${v}\t${!v}"
  done
} > "$ARGFILE"

# ---- 1) Mask regions from FASTA ----
if [[ ! -s "${MASK_PREFIX}.normalACGT.bed" ]]; then
  echo "[$(timestamp)] Running mask_regions_from_fasta.py..."
  python3 "${MASK_SCRIPT}" -i "${REF}" -p "${MASK_PREFIX}"
else
  echo "[$(timestamp)] [SKIP] Mask BEDs already exist"
fi

# ---- 2) Build ANGSD callable sites (0-based BED → 1-based ANGSD) ----
if [[ ! -s "${CALLABLE_SITES}" ]]; then
  echo "[$(timestamp)] Converting to ANGSD 1-based sites..."
  awk 'BEGIN{OFS="\t"} {print $1, $2+1, $3}' "${MASK_PREFIX}.normalACGT.bed" > "${CALLABLE_SITES}"
  angsd sites index "${CALLABLE_SITES}"
else
  echo "[$(timestamp)] [SKIP] ANGSD sites already indexed"
fi

# ---- 3) Create per-chromosome RF files ----
CHUNK_LIST="${OUTDIR}/chunk_rf.list"
RF_DIR="${OUTDIR}/rf_files"
if [[ ! -s "${CHUNK_LIST}" ]]; then
  echo "[$(timestamp)] Creating RF chunk files..."
  mkdir -p "${RF_DIR}"
  awk '/^>/ { name = substr($1, 2); if (name != "") print name }' "${REF}" \
  | while read -r chr; do
    printf "%s\n" "$chr" > "${RF_DIR}/${chr}.rf.txt"
  done
  find "${RF_DIR}" -maxdepth 1 -name "*.rf.txt" | sort > "${CHUNK_LIST}"
  echo "[$(timestamp)] Created $(wc -l < "${CHUNK_LIST}") chunks"
else
  echo "[$(timestamp)] [SKIP] Chunk list exists"
fi

N_CHUNKS=$(wc -l < "${CHUNK_LIST}")

# ---- 4) Per-chunk SAF (SLURM array) ----
echo "[$(timestamp)] === SAF per chunk ==="
echo "[INFO] Submit as SLURM array:"
echo "  sbatch --array=0-$((N_CHUNKS-1))%8 slurm/SLURM_A01a_saf_chunk.sh ${CHUNK_LIST}"
echo ""
echo "[INFO] Or run sequentially (slow):"
echo "  for i in \$(seq 0 $((N_CHUNKS-1))); do SLURM_ARRAY_TASK_ID=\$i bash slurm/SLURM_A01a_saf_chunk.sh ${CHUNK_LIST}; done"

# ---- 5) Merge chunk SAFs ----
GLOBAL_SAF="${OUTDIR}/global_sfs/catfish.global.saf.idx"
if [[ -s "${GLOBAL_SAF}" ]]; then
  echo "[$(timestamp)] [SKIP] Global SAF already exists"
else
  echo "[$(timestamp)] === Merge SAFs ==="
  echo "[INFO] After SAF chunks complete, run:"
  echo "  sbatch slurm/SLURM_A01b_merge_sfs.sh"
fi

# ---- 6) Folded SFS + mean pest ----
PEST="${OUTDIR}/global_sfs/catfish.global.folded.mean.pest"
if [[ -s "${PEST}" ]]; then
  echo "[$(timestamp)] [SKIP] Pest prior already exists"
else
  echo "[$(timestamp)] === Folded SFS ==="
  echo "[INFO] After SAF merge, run:"
  echo "  sbatch slurm/SLURM_A01b_merge_sfs.sh"
fi

echo "[$(timestamp)] [DONE] 01_masks_sfs"
