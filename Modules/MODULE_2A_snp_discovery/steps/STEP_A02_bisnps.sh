#!/usr/bin/env bash
###############################################################################
# STEP_A02_bisnps.sh — biSNP discovery per chunk using pest prior
# Combines old: S06_run_bisnp_discovery_chunks
# Called by: run_step1.sh call_bisnps
###############################################################################
set -euo pipefail
source "$(dirname "$0")/../config.sh"

CHUNK_LIST="${GLOBAL_DIR}/chunk_rf.list"
PEST="${GLOBAL_DIR}/global_sfs/catfish.global.folded.mean.pest"

[[ -s "$CHUNK_LIST" ]] || { echo "[ERROR] Missing $CHUNK_LIST — run masks_sfs first" >&2; exit 1; }
[[ -s "$PEST" ]] || { echo "[ERROR] Missing pest prior — run masks_sfs first" >&2; exit 1; }

N=$(wc -l < "$CHUNK_LIST")
echo "[$(timestamp)] biSNP discovery: $N chunks"
echo "[INFO] Submit: sbatch --array=0-$((N-1))%8 slurm/SLURM_A02_bisnp_chunk.sh ${CHUNK_LIST}"
