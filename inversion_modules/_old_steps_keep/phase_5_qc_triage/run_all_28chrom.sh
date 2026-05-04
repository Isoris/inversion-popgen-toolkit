#!/bin/bash
# =============================================================================
# run_all_28chrom.sh — end-to-end pipeline + registry population
# =============================================================================
# Runs the full Q01→Q10 pipeline across ALL 28 chromosomes (via SLURM array),
# waits for completion, stitches the per-chrom outputs into global tables,
# and prints a registry summary.
#
# This is the "push button, find all inversions in the cohort, database it"
# script. Intended to be run ONCE to populate the registry for a new genome,
# or RE-RUN incrementally after pipeline updates.
#
# Usage:
#   bash run_all_28chrom.sh [--no-array] [--resume] [--query-only]
#
# Flags:
#   --no-array     run sequentially instead of SLURM array (for login node
#                  debugging; NOT recommended for 28 chroms = ~10 hrs wall)
#   --resume       skip chromosomes that already have Q10 registry entries
#                  (useful for re-runs after fixing something)
#   --query-only   skip all compute, just print the registry summary
# =============================================================================
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/00_config.sh"

NO_ARRAY=0
RESUME=0
QUERY_ONLY=0
for arg in "$@"; do
  case "${arg}" in
    --no-array)   NO_ARRAY=1 ;;
    --resume)     RESUME=1 ;;
    --query-only) QUERY_ONLY=1 ;;
  esac
done

if [[ "${QUERY_ONLY}" == "1" ]]; then
  qc_log "=== REGISTRY QUERY ONLY ==="
  bash "${here}/scripts/registry_query.sh" summary
  exit 0
fi

# ---- Phase 1: Discover chromosomes ------------------------------------------
CHR_LIST="${BEAGLE_DIR}/chr.list"
[[ -f "${CHR_LIST}" ]] || qc_die "chr.list not found at ${CHR_LIST}"
n_chroms=$(wc -l < "${CHR_LIST}")
qc_log "Found ${n_chroms} chromosomes in ${CHR_LIST}"

# ---- Phase 2: Dispatch compute ---------------------------------------------
if [[ "${NO_ARRAY}" == "1" ]]; then
  qc_log "=== Sequential mode (no SLURM array) ==="
  while IFS= read -r chr; do
    [[ -z "${chr}" ]] && continue
    # Skip if already in registry (resume mode)
    if [[ "${RESUME}" == "1" ]]; then
      if grep -q "^inv_${chr}_" "${QC_LOGS}/q10_${chr}.log" 2>/dev/null; then
        qc_log "  SKIP ${chr}: already registered (resume mode)"
        continue
      fi
    fi
    qc_log "--- ${chr} ---"
    bash "${here}/run_chrom.sh" "${chr}"
  done < "${CHR_LIST}"
else
  qc_log "=== SLURM array mode ==="
  mkdir -p "${here}/logs"
  job_info=$(sbatch --parsable "${here}/slurm/array_28chrom.sh")
  qc_log "Submitted SLURM array: ${job_info}"
  qc_log "Monitor:   squeue -j ${job_info}"
  qc_log "Wait for all tasks to complete, then re-run:"
  qc_log "  bash run_all_28chrom.sh --query-only"
  exit 0
fi

# ---- Phase 3: Global stitching ---------------------------------------------
qc_log "=== Stitching genome-wide gap table ==="
bash "${here}/STEP_Q09_gap_characterization.sh" ALL

qc_log "=== Batch-registering any missed candidates ==="
bash "${here}/STEP_Q10_register.sh" ALL

# ---- Phase 4: Summary ------------------------------------------------------
qc_log ""
qc_log "=== DONE. REGISTRY SUMMARY: ==="
bash "${here}/scripts/registry_query.sh" summary

qc_log ""
qc_log "Inspect individual candidates:"
qc_log "  bash scripts/registry_query.sh list_candidates"
qc_log "  bash scripts/registry_query.sh describe <cid>"
qc_log "  bash scripts/registry_query.sh fst <cid>"
qc_log "  bash scripts/registry_query.sh karyotypes <cid>"
