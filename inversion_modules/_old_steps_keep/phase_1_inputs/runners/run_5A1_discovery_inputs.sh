#!/usr/bin/env bash
# =============================================================================
# run_5A1_discovery_inputs.sh — Runner for MODULE_5A1 (STEP01–STEP07)
#
# Upstream preparation: masking, ANGSD SAF/SFS, SNP calling, BEAGLE.
# Everything that prepares the data needed for local structure discovery.
#
# Usage:
#   bash run_5A1_discovery_inputs.sh                 # run all STEP01–07
#   bash run_5A1_discovery_inputs.sh --step 03       # run STEP03 only
#   bash run_5A1_discovery_inputs.sh --from 04 --to 07
# =============================================================================

set -euo pipefail
source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../00_inversion_config.sh"
inv_init_dirs

# Resolve 5A1 paths
A1_STEPS="${SCRIPT_DIR}/steps"
A1_LAUNCHERS="${SCRIPT_DIR}/launchers"

# ── Parse arguments ──────────────────────────────────────────────────────
STEP_ONLY=""
FROM_STEP=""
TO_STEP=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --step)  STEP_ONLY="$2"; shift 2 ;;
    --from)  FROM_STEP="$2"; shift 2 ;;
    --to)    TO_STEP="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

A1_STEPS_LIST=(01 02 03 04 05 06 07)

should_run() {
  local step="$1"
  if [[ -n "${STEP_ONLY}" ]]; then
    [[ "${step}" == "${STEP_ONLY}" ]]; return $?
  fi
  if [[ -n "${FROM_STEP}" || -n "${TO_STEP}" ]]; then
    local started=false
    for s in "${A1_STEPS_LIST[@]}"; do
      [[ "$s" == "${FROM_STEP:-01}" ]] && started=true
      if $started && [[ "$s" == "$step" ]]; then return 0; fi
      [[ "$s" == "${TO_STEP:-07}" ]] && $started && return 1
    done
    return 1
  fi
  return 0
}

echo "================================================================"
echo "  MODULE_5A1 — Discovery Inputs (STEP01–07)"
echo "  Started: $(date)"
echo "================================================================"

# ── STEP01+02 ─────────────────────────────────────────────────────────
MASK_PREFIX="${BASE}/fClaHyb_Gar_LG.mask_regions"
if should_run 01 || should_run 02; then
  if [[ ! -f "${MASK_PREFIX}.inversion_acgt_allcase.1-based.angsd.idx" ]]; then
    inv_log "STEP01+02: Build inversion callable mask"
    bash "${A1_STEPS}/STEP02_make_inversion_callable_angsd_mask.sh" \
        --fasta "${REF}" \
        --mask-script "${A1_STEPS}/STEP01_mask_regions_from_fasta.py" \
        --prefix "${MASK_PREFIX}"
  else
    inv_log "[SKIP] Mask already exists"
  fi
fi

# ── STEP03 ─────────────────────────────────────────────────────────────
if should_run 03; then
  inv_log "STEP03: Depth/mapQ stats"
  python3 "${A1_STEPS}/STEP03_mask_depth_mapq_stats.py" "${BAMLIST}" "${INVDIR}/03_sites"
fi

# ── STEP04–07 (SLURM array jobs) ──────────────────────────────────────
CHUNK_RF_LIST="${INVDIR}/chunk_rf.list"

if should_run 04; then
  inv_log "STEP04: Submit SAF chunks via sbatch:"
  N_CHUNKS=$(wc -l < "${CHUNK_RF_LIST}" 2>/dev/null || echo "?")
  echo "  sbatch --array=0-$((N_CHUNKS-1))%8 ${A1_LAUNCHERS}/LAUNCH_STEP04_run_angsd_saf_chunks_inversion.slurm ${CHUNK_RF_LIST}"
fi

if should_run 05; then
  inv_log "STEP05: Submit merge SAFs via sbatch:"
  echo "  sbatch ${A1_LAUNCHERS}/LAUNCH_STEP05_merge_chunk_saf_to_global_inversion.slurm"
fi

if should_run 06; then
  inv_log "STEP06: Submit folded SFS via sbatch:"
  echo "  sbatch ${A1_LAUNCHERS}/LAUNCH_STEP06_make_global_folded_sfs_and_mean_pest_inversion.slurm"
fi

if should_run 07; then
  inv_log "STEP07: Submit SNP calling via sbatch:"
  echo "  sbatch --array=0-\$((N_CHUNKS-1))%8 ${A1_LAUNCHERS}/LAUNCH_STEP07_angsd_call_snps_and_beagle_inversion.slurm ${CHUNK_RF_LIST}"
fi

echo ""
echo "  MODULE_5A1 complete — $(date)"
