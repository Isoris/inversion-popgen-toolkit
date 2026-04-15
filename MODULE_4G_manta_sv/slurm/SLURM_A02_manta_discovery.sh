#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 127
#SBATCH --mem=237G
#SBATCH -t 2-00:00:00
#SBATCH -J manta_disc
#SBATCH -o logs/02_discovery.%j.out
#SBATCH -e logs/02_discovery.%j.err
# =============================================================================
# 02_manta_discovery.sh — Per-sample Manta SV calling
# =============================================================================
# Manta is a two-step caller:
#   1. configManta.py  — configure the run (inputs, reference, options)
#   2. runWorkflow.py  — execute (produces diploidSV.vcf.gz per sample)
#
# Unlike DELLY, Manta calls ALL SV types in one pass (DEL, DUP, INS, INV/BND).
# There is no separate merge/regenotype workflow — each sample is independent.
# =============================================================================
set -euo pipefail
source ~/.bashrc
mamba activate manta_py2
export PATH="/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin:$PATH"

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${SCRIPT_DIR}/../00_module4g_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
set -a
source "${CONFIG}"
set +a

mv_init_dirs
mv_log "=== STEP 2: Per-sample Manta SV discovery ==="

# Validate
mv_check_file "$REF" "Reference"
mv_check_file "$REF_FAI" "Reference index"
mv_check_file "$CALL_REGIONS_BED_GZ" "Call regions BED"
mv_check_file "${CALL_REGIONS_BED_GZ}.tbi" "Call regions tabix index"
mv_check_file "$MANTA_CUSTOM_INI" "Custom Manta ini"
mv_check_file "$MANTA_CONFIG" "configManta.py"
mv_check_file "$SAMPLES_ALL" "Sample list"

N_SAMPLES=$(wc -l < "$SAMPLES_ALL")
mv_log "Samples: ${N_SAMPLES}"
mv_log "Manta threads/call: ${MANTA_THREADS_PER_CALL}, parallel: ${MANTA_PARALLEL}"
mv_log "Custom ini: ${MANTA_CUSTOM_INI}"

# ─────────────────────────────────────────────────────────────────────────────
# Per-sample Manta run
# ─────────────────────────────────────────────────────────────────────────────
run_manta_sample() {
  local sid="$1"
  local bam="${DIR_MARKDUP}/${sid}.markdup.bam"
  local rundir="${DIR_PERSAMPLE}/${sid}"
  local final_vcf="${rundir}/results/variants/diploidSV.vcf.gz"
  local logf="${DIR_LOGS}/manta_${sid}.log"

  # Skip if final output exists
  if [[ -f "${final_vcf}" ]]; then
    echo "[$(date '+%T')] ${sid}: diploidSV.vcf.gz exists, skipping"
    return 0
  fi

  if [[ ! -f "${bam}" ]]; then
    echo "[$(date '+%T')] ${sid}: ERROR — markdup BAM missing: ${bam}" >&2
    return 1
  fi

  echo "[$(date '+%T')] ${sid}: configuring Manta..."

  # Clean up any partial run
  rm -rf "${rundir}"

  # Step 1: Configure
  python2 "${MANTA_CONFIG}" \
    --bam "${bam}" \
    --referenceFasta "${REF}" \
    --callRegions "${CALL_REGIONS_BED_GZ}" \
    --config "${MANTA_CUSTOM_INI}" \
    --runDir "${rundir}" \
    2>>"${logf}"

  # Step 2: Execute (local mode, limited threads per sample)
  echo "[$(date '+%T')] ${sid}: running Manta workflow..."
  python2 "${rundir}/runWorkflow.py" \
    -m local \
    -j "${MANTA_THREADS_PER_CALL}" \
    2>>"${logf}"

  # Verify output
  if [[ ! -f "${final_vcf}" ]]; then
    echo "[$(date '+%T')] ${sid}: ERROR — diploidSV.vcf.gz not produced" >&2
    return 1
  fi

  local n_sv
  n_sv=$(bcftools view -H "${final_vcf}" 2>/dev/null | wc -l)
  echo "[$(date '+%T')] ${sid}: ${n_sv} SVs in diploidSV.vcf.gz"
}

export -f run_manta_sample
export DIR_MARKDUP REF REF_FAI CALL_REGIONS_BED_GZ MANTA_CONFIG MANTA_CUSTOM_INI
export DIR_PERSAMPLE DIR_LOGS MANTA_THREADS_PER_CALL

parallel -j "${MANTA_PARALLEL}" --joblog "${DIR_LOGS}/discovery_parallel.log" \
  run_manta_sample {} \
  :::: "${SAMPLES_ALL}"

# ── Verify all samples ──────────────────────────────────────────────────────
mv_log "Verifying outputs..."
FAIL=0
while IFS= read -r sid; do
  vcf="${DIR_PERSAMPLE}/${sid}/results/variants/diploidSV.vcf.gz"
  [[ -f "${vcf}" ]] || { mv_err "Missing: ${sid}"; ((FAIL++)) || true; }
done < "$SAMPLES_ALL"
[[ ${FAIL} -eq 0 ]] || mv_die "${FAIL} samples failed. Check ${DIR_LOGS}/manta_*.log"

# ── Summary counts ──────────────────────────────────────────────────────────
{
  echo -e "sample\tn_SV"
  while IFS= read -r sid; do
    vcf="${DIR_PERSAMPLE}/${sid}/results/variants/diploidSV.vcf.gz"
    n=$(bcftools view -H "${vcf}" 2>/dev/null | wc -l)
    echo -e "${sid}\t${n}"
  done < "${SAMPLES_ALL}"
} | tee "${DIR_LOGS}/discovery_SV_counts.tsv"

mv_log "=== STEP 2 COMPLETE ==="
