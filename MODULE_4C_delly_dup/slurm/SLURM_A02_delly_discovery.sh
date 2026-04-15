#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 127
#SBATCH --mem=237G
#SBATCH -t 2-00:00:00
#SBATCH -J delly_dup_disc
#SBATCH -o logs/02_discovery.%j.out
#SBATCH -e logs/02_discovery.%j.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${SCRIPT_DIR}/../00_module4c_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }

set -a
source "${CONFIG}"
set +a

dv_init_dirs
dv_log "=== STEP 2: Per-sample DELLY ${SVTYPE} discovery ==="

dv_check_file "$REF" "Reference"
dv_check_file "$EXCL_BED" "Exclusion BED"
dv_check_file "$SAMPLES_ALL" "Sample list"
[[ -x "${DELLY_BIN}" ]] || dv_die "DELLY not executable: ${DELLY_BIN}"

run_discovery() {
  local sid="$1"
  local bam="${DIR_MARKDUP}/${sid}.markdup.bam"
  local out_bcf="${DIR_DISC}/${sid}.disc.${SVTYPE}.bcf"
  local logf="${DIR_LOGS}/disc_${sid}.log"

  if [[ -f "${out_bcf}" && -f "${out_bcf}.csi" ]]; then
    echo "[$(date '+%T')] ${sid}: already done, skipping"
    return 0
  fi

  [[ -f "${bam}" ]] || { echo "[$(date '+%T')] ${sid}: missing BAM ${bam}" >&2; return 1; }

  "${DELLY_BIN}" call \
    -t "${SVTYPE}" \
    -g "${REF}" \
    -x "${EXCL_BED}" \
    -h "${DELLY_THREADS_PER_CALL}" \
    -o "${out_bcf}" \
    "${bam}" \
    2>>"${logf}"

  bcftools index -f "${out_bcf}" 2>>"${logf}"

  local n
  n=$(bcftools view -H "${out_bcf}" 2>/dev/null | wc -l)
  echo "[$(date '+%T')] ${sid}: ${n} ${SVTYPE}"
}

export -f run_discovery
export DIR_MARKDUP REF EXCL_BED DELLY_BIN DELLY_THREADS_PER_CALL DIR_DISC DIR_LOGS SVTYPE

parallel -j "${DELLY_PARALLEL}" --joblog "${DIR_LOGS}/discovery_parallel.log" \
  run_discovery {} \
  :::: "${SAMPLES_ALL}"

FAIL=0
while IFS= read -r sid; do
  [[ -f "${DIR_DISC}/${sid}.disc.${SVTYPE}.bcf" ]] || { dv_err "Missing discovery output for ${sid}"; ((FAIL++)) || true; }
done < "$SAMPLES_ALL"

[[ ${FAIL} -eq 0 ]] || dv_die "${FAIL} samples failed discovery."

{
  echo -e "sample\tn_${SVTYPE}"
  for bcf in "${DIR_DISC}"/*.disc."${SVTYPE}".bcf; do
    sid=$(basename "$bcf" .disc."${SVTYPE}".bcf)
    n=$(bcftools view -H "$bcf" 2>/dev/null | wc -l)
    echo -e "${sid}\t${n}"
  done
} > "${DIR_LOGS}/discovery_${SVTYPE}_counts.tsv"

dv_log "=== STEP 2 COMPLETE ==="
