#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=32G
#SBATCH -t 0-02:00:00
#SBATCH -J delly_dup_prep
#SBATCH -o logs/01_prep.%j.out
#SBATCH -e logs/01_prep.%j.err

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
dv_log "=== STEP 1: Reuse existing DEL resources for DUP pipeline ==="

dv_check_file "$REF" "Reference"
dv_check_file "$REF_FAI" "Reference index"
dv_check_file "$EXCL_BED" "Exclude BED"
dv_check_file "$SAMPLES_ALL" "All-sample list"
dv_check_file "$SAMPLES_UNRELATED" "Unrelated sample list"
[[ -x "${DELLY_BIN}" ]] || dv_die "DELLY binary not executable: ${DELLY_BIN}"

N_ALL=$(wc -l < "$SAMPLES_ALL")
N_UNREL=$(wc -l < "$SAMPLES_UNRELATED")

dv_log "All samples: ${N_ALL}"
dv_log "Unrelated samples: ${N_UNREL}"
dv_log "Checking markdup BAMs in ${DIR_MARKDUP}"

FAIL=0
while IFS= read -r sid; do
  bam="${DIR_MARKDUP}/${sid}.markdup.bam"
  bai="${bam}.bai"
  if [[ ! -f "${bam}" ]]; then
    dv_err "Missing BAM: ${bam}"
    ((FAIL++)) || true
    continue
  fi
  if [[ ! -f "${bai}" ]]; then
    dv_err "Missing BAI: ${bai}"
    ((FAIL++)) || true
    continue
  fi
done < "$SAMPLES_ALL"

[[ ${FAIL} -eq 0 ]] || dv_die "${FAIL} markdup BAMs or indexes missing."

dv_log "All reused BAMs/indexes found."
dv_log "No duplicate marking was run."
dv_log "=== STEP 1 COMPLETE ==="
