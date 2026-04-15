#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=32G
#SBATCH -t 0-02:00:00
#SBATCH -J delly_dup_sum
#SBATCH -o logs/05_summary.%j.out
#SBATCH -e logs/05_summary.%j.err

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
dv_log "=== FINAL SUMMARY (${SVTYPE}) ==="

REPORT="${OUTDIR}/summary_${SVTYPE}.txt"

FINAL_226_VCF="${DIR_FINAL}/catalog_226.${SVTYPE}.vcf.gz"
FINAL_81_PASS=$(ls "${DIR_FINAL}"/catalog_81."${SVTYPE}".*.PASS.vcf.gz 2>/dev/null | head -1)

N_226=0
[[ -f "${FINAL_226_VCF}" ]] && N_226=$(bcftools view -H "${FINAL_226_VCF}" | wc -l)

N_81=0
[[ -n "${FINAL_81_PASS}" && -f "${FINAL_81_PASS}" ]] && N_81=$(bcftools view -H "${FINAL_81_PASS}" | wc -l)

{
  echo "============================================================"
  echo "DELLY ${SVTYPE} SUMMARY"
  echo "Generated: $(date '+%F %T')"
  echo "============================================================"
  echo "226 strict catalog: ${N_226}"
  echo "81 strict catalog:  ${N_81}"
  echo
  echo "Files:"
  echo "  ${FINAL_226_VCF}"
  echo "  ${FINAL_81_PASS}"
  echo "============================================================"
} | tee "${REPORT}"

dv_log "=== SUMMARY COMPLETE ==="
