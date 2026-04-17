#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=64G
#SBATCH -t 0-06:00:00
#SBATCH -J pop_regeno
#SBATCH -o logs/A03b_pop_regenotype.%j.out
#SBATCH -e logs/A03b_pop_regenotype.%j.err
# =============================================================================
# SLURM_A03b_population_regenotype.sh
# =============================================================================
# Population-prior regenotyping for DELLY2 INV and BND VCFs.
#
# PIPELINE POSITION:
#   After:  SLURM_A03 (merge + regenotype → cohort_226 raw VCF)
#   Before: Strict filtering (which would otherwise use the dropout-biased GTs)
#
# WHAT IT DOES:
#   Reads raw FORMAT fields (DR/DV/RR/RV for DELLY, PR/SR for Manta) from
#   the merged 226-sample VCF. Estimates cohort allele frequency from
#   high-evidence samples, then applies a Bayesian prior to rescue
#   low-evidence samples that were misgenotyped as 0/0.
#
# USAGE:
#   cd MODULE_4D_delly_inv/slurm/   # or MODULE_4E_delly_bnd/slurm/
#   sbatch SLURM_A03b_population_regenotype.sh
#
# To run for Manta INV (MODULE_4G):
#   Set CALLER=manta and point RAW_VCF at the Manta merged VCF.
# =============================================================================
set -euo pipefail
source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"

# ── Auto-detect module ────────────────────────────────────────────────────
# Try MODULE_4D (INV) first, then MODULE_4E (BND)
if [[ -f "${SCRIPT_DIR}/../00_module4d_config.sh" ]]; then
  CONFIG="${SCRIPT_DIR}/../00_module4d_config.sh"
  SV_LABEL="INV"
  CALLER="delly"
elif [[ -f "${SCRIPT_DIR}/../00_module4e_config.sh" ]]; then
  CONFIG="${SCRIPT_DIR}/../00_module4e_config.sh"
  SV_LABEL="BND"
  CALLER="delly"
elif [[ -f "${SCRIPT_DIR}/../00_module4g_config.sh" ]]; then
  CONFIG="${SCRIPT_DIR}/../00_module4g_config.sh"
  SV_LABEL="INV"
  CALLER="manta"
else
  echo "ERROR: Cannot find module config in ${SCRIPT_DIR}/.." >&2
  exit 1
fi

set -a
source "${CONFIG}"
set +a

# ── Locate the regenotyping script ────────────────────────────────────────
REGENO_PY="${SCRIPT_DIR}/population_regenotype.py"
if [[ ! -f "${REGENO_PY}" ]]; then
  # Try parent directory
  REGENO_PY="${SCRIPT_DIR}/../utils/population_regenotype.py"
fi
[[ -f "${REGENO_PY}" ]] || { echo "ERROR: population_regenotype.py not found" >&2; exit 1; }

echo "========================================"
echo "POPULATION-PRIOR REGENOTYPING"
echo "========================================"
echo "Module:   ${SV_LABEL} (${CALLER})"
echo "Config:   ${CONFIG}"
echo "Script:   ${REGENO_PY}"
echo ""

# ── Locate raw merged VCF ─────────────────────────────────────────────────
# This is the 226-sample merged VCF BEFORE strict filtering
if [[ "${CALLER}" == "delly" ]]; then
  RAW_VCF="${DIR_FINAL}/catalog_226.${SV_LABEL}.raw.vcf.gz"
  # Fallback: the merged BCF
  if [[ ! -f "${RAW_VCF}" ]]; then
    RAW_VCF="${DIR_MERGED}/cohort_226.${SV_LABEL}.merged.bcf"
  fi
elif [[ "${CALLER}" == "manta" ]]; then
  RAW_VCF="${DIR_SPLIT}/catalog_226.${SV_LABEL}.vcf.gz"
  if [[ ! -f "${RAW_VCF}" ]]; then
    RAW_VCF="${DIR_MERGED}/cohort_226.ALL.inv_converted.vcf.gz"
  fi
fi

[[ -f "${RAW_VCF}" ]] || { echo "ERROR: Raw VCF not found: ${RAW_VCF}" >&2; exit 1; }
echo "Raw VCF:  ${RAW_VCF}"

# ── Output directory ──────────────────────────────────────────────────────
REGENO_DIR="${OUTDIR}/A03b_population_regenotype"
mkdir -p "${REGENO_DIR}"
echo "Output:   ${REGENO_DIR}"
echo ""

# ── Run regenotyping ──────────────────────────────────────────────────────
echo "--- Running population_regenotype.py ---"

python3 "${REGENO_PY}" \
  --vcf "${RAW_VCF}" \
  --caller "${CALLER}" \
  --sv_type "${SV_LABEL}" \
  --outdir "${REGENO_DIR}" \
  --min_confident_reads 4 \
  --min_confident_gq 10 \
  --prior_floor 0.005 \
  --posterior_threshold 0.5 \
  --error_rate 0.01 \
  --write_posteriors \
  2>&1 | tee "${REGENO_DIR}/run.log"

echo ""
echo "--- Verifying outputs ---"
for f in \
  "population_regenotype.GT_matrix.tsv" \
  "population_regenotype.rescue_summary.tsv" \
  "population_regenotype.posteriors.tsv.gz"; do
  if [[ -f "${REGENO_DIR}/${f}" ]]; then
    echo "  OK: ${f} ($(wc -l < "${REGENO_DIR}/${f}" 2>/dev/null || echo '?') lines)"
  else
    echo "  MISSING: ${f}" >&2
  fi
done

# ── Compare original vs corrected carrier counts ─────────────────────────
echo ""
echo "--- Rescue summary ---"
if [[ -f "${REGENO_DIR}/population_regenotype.rescue_summary.tsv" ]]; then
  N_SITES=$(tail -n +2 "${REGENO_DIR}/population_regenotype.rescue_summary.tsv" | wc -l)
  N_RESCUED=$(awk -F'\t' 'NR>1 && $15>0' "${REGENO_DIR}/population_regenotype.rescue_summary.tsv" | wc -l)
  TOTAL_RESCUED=$(awk -F'\t' 'NR>1{s+=$15}END{print s}' "${REGENO_DIR}/population_regenotype.rescue_summary.tsv")
  echo "  Sites:           ${N_SITES}"
  echo "  Sites with rescue: ${N_RESCUED}"
  echo "  Total rescued:   ${TOTAL_RESCUED}"

  # Top 10 most-rescued sites
  echo ""
  echo "  Top 10 most-rescued sites:"
  tail -n +2 "${REGENO_DIR}/population_regenotype.rescue_summary.tsv" \
    | sort -t$'\t' -k15 -rn | head -10 \
    | awk -F'\t' '{printf "    %s  %s:%s  freq=%.3f  orig=%s→new=%s  rescued=%s (+%.0f%%)\n", $1,$2,$3,$5,$13,$14,$15,$16}'
fi

echo ""
echo "========================================"
echo "COMPLETE"
echo "========================================"
echo ""
echo "NEXT STEPS:"
echo "  1. Inspect rescue_summary.tsv for sites with large carrier increases"
echo "  2. Use population_regenotype.GT_matrix.tsv as input for phase_3_refine"
echo "  3. For BND module: run on both INV and BND, then re-run STEP06"
echo "     with corrected carrier sets for improved BND junction concordance"
