#!/usr/bin/env bash
# =============================================================================
# 00_extra_plots_config.sh — config for phase_12_cargo/extra_plots
#
# Sources the cargo config (which sources the main inversion config) and adds
# extras-specific output paths.
# =============================================================================
set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../00_cargo_config.sh"

# ── Output root for extras ───────────────────────────────────────────────────
export EXTRAS_DIR="${EXTRAS_DIR:-${CARGO_DIR}/extras}"
export EXTRAS_FIG_DIR="${EXTRAS_FIG_DIR:-${EXTRAS_DIR}/figures}"
export EXTRAS_TBL_DIR="${EXTRAS_TBL_DIR:-${EXTRAS_DIR}/tables}"
export EXTRAS_LOGS_DIR="${EXTRAS_LOGS_DIR:-${EXTRAS_DIR}/logs}"
mkdir -p "${EXTRAS_FIG_DIR}" "${EXTRAS_TBL_DIR}" "${EXTRAS_LOGS_DIR}"

# ── Optional inputs we'll attempt to find ────────────────────────────────────
# ROH bed (built by MODULE_3 het_roh)
export ROH_BED="${ROH_BED:-${HETDIR:-${BASE}/het_roh}/04_roh/all_samples_roh.bed.gz}"
# DELLY SV calls (per type)
export DELLY_DEL_VCF="${DELLY_DEL_VCF:-${BASE}/sv_calling/delly_del/delly_del.merged.vcf.gz}"
export DELLY_DUP_VCF="${DELLY_DUP_VCF:-${BASE}/sv_calling/delly_dup/delly_dup.merged.vcf.gz}"
export DELLY_INV_VCF="${DELLY_INV_VCF:-${BASE}/sv_calling/delly_inv/delly_inv.merged.vcf.gz}"
export DELLY_BND_VCF="${DELLY_BND_VCF:-${BASE}/sv_calling/delly_bnd/delly_bnd.merged.vcf.gz}"
# Manta combined output
export MANTA_VCF="${MANTA_VCF:-${BASE}/sv_calling/manta/manta.cohort.vcf.gz}"
# NGSadmix Q at canonical K (groups for "across-group" plots)
export NGSADMIX_Q_FILE="${NGSADMIX_Q_FILE:-${NGSADMIX_DIR}/runs_thin500/K${CANONICAL_K:-8}/best.qopt}"

# ── Defaults ─────────────────────────────────────────────────────────────────
export EXTRAS_TOPN_RECURRENT="${EXTRAS_TOPN_RECURRENT:-25}"
export EXTRAS_TOPN_SHARING="${EXTRAS_TOPN_SHARING:-50}"
# For PLOT_07 group exclusivity: only show variants exclusive to a single group
export EXTRAS_EXCLUSIVITY_THRESHOLD="${EXTRAS_EXCLUSIVITY_THRESHOLD:-0.9}"

echo "[extras_config] EXTRAS_DIR=${EXTRAS_DIR}"
