#!/bin/bash
# =============================================================================
# 00_config.sh — MODULE_QC_ShelfDiagnosis
# =============================================================================
# Centralized paths, sample list, tool bins. Source this at the top of every
# script in this module. Override any variable via env before sourcing.
# =============================================================================

# ---- Project root ----
: "${BASE:=/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"

# ---- Inputs ----
: "${BEAGLE_DIR:=${BASE}/inversion_localpca_v7/02_snps_beagle}"
: "${PRECOMP_DIR:=${BASE}/inversion_localpca_v7/06_mds_candidates/snake_regions_multiscale/precomp}"
: "${BAM_DIR:=${BASE}/MODULE_1_Reads}"        # adjust if BAMs live elsewhere
: "${BAM_PATTERN:=*.bam}"
: "${MOSDEPTH_EXISTING:=${BASE}/het_roh/mosdepth_output}"  # if already run

# ---- Heterozygosity / theta ----
: "${HET_DIR:=${BASE}/het_roh}"
: "${THETA_MAIN_DIR:=${HET_DIR}/02_heterozygosity/03_theta}"
: "${THETA_MULTI_DIR:=${THETA_MAIN_DIR}/multiscale}"

# ---- Unified ancestry (Engine B) ----
: "${UNIFIED_ANCESTRY_DIR:=${BASE}/unified_ancestry}"
: "${LOCAL_Q_DIR:=${UNIFIED_ANCESTRY_DIR}/local_Q}"

# ---- Engine F (region_popstats — Fst/dXY/theta/Tajima in one C binary) ----
: "${REGION_POPSTATS_BIN:=${UNIFIED_ANCESTRY_DIR}/engines/fst_dxy/region_popstats}"
: "${SAMPLE_LIST_POPSTATS:=${BASE}/het_roh/01_inputs_check/samples.ind}"

# ---- Outputs ----
: "${QC_OUT:=${BASE}/inversion_modules/phase_qc_shelf/results}"
: "${QC_TRACKS:=${QC_OUT}/tracks}"       # per-chrom binned tsv tracks
: "${QC_FIGS:=${QC_OUT}/figures}"        # diagnostic PDFs
: "${QC_LOGS:=${QC_OUT}/logs}"

mkdir -p "${QC_TRACKS}" "${QC_FIGS}" "${QC_LOGS}"

# ---- Sample list ----
: "${SAMPLE_LIST:=${BASE}/het_roh/01_inputs_check/samples.ind}"

# ---- Binning ----
: "${BIN_MB:=0.05}"                      # Mb per genome-wide track bin

# ---- Tool binaries ----
: "${RSCRIPT_BIN:=Rscript}"
: "${MOSDEPTH_BIN:=mosdepth}"
: "${SAMTOOLS_BIN:=samtools}"
: "${TABIX_BIN:=tabix}"
: "${BGZIP_BIN:=bgzip}"
: "${ZCAT_BIN:=zcat}"

# ---- Logging helpers ----
qc_log() { printf '[%s] %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$*" >&2; }
qc_die() { qc_log "FATAL: $*"; exit 1; }

# ---- Sanity check on source ----
_check_qc_config() {
  [[ -d "${BEAGLE_DIR}" ]] || qc_log "WARNING: BEAGLE_DIR not found: ${BEAGLE_DIR}"
  [[ -d "${PRECOMP_DIR}" ]] || qc_log "WARNING: PRECOMP_DIR not found: ${PRECOMP_DIR}"
  [[ -r "${SAMPLE_LIST}" ]] || qc_log "WARNING: SAMPLE_LIST not found: ${SAMPLE_LIST}"
}
_check_qc_config
