#!/usr/bin/env bash
# =============================================================================
# 00_delly_config.sh — Central configuration for DELLY DUP pipeline
# =============================================================================

module load HTSlib/1.17-cpeGNU-23.03
module load Boost/1.81.0-cpeGNU-23.03

DELLY_BIN="/project/lt200308-agbsci/01-catfish_assembly/delly/src/delly"

DELLY_PROJECT="${DELLY_PROJECT:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"

REF="${DELLY_PROJECT}/00-samples/fClaHyb_Gar_LG.fa"
REF_FAI="${REF}.fai"

# Reuse old DEL resources
EXCL_BED="${DELLY_PROJECT}/delly_sv/exclude.minimal.bed"
DIR_MARKDUP="${DELLY_PROJECT}/delly_sv/00_markdup"
SAMPLES_ALL="${DELLY_PROJECT}/delly_sv/samples_all_226.txt"
SAMPLES_UNRELATED="${DELLY_PROJECT}/delly_sv/samples_unrelated_81.txt"

# This pipeline's own output directory
OUTDIR="${DELLY_PROJECT}/delly_sv_DUP"

DIR_DISC="${OUTDIR}/01_discovery"
DIR_SITES="${OUTDIR}/02_merged_sites"
DIR_GENO="${OUTDIR}/03_genotyped"
DIR_MERGED="${OUTDIR}/04_merged_cohort"
DIR_SUBSET81="${OUTDIR}/05_subset_81"
DIR_FILTERED="${OUTDIR}/06_germline_filtered"
DIR_FINAL="${OUTDIR}/07_final_catalogs"
DIR_LOGS="${OUTDIR}/logs"

THREADS="${SLURM_CPUS_PER_TASK:-127}"
DELLY_THREADS_PER_CALL=4
DELLY_PARALLEL=30

# Variant type
SVTYPE="DUP"

# Conservative first-pass strict filter
STRICT_REQUIRE_PASS=1
STRICT_REQUIRE_PRECISE=1
STRICT_MIN_QUAL=500
STRICT_MIN_PE=5

#STRICT_REQUIRE_PASS=1
#STRICT_REQUIRE_PRECISE=1
#STRICT_MIN_QUAL=300
#STRICT_MIN_PE=3

dv_timestamp() { date '+%F %T'; }
dv_log() { echo "[$(dv_timestamp)] [DELLY-${SVTYPE}] $*"; }
dv_err() { echo "[$(dv_timestamp)] [DELLY-${SVTYPE}] [ERROR] $*" >&2; }
dv_die() { dv_err "$@"; exit 1; }

dv_check_file() {
  local f="$1" label="${2:-file}"
  [[ -e "$f" ]] || dv_die "Missing ${label}: ${f}"
}

dv_check_cmd() {
  command -v "$1" &>/dev/null || dv_die "Command not found: $1"
}

dv_init_dirs() {
  mkdir -p \
    "$OUTDIR" \
    "$DIR_DISC" "$DIR_SITES" "$DIR_GENO" \
    "$DIR_MERGED" "$DIR_SUBSET81" "$DIR_FILTERED" "$DIR_FINAL" \
    "$DIR_LOGS"
}
