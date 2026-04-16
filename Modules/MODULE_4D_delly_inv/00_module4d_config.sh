#!/usr/bin/env bash
# =============================================================================
# 00_delly_config.sh — Central configuration for DELLY INV pipeline
# =============================================================================
# Adapted from DEL pipeline. SV type: Inversions
# DELLY call type: -t INV
# =============================================================================

# ── Module loads ────────────────────────────────────────────────────────────
module load HTSlib/1.17-cpeGNU-23.03
module load Boost/1.81.0-cpeGNU-23.03

# ── DELLY binary ────────────────────────────────────────────────────────────
DELLY_BIN="/project/lt200308-agbsci/01-catfish_assembly/delly/src/delly"

# ── Project root ────────────────────────────────────────────────────────────
BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
DELLY_PROJECT="${BASE}"

# ── Reference genome ────────────────────────────────────────────────────────
REF="${BASE}/00-samples/fClaHyb_Gar_LG.fa"
REF_FAI="${REF}.fai"

# ── SV type ─────────────────────────────────────────────────────────────────
SV_TYPE="INV"
DELLY_CALL_TYPE="INV"

# ── Reuse shared inputs from DEL pipeline ───────────────────────────────────
EXCL_BED="${BASE}/delly_sv/exclude.minimal.bed"
DIR_MARKDUP="${BASE}/delly_sv/00_markdup"
BAM_MANIFEST="${BASE}/pa_roary_results/00_manifests/sample_bam_minimap2_vs_P99TLENMAPQ30.tsv"

# ── Sample lists (reuse from DEL pipeline) ─────────────────────────────────
SAMPLES_ALL="${BASE}/delly_sv/samples_all_226.txt"
SAMPLES_UNRELATED="${BASE}/delly_sv/samples_unrelated_81.txt"
NATORA_KEEP="${BASE}/popstruct_thin/05_ngsrelate/catfish_first_degree_pairwise_toKeep.txt"

# ── Annotation files ────────────────────────────────────────────────────────
GFF3="${BASE}/00-samples/fClaHyb_Gar_LG.from_CGAR.gff3_polished.sorted.gff3"
ANNOT_DIR="${BASE}/pa_roary_results/annotation_flag_beds"
GENE_BED="${ANNOT_DIR}/features/GENE.bed"
EXON_BED="${ANNOT_DIR}/features/EXON.bed"
CDS_BED="${ANNOT_DIR}/features/CDS.bed"
REPEAT_BED="${BASE}/fClaHyb_Gar_LG.mask_regions.softacgt.renamed.bed"

# ── Output directories ──────────────────────────────────────────────────────
OUTDIR="${BASE}/MODULE_4D_INV_Delly"
DIR_DISC="${OUTDIR}/01_discovery"
DIR_SITES="${OUTDIR}/02_merged_sites"
DIR_GENO="${OUTDIR}/03_genotyped"
DIR_MERGED="${OUTDIR}/04_merged_cohort"
DIR_SUBSET81="${OUTDIR}/05_subset_81"
DIR_FILTERED="${OUTDIR}/06_germline_filtered"
DIR_FINAL="${OUTDIR}/07_final_catalogs"
DIR_ANNOT="${OUTDIR}/08_annotation"
DIR_DEPTH="${OUTDIR}/09_depth_support"
DIR_MATDIST="${OUTDIR}/10_mate_distance_qc"
DIR_SUMMARY="${OUTDIR}/11_summary"
DIR_PLOTS="${OUTDIR}/12_plots"
DIR_LOGS="${OUTDIR}/logs"

# ── Performance ─────────────────────────────────────────────────────────────
THREADS="${SLURM_CPUS_PER_TASK:-127}"
DELLY_THREADS_PER_CALL=4
DELLY_PARALLEL=30

# ── Strict filter thresholds ────────────────────────────────────────────────
STRICT_REQUIRE_PASS=1
STRICT_REQUIRE_PRECISE=1
STRICT_MIN_QUAL=300
STRICT_MIN_PE=3

# ── Mate distance / SVLEN thresholds ────────────────────────────────────────
MATE_WARN_KB=20
MATE_SUSPICIOUS_KB=50
MATE_EXTREME_KB=100

# ── Ancestry / admixture (for plotting) ────────────────────────────────────
PA_BEST_SEED_TABLE="${BASE}/popstruct_thin/05_ngsadmix_global/best_seed_by_K.tsv"
PA_NGSADMIX_DIR="${BASE}/popstruct_thin/05_ngsadmix_global/runs_thin500"

# =============================================================================
# Helper functions
# =============================================================================
dv_timestamp() { date '+%F %T'; }
dv_log() { echo "[$(dv_timestamp)] [DELLY-INV] $*"; }
dv_err() { echo "[$(dv_timestamp)] [DELLY-INV] [ERROR] $*" >&2; }
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
    "$DIR_ANNOT" "$DIR_DEPTH" "$DIR_MATDIST" "$DIR_SUMMARY" \
    "$DIR_PLOTS" "$DIR_LOGS"
}
