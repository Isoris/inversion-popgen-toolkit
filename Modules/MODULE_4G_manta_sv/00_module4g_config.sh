#!/usr/bin/env bash
# =============================================================================
# 00_manta_config.sh — Central configuration for Manta SV pipeline
# =============================================================================
#
# Manta calls ALL SV types in a single run (DEL, DUP, INS, INV/BND).
# Post-processing splits results into per-type catalogs.
#
# KEY DIFFERENCES FROM DELLY:
#   - Manta runs per-sample (configManta.py + runWorkflow.py), no merge/regenotype
#   - Uses custom config ini with minCandidateVariantSize = 50
#   - Inversions are encoded as BND pairs with INV3/INV5 tags; converted via
#     convertInversion.py from Manta's libexec
#   - Insertions split into two classes:
#       * Small INS: fully assembled, SVLEN present, 50 bp ≤ SVLEN ≤ ~200 bp
#       * Large INS: incompletely assembled, LEFT_SVINSSEQ / RIGHT_SVINSSEQ
#         present, no SVLEN (unknown total length)
#
# HOW TO SOURCE THIS (in SLURM scripts):
#   SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
#   CONFIG="${SCRIPT_DIR}/00_manta_config.sh"
#   [[ -f "${CONFIG}" ]] || { echo "Missing config" >&2; exit 1; }
#   set -a
#   source "${CONFIG}"
#   set +a
# =============================================================================

# ── Manta installation ──────────────────────────────────────────────────────
# Manta installed in manta_py2 environment
MANTA_INSTALL="/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/manta_py2"
MANTA_CONFIG="${MANTA_INSTALL}/bin/configManta.py"
MANTA_CONVERT_INV="${MANTA_INSTALL}/bin/convertInversion.py"

# ── Project root ────────────────────────────────────────────────────────────
MANTA_PROJECT="${MANTA_PROJECT:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"

# ── Reference genome ────────────────────────────────────────────────────────
REF="${MANTA_PROJECT}/00-samples/fClaHyb_Gar_LG.fa"
REF_FAI="${REF}.fai"

# ── BAM manifest ────────────────────────────────────────────────────────────
# TSV format: sample_id <tab> unfiltered_bam <tab> filtered_bam
# Column 1 = sample ID
# Column 2 = unfiltered/raw BAM path → used for Manta (markdup copies)
# Column 3 = P99/TLEN/MAPQ30 filtered BAM → IGNORE for Manta
BAM_MANIFEST="${MANTA_PROJECT}/pa_roary_results/00_manifests/sample_bam_minimap2_vs_P99TLENMAPQ30.tsv"

# ── Duplicate-marked BAM directory ─────────────────────────────────────────
# Reuse DELLY markdup BAMs
DIR_MARKDUP="${MANTA_PROJECT}/MODULE_4B_DEL_Delly/00_markdup"

# ── Sample lists ────────────────────────────────────────────────────────────
SAMPLES_ALL="${MANTA_PROJECT}/MODULE_4G_ALL_Manta/samples_all_226.txt"
SAMPLES_UNRELATED="${MANTA_PROJECT}/MODULE_4G_ALL_Manta/samples_unrelated_81.txt"
NATORA_KEEP="${MANTA_PROJECT}/popstruct_thin/05_ngsrelate/catfish_first_degree_pairwise_toKeep.txt"

# ── Call regions BED (bgzipped + tabix-indexed) ────────────────────────────
CALL_REGIONS_BED_GZ="${MANTA_PROJECT}/MODULE_4G_ALL_Manta/call_regions.bed.gz"

# ── Exclude BED used to build call regions ─────────────────────────────────
EXCL_BED="${MANTA_PROJECT}/MODULE_4G_ALL_Manta/exclude.minimal.bed"

# ── Annotation files ───────────────────────────────────────────────────────
GFF3="${MANTA_PROJECT}/00-samples/fClaHyb_Gar_LG.from_CGAR.gff3_polished.sorted.gff3"
ANNOT_DIR="${MANTA_PROJECT}/pa_roary_results/annotation_flag_beds"
GENE_BED="${ANNOT_DIR}/features/GENE.bed"
EXON_BED="${ANNOT_DIR}/features/EXON.bed"
CDS_BED="${ANNOT_DIR}/features/CDS.bed"
REPEAT_BED="${MANTA_PROJECT}/fClaHyb_Gar_LG.mask_regions.softacgt.renamed.bed"

# ── Custom Manta config ini ────────────────────────────────────────────────
MANTA_CUSTOM_INI="${MANTA_PROJECT}/MODULE_4G_ALL_Manta/configManta.custom.ini"

# ── Output directories ─────────────────────────────────────────────────────
OUTDIR="${MANTA_PROJECT}/MODULE_4G_ALL_Manta"
DIR_PERSAMPLE="${OUTDIR}/01_per_sample"
DIR_CONVERTED="${OUTDIR}/01b_per_sample_converted"
DIR_MERGED="${OUTDIR}/02_merged_cohort"
DIR_SUBSET81="${OUTDIR}/03_subset_81"
DIR_SPLIT="${OUTDIR}/04_split_by_type"
DIR_FINAL="${OUTDIR}/05_final_catalogs"
DIR_ANNOT="${OUTDIR}/06_annotation"
DIR_DEPTH="${OUTDIR}/07_depth_support"
DIR_SUMMARY="${OUTDIR}/08_summary"
DIR_PLOTS="${OUTDIR}/09_plots"
DIR_LOGS="${OUTDIR}/logs"

# ── Performance ─────────────────────────────────────────────────────────────
THREADS="${SLURM_CPUS_PER_TASK:-127}"
MANTA_THREADS_PER_CALL=4
MANTA_PARALLEL=30

# ── Depth support (mosdepth) ───────────────────────────────────────────────
DEPTH_WINDOW=500
DEPTH_MAPQ=30
DEPTH_THREADS=4

# ── Strict final catalog filters ───────────────────────────────────────────
STRICT_REQUIRE_PASS=1
STRICT_MIN_QUAL=20

# ── Ancestry / admixture (for plotting) ────────────────────────────────────
PA_BEST_SEED_TABLE="${MANTA_PROJECT}/popstruct_thin/05_ngsadmix_global/best_seed_by_K.tsv"
PA_NGSADMIX_DIR="${MANTA_PROJECT}/popstruct_thin/05_ngsadmix_global/runs_thin500"

# =============================================================================
# Helper functions
# =============================================================================
mv_timestamp() { date '+%F %T'; }
mv_log() { echo "[$(mv_timestamp)] [MANTA] $*"; }
mv_err() { echo "[$(mv_timestamp)] [MANTA] [ERROR] $*" >&2; }
mv_die() { mv_err "$@"; exit 1; }

mv_check_file() {
  local f="$1" label="${2:-file}"
  [[ -e "$f" ]] || mv_die "Missing ${label}: ${f}"
}

mv_check_cmd() {
  command -v "$1" &>/dev/null || mv_die "Command not found: $1"
}

mv_init_dirs() {
  mkdir -p \
    "$OUTDIR" \
    "$DIR_PERSAMPLE" "$DIR_CONVERTED" "$DIR_MERGED" "$DIR_SUBSET81" "$DIR_SPLIT" \
    "$DIR_FINAL" "$DIR_ANNOT" "$DIR_DEPTH" "$DIR_SUMMARY" \
    "$DIR_PLOTS" "$DIR_LOGS"
}

mv_read_manifest() {
  [[ -s "$BAM_MANIFEST" ]] || mv_die "BAM manifest not found or empty: $BAM_MANIFEST"

  MANIFEST_SIDS=()
  MANIFEST_RAW_BAMS=()

  {
    read -r _ _ _
    while IFS=$'\t' read -r sid raw_bam _ignore; do
      [[ -z "$sid" ]] && continue
      MANIFEST_SIDS+=("$sid")
      MANIFEST_RAW_BAMS+=("$raw_bam")
    done
  } < "$BAM_MANIFEST"

  MANIFEST_N_SAMPLES=${#MANIFEST_SIDS[@]}
  [[ $MANIFEST_N_SAMPLES -gt 0 ]] || mv_die "Zero samples loaded from manifest: $BAM_MANIFEST"

  mv_log "Manifest loaded: ${MANIFEST_N_SAMPLES} samples from ${BAM_MANIFEST}"
}

mv_markdup_bam() {
  local sid="$1"
  echo "${DIR_MARKDUP}/${sid}.markdup.bam"
}
