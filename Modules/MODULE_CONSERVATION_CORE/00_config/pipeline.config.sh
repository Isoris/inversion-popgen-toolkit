#!/usr/bin/env bash
# =============================================================================
# 00_conservation_config.sh — Central configuration for MODULE_CONSERVATION
# CORE-ONLY VERSION: No cross-species alignment. Three scoring tools only:
#   bcftools csq + SIFT4G + VESM + splice module + hatchery risk
#
# Source this at the top of every shell script.
# =============================================================================

# ── Project root ─────────────────────────────────────────────────────────────
BASE="${BASE:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"
MODCONS="${BASE}/MODULE_CONSERVATION"

# ── SLURM defaults ───────────────────────────────────────────────────────────
SLURM_ACCOUNT="lt200308"
SLURM_PARTITION="compute"

# ── Reference genome (Gar haplotype, primary) ────────────────────────────────
REF_FASTA="${BASE}/00-samples/fClaHyb_Gar_LG.fa"
REF_FAI="${BASE}/00-samples/fClaHyb_Gar_LG.fa.fai"
REF_GTF="${BASE}/00-samples/fClaHyb_Gar_LG.gtf"
REF_GFF3="${MODCONS}/01_reference/fClaHyb_Gar_LG.gff3"

# ── Derived annotation (built by STEP_00_setup.sh) ──────────────────────────
REF_CDS="${MODCONS}/01_reference/fClaHyb_Gar_LG.cds.fa"
REF_PEP="${MODCONS}/01_reference/fClaHyb_Gar_LG.pep.fa"
REF_CDNA="${MODCONS}/01_reference/fClaHyb_Gar_LG.cdna.fa"
GENE_MAP="${MODCONS}/01_reference/gene_transcript_protein_map.tsv"
CANONICAL_TX="${MODCONS}/01_reference/canonical_transcripts.txt"

# ── Chromosomes ──────────────────────────────────────────────────────────────
CHROMS=(C_gar_LG{01..28})
N_CHROMS=28

# ── Sample manifests ─────────────────────────────────────────────────────────
BAM_MANIFEST="${BASE}/pa_roary_results/00_manifests/sample_bam_minimap2_vs_P99TLENMAPQ30.tsv"
SAMPLES_ALL="${BASE}/pa_roary_results/00_manifests/sample_list_226.txt"
SAMPLES_UNRELATED="${BASE}/pa_roary_results/00_manifests/samples_unrelated_81.txt"
N_SAMPLES_ALL=226
N_SAMPLES_UNRELATED=81

# ── Upstream modules (variant sources) ───────────────────────────────────────
CLAIR3_DIR="${BASE}/MODULE_4A_SNP_INDEL50_Clair3"
ROH_DIR="${BASE}/MODULE_3_HetROH"
NGSADMIX_DIR="${BASE}/pa_roary_results/03_ngsadmix"

# ── Programs ─────────────────────────────────────────────────────────────────
PROG="/scratch/lt200308-agbsci/13-programs"
RSCRIPT="${PROG}/mambaforge/envs/assembly/bin/Rscript"
CONDA_BASE="/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge"

# ── Tools (CORE only: SnpEff + SIFT4G) ──────────────────────────────────────
SNPEFF_DIR="${PROG}/snpEff"
SNPEFF_JAR="${SNPEFF_DIR}/snpEff.jar"
SNPSIFT_JAR="${SNPEFF_DIR}/SnpSift.jar"
SIFT4G_DIR="${PROG}/sift4g"
SIFT4G_BIN="${SIFT4G_DIR}/bin/sift4g"

# ── Results subdirectories (CORE only) ───────────────────────────────────────
DIR_REF="${MODCONS}/01_reference"
DIR_ANNO="${MODCONS}/02_annotation"
DIR_VARIANTS="${MODCONS}/03_variants"
DIR_SNPEFF="${MODCONS}/04_snpeff"
DIR_SIFT4G="${MODCONS}/05_sift4g"
DIR_VESM="${MODCONS}/05B_vesm"
DIR_MERGED="${MODCONS}/16_merged_variant_tables"
DIR_BURDEN="${MODCONS}/17_burden_tables"
DIR_ROH="${MODCONS}/18_roh_overlap"
DIR_RESULTS="${MODCONS}/20_results"
DIR_LOGS="${MODCONS}/21_logs"

# ── Scoring config ───────────────────────────────────────────────────────────
SCORE_CONFIG="${MODCONS}/00_config/scoring_weights.tsv"

# ── Scoring thresholds (configurable) ────────────────────────────────────────
# Max possible del_evidence = 8(csq) + 4(SIFT) + 5(VESM) + 5(splice) = 22
# Max possible hatchery     = 3(MAF) + 3(hom) + 3(ROH) = 9
# Max total                 = 31
# Threshold rationale:
#   Class A ≥ 16 : severe consequence + damaging (SIFT or VESM) + high exposure
#   Class B ≥ 10 : protein-changing + at least one damaging predictor + some exposure
#   Class C ≥ 5  : some evidence, possibly ambiguous
#   Class D < 5  : low priority
PRIORITY_THRESHOLD_A=16
PRIORITY_THRESHOLD_B=10
PRIORITY_THRESHOLD_C=5

# ── Minimum genotype depth for burden counts (configurable) ──────────────────
MIN_GT_DP=3

# ── Helper functions ─────────────────────────────────────────────────────────
cons_timestamp() { date '+%F %T'; }
cons_log()  { echo "[$(cons_timestamp)] [CONS] $*"; }
cons_err()  { echo "[$(cons_timestamp)] [CONS] [ERROR] $*" >&2; }
cons_die()  { cons_err "$@"; exit 1; }
cons_check_file() {
    local f="$1" label="${2:-file}"
    [[ -e "$f" ]] || cons_die "Missing ${label}: ${f}"
}
cons_init_dirs() {
    mkdir -p \
        "${MODCONS}/slurm_logs" \
        "${DIR_REF}" "${DIR_ANNO}" "${DIR_VARIANTS}" \
        "${DIR_SNPEFF}" "${DIR_SIFT4G}" "${DIR_VESM}" \
        "${DIR_MERGED}" "${DIR_BURDEN}" "${DIR_ROH}" \
        "${DIR_RESULTS}" "${DIR_LOGS}"
}

# ── Export everything ────────────────────────────────────────────────────────
export BASE MODCONS SLURM_ACCOUNT SLURM_PARTITION
export REF_FASTA REF_FAI REF_GTF REF_GFF3 REF_CDS REF_PEP REF_CDNA GENE_MAP CANONICAL_TX
export CHROMS N_CHROMS
export BAM_MANIFEST SAMPLES_ALL SAMPLES_UNRELATED N_SAMPLES_ALL N_SAMPLES_UNRELATED
export CLAIR3_DIR ROH_DIR
export PROG RSCRIPT CONDA_BASE
export SNPEFF_DIR SNPEFF_JAR SNPSIFT_JAR SIFT4G_DIR SIFT4G_BIN
export DIR_REF DIR_ANNO DIR_VARIANTS DIR_SNPEFF DIR_SIFT4G DIR_VESM
export DIR_MERGED DIR_BURDEN DIR_ROH DIR_RESULTS DIR_LOGS
export SCORE_CONFIG
export PRIORITY_THRESHOLD_A PRIORITY_THRESHOLD_B PRIORITY_THRESHOLD_C
export MIN_GT_DP
export -f cons_timestamp cons_log cons_err cons_die cons_check_file cons_init_dirs
