#!/usr/bin/env bash
# =============================================================================
# 00_module2b_config.sh — Central config for MODULE_2B Population Structure
#
# v9.0 REWIRED:
#   - Added REGISTRY_DIR + LOAD_BRIDGE paths
#   - Added SAMPLES_IND (canonical samples.ind path)
#   - Added ANCESTRY_CONFIG for cross-module wiring
#   - Palette unified with theme_systems_plate.R
#   - THIN_FINE defined (was missing in old config)
#   - Flat output directory convention
#
# This config can coexist with 00_inversion_config.sh and 00_ancestry_config.sh.
# All three share BASE but use module-prefixed variables for their outputs.
# =============================================================================

# ── Project root ─────────────────────────────────────────────────────────────
export BASE="${BASE:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"

# ── Reference genome ─────────────────────────────────────────────────────────
export REF="${BASE}/00-samples/fClaHyb_Gar_LG.fa"
export REF_FAI="${REF}.fai"

# ── Programs ─────────────────────────────────────────────────────────────────
export CONDA_ENV="assembly"
export RSCRIPT_BIN="${RSCRIPT_BIN:-/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript}"
export EVALADMIX_DIR="${BASE}/evalAdmix"
export EVALADMIX_BIN="${EVALADMIX_DIR}/evalAdmix"

# ── Genome ───────────────────────────────────────────────────────────────────
export GENOME_SIZE=963905721
export N_SAMPLES=226
export N_CHROMOSOMES=28

# ── Sample lists ─────────────────────────────────────────────────────────────
export BAMLIST="${BASE}/bamlist.pp.samechr.tlenP99.filtered.txt"
export SAMPLE_LIST="${BASE}/popstruct_thin/list_of_samples_one_per_line_same_bamfile_list.tsv"
export SAMPLES_IND="${BASE}/het_roh/01_inputs_check/samples.ind"
export PRUNED_LIST="${BASE}/popstruct_thin/06_relatedness/pruned_samples.txt"
export N_PRUNED=81

# ── Cross-module wiring (v9.0) ──────────────────────────────────────────────
export REGISTRY_DIR="${BASE}/sample_registry"
export LOAD_BRIDGE="${BASE}/inversion_codebase_v8.5/utils/load_bridge.R"
export ANCESTRY_CONFIG="${BASE}/unified_ancestry/00_ancestry_config.sh"
export INVERSION_CONFIG="${BASE}/inversion_codebase_v8.5/00_inversion_config.sh"

# ── Upstream inputs (from STEP1 biSNP/BEAGLE pipeline) ──────────────────────
export THIN_DIR="${BASE}/popstruct_thin"
export BEAGLE_MERGED_DIR="${THIN_DIR}/03_merged_beagle"
export BEAGLE_BYCHR_DIR="${THIN_DIR}/04_beagle_byRF_majmin"
export BEAGLE_DIR="${BEAGLE_BYCHR_DIR}"
export SITES_DIR="${THIN_DIR}/03_sites"
export RELATE_DIR="${THIN_DIR}/06_relatedness"

# ── Thinning panels ─────────────────────────────────────────────────────────
export THIN_DISTANCES=(200 500 1000)
export THIN_FINE=("${THIN_DISTANCES[@]}")
export THIN_PANELS=(500 1000)

# ── NGSadmix parameters ─────────────────────────────────────────────────────
export K_MIN=2
export K_MAX=20
export SEEDS=(1 2 3 4 5)
export NGSADMIX_MINMAF=0.05

# ── evalAdmix parameters ────────────────────────────────────────────────────
export EVALADMIX_NITS=5
export EVALADMIX_MISTOL=0.05
export EVALADMIX_MINMAF=0.05

# ── PCAngsd parameters ──────────────────────────────────────────────────────
export PCANGSD_EIG=6
export PCANGSD_MAF=0.05
export PCANGSD_ITER=250

# ── Relatedness thresholds ───────────────────────────────────────────────────
export THETA_DUP_MZ=0.354
export THETA_FIRST=0.177
export THETA_SECOND=0.0884
export THETA_THIRD=0.0442

# ── NAToRA parameters ───────────────────────────────────────────────────────
export NATORA_MIN_V=0.0
export NATORA_MAX_V=0.5
export NATORA_HEURISTIC_E=1
export NATORA_DUP_MZ_LOW=0.3540
export NATORA_DUP_MZ_HIGH=0.5000
export NATORA_FIRST_LOW=0.1770
export NATORA_FIRST_HIGH=0.3540
export NATORA_SECOND_LOW=0.0884
export NATORA_SECOND_HIGH=0.1770
export NATORA_THIRD_LOW=0.0442
export NATORA_THIRD_HIGH=0.0884

# ── Best-seed selection ─────────────────────────────────────────────────────
export PALETTE_NAME="catfish_ngsadmix_v1"
export BEST_SEED_RULE="highest_loglik_then_lowest_mean_abs_resid"

# =============================================================================
# OUTPUT DIRECTORY CONVENTION (flat)
# =============================================================================
export MODULE2B_RESULTS="${THIN_DIR}/structure_results"
export MODULE2B_EVALADMIX="${MODULE2B_RESULTS}/evaladmix"
export MODULE2B_PCANGSD="${MODULE2B_RESULTS}/pcangsd"
export MODULE2B_FIGURES="${MODULE2B_RESULTS}/figures"

build_run_tag() {
  local scope="${1:-wholegenome}"
  local thin="${2:-500}"
  local sample_set="${3:-all}"
  local n_samp="${4:-${N_SAMPLES}}"
  echo "${scope}_thin${thin}_${sample_set}${n_samp}"
}
export -f build_run_tag

# ── Theme file ───────────────────────────────────────────────────────────────
export THEME_FILE="${BASE}/unified_ancestry/utils/theme_systems_plate.R"

# ── Relatedness input files ──────────────────────────────────────────────────
export RELATEDNESS_RES="${RELATE_DIR}/catfish_${N_SAMPLES}_relatedness.res"
export NATORA_SUMMARY="${RELATE_DIR}/natora_cutoff_summary.tsv"

# ── Palette ──────────────────────────────────────────────────────────────────
export ANCESTRY_PALETTE=(
  "#3B6FA0" "#CF6839" "#C44E52" "#6BA08E" "#5A8F4A"
  "#C9A83E" "#8B6DAD" "#E8919C" "#5C7A3A" "#B07850"
  "#4A8C9F" "#9E6B8A" "#7A8B3C" "#C47A5E" "#5E7FAA"
  "#A06B4F" "#6B9E7A" "#B0855A" "#7E6EA0" "#A0856B"
)

# ── SLURM defaults ───────────────────────────────────────────────────────────
export SLURM_ACCOUNT="${SLURM_ACCOUNT:-lt200308}"
export SLURM_PARTITION="${SLURM_PARTITION:-compute}"
export SLURM_DEFAULT_TIME="1-00:00:00"
export SLURM_DEFAULT_MEM="237G"
export SLURM_DEFAULT_CPUS=80

# ── Convenience ──────────────────────────────────────────────────────────────
timestamp(){ date '+%F %T'; }
export -f timestamp

m2b_log()  { echo "[$(timestamp)] [module2b] $*"; }
m2b_die()  { echo "[$(timestamp)] [module2b] FATAL: $*" >&2; exit 1; }
m2b_check_file() { [[ -s "$1" ]] || m2b_die "Missing or empty: $1"; }
m2b_init_dirs() {
  mkdir -p "${MODULE2B_RESULTS}" "${MODULE2B_EVALADMIX}" \
           "${MODULE2B_PCANGSD}" "${MODULE2B_FIGURES}" \
           "${REGISTRY_DIR}/groups" "${REGISTRY_DIR}/backups"
}
export -f m2b_log m2b_die m2b_check_file m2b_init_dirs
