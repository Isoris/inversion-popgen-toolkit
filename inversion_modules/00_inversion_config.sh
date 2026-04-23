#!/usr/bin/env bash
# =============================================================================
# 00_inversion_config.sh — Central configuration for inversion pipeline
#
# Source this file at the top of every shell script and SLURM launcher.
# All paths resolve from BASE. Nothing is hardcoded in individual scripts.
#
# Usage:
#   source "$(dirname "${BASH_SOURCE[0]}")/00_inversion_config.sh"
#   — or —
#   source /path/to/00_inversion_config.sh
#
# =============================================================================
# 2026-04-17 rewrite (chat 15):
#   - Flattened codebase layout. The old CODEBASE=inversion_codebase_v8.5
#     pointer is kept as a fallback, but auto-discovery prefers the flattened
#     inversion_modules/ layout at BASE root.
#   - Added LOCAL_Q_DIR, STATS_CACHE_DIR, CANONICAL_K, K_SWEEP, SV_PRIOR_DIR,
#     GHSL_DIR. All data dirs live OUTSIDE the code tree on scratch.
#   - Auto-discovery fallbacks for LOAD_BRIDGE / UTILS_DIR so the same config
#     works under both the flattened and legacy layouts during the transition.
#   - NGSADMIX_DIR now matches what the ancestry config uses
#     (05_ngsadmix_global/runs_thin500 via BEST_SEED_BY_K).
#   - SV_PRIOR_DIR / GHSL_DIR declared here so inversion launchers don't
#     need to source the ancestry config.
# =============================================================================

# ── Project root ─────────────────────────────────────────────────────────────
BASE="${BASE:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"

# ── Results root (where pipeline outputs live) ───────────────────────────────
INVDIR="${INVDIR:-${BASE}/inversion_localpca_v7}"

# ── Code root — auto-discovery ──────────────────────────────────────────────
# Prefer flattened layout; fall back to legacy v8.5.
if [[ -z "${CODEBASE:-}" ]]; then
  for _try in \
      "${BASE}" \
      "${BASE}/inversion_modules" \
      "${BASE}/inversion_codebase_v8.5"; do
    if [[ -d "${_try}/inversion_modules" || -d "${_try}/phase_2_discovery" ]]; then
      CODEBASE="$_try"
      break
    fi
  done
  : "${CODEBASE:=${BASE}/inversion_modules}"
  unset _try
fi
export CODEBASE

# Legacy module-split aliases — preserved for launchers that still reference
# them. Will point at non-existent paths under the flattened layout; that's
# fine as long as nothing under inversion_modules/ references them.
DISCOVERYDIR="${CODEBASE}/MODULE_5A2_Discovery_Core"
FOLLOWUPDIR="${CODEBASE}/MODULE_5B_Inversion_Followup"
LDDIR="${CODEBASE}/MODULE_5C_Inversion_LD"
FSTDIR="${CODEBASE}/MODULE_5D_Inversion_FST"
# HOBSDIR removed 2026-04-24: MODULE_5E archived to
# inversion_modules/_archive_superseded/MODULE_5E_Inversion_HOBS_superseded_by_Q07b/
# (superseded by 4b_qc_triage/STEP_Q07b + Q07c per-group Hobs).
DISCOVERY1DIR="${CODEBASE}/MODULE_5A1_Discovery_Inputs"
DISCOVERY2DIR="${CODEBASE}/MODULE_5A2_Discovery_Core"
DISCOVERY3DIR="${CODEBASE}/MODULE_5A3_Discovery_Postprocessing"

# Current phase-based layout (flattened). Preferred over the MODULE_5A/B/C/D/E
# aliases for any new code.
PHASE2_DIR="${CODEBASE}/inversion_modules/phase_2_discovery"
[[ -d "$PHASE2_DIR" ]] || PHASE2_DIR="${CODEBASE}/phase_2_discovery"
PHASE3_DIR="${CODEBASE}/inversion_modules/phase_3_refine"
[[ -d "$PHASE3_DIR" ]] || PHASE3_DIR="${CODEBASE}/phase_3_refine"
PHASE4_DIR="${CODEBASE}/inversion_modules/phase_4_postprocessing"
[[ -d "$PHASE4_DIR" ]] || PHASE4_DIR="${CODEBASE}/phase_4_postprocessing"
export PHASE2_DIR PHASE3_DIR PHASE4_DIR

# ── Rscript binary ──────────────────────────────────────────────────────────
RSCRIPT_BIN="${RSCRIPT_BIN:-/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript}"

# ── Stable external inputs ───────────────────────────────────────────────────
HETDIR="${BASE}/het_roh"
REF="${BASE}/00-samples/fClaHyb_Gar_LG.fa"
REF_FAI="${REF}.fai"
BAMLIST="${HETDIR}/01_inputs_check/bamlist_qcpass.txt"
SAMPLES_IND="${HETDIR}/01_inputs_check/samples.ind"

# ── Utils + bridge wiring — auto-discovery ─────────────────────────────────
if [[ -z "${UTILS_DIR:-}" ]]; then
  for _try in \
      "${BASE}/utils" \
      "${CODEBASE}/utils" \
      "${BASE}/inversion_modules/utils" \
      "${BASE}/inversion_codebase_v8.5/utils"; do
    if [[ -f "${_try}/load_bridge.R" ]]; then
      UTILS_DIR="$_try"; break
    fi
  done
  : "${UTILS_DIR:=${BASE}/utils}"
  unset _try
fi
export UTILS_DIR
export LOAD_BRIDGE="${LOAD_BRIDGE:-${UTILS_DIR}/load_bridge.R}"
SAMPLE_MAP_R="${UTILS_DIR}/sample_map.R"
SAMPLE_REGISTRY_R="${UTILS_DIR}/sample_registry.R"
THEME_FILE="${THEME_FILE:-${BASE}/unified_ancestry/utils/theme_systems_plate.R}"

# ── Sample registry ─────────────────────────────────────────────────────────
export REGISTRY_DIR="${REGISTRY_DIR:-${BASE}/sample_registry}"

# ── Registries tree (chat 15, rewired chat 16) ─────────────────────────────
export REGISTRIES_DATA_DIR="${REGISTRIES_DATA_DIR:-${BASE}/registries/data}"
# Results registry (chat 16, replaces chat-15 stats_cache). First-class
# fourth registry persisting FST/dxy tracks, per-candidate Q/F matrices,
# and per-interval ancestry summaries. See registries/DATABASE_DESIGN.md.
# A legacy STATS_CACHE_DIR export is preserved for scripts still using the
# chat-15 name; both resolve to the same directory so only one physical
# store exists.
export RESULTS_REGISTRY_DIR="${RESULTS_REGISTRY_DIR:-${REGISTRIES_DATA_DIR}/results_registry}"
export STATS_CACHE_DIR="${STATS_CACHE_DIR:-$RESULTS_REGISTRY_DIR}"

# Sample registry path (used by launchers for FK validation before writing).
export SAMPLE_REGISTRY="${SAMPLE_REGISTRY:-${REGISTRIES_DATA_DIR}/sample_registry}"

# Active sample group for results_registry writes. Default "all_226" — the
# full 226-sample cohort, registered by load_bridge.R STEP 6.5. Override
# for subset analyses (e.g. SAMPLE_GROUP=unrelated_81).
export SAMPLE_GROUP="${SAMPLE_GROUP:-all_226}"

# ── Population structure (MODULE_2B) ────────────────────────────────────────
POPSTRUCT_DIR="${BASE}/popstruct_thin"
export BEAGLE_DIR="${POPSTRUCT_DIR}/04_beagle_byRF_majmin"
SAMPLE_LIST="${SAMPLE_LIST:-${POPSTRUCT_DIR}/list_of_samples_one_per_line_same_bamfile_list.tsv}"
export PRUNED_LIST="${PRUNED_LIST:-${POPSTRUCT_DIR}/06_relatedness/pruned_samples.txt}"
export RELATEDNESS_RES="${POPSTRUCT_DIR}/06_relatedness/catfish_226_relatedness.res"

# NGSadmix results (ancestry config uses 05_ngsadmix_global/runs_thin500/)
NGSADMIX_DIR="${POPSTRUCT_DIR}/05_ngsadmix_global"
BEST_SEED_BY_K="${NGSADMIX_DIR}/runs_thin500/best_seed_by_K.tsv"
# Legacy-capital-N alias fallback
[[ -d "${POPSTRUCT_DIR}/05_NGSadmix" && ! -d "$NGSADMIX_DIR" ]] && \
  NGSADMIX_DIR="${POPSTRUCT_DIR}/05_NGSadmix"

# ── Unified ancestry (instant Q, dispatchers) ──────────────────────────────
ANCESTRY_DIR="${BASE}/unified_ancestry"
export ANCESTRY_CONFIG="${ANCESTRY_DIR}/00_ancestry_config.sh"
export INSTANT_Q_R="${ANCESTRY_DIR}/wrappers/instant_q.R"
export DISPATCHER_R="${ANCESTRY_DIR}/dispatchers/region_stats_dispatcher.R"

# Ancestry cache — moved OUT of the code tree (chat 15). The K sweep writes
# to <LOCAL_Q_DIR>/K<NN>/. Canonical K's summary gets flattened into precomp RDS.
export LOCAL_Q_DIR="${LOCAL_Q_DIR:-${BASE}/ancestry_cache}"
export CANONICAL_K="${CANONICAL_K:-${DEFAULT_K:-8}}"
export K_SWEEP="${K_SWEEP:-2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20}"

# ── MODULE_2B config ────────────────────────────────────────────────────────
MODULE2B_DIR="${BASE}/MODULE_2B_Structure"
export MODULE2B_CONFIG="${MODULE2B_DIR}/00_module2b_config.sh"

# ── Results subdirectories ───────────────────────────────────────────────────
DOSAGE_DIR="${INVDIR}/04_dosage_by_chr"
LPCA_DIR="${INVDIR}/05_local_pca"
MDS_DIR="${INVDIR}/06_mds_candidates"
BRIDGE_DIR="${INVDIR}/07_het_overlap"
REGIONAL_DIR="${INVDIR}/08_regional_pca"
COMBINED_DIR="${INVDIR}/09_combined_plots"
FOLLOWUP_RESULTS="${INVDIR}/20_candidate_followup"
PLOTS_RESULTS="${INVDIR}/21_candidate_plots"

# ── Shared discovery prefix ──────────────────────────────────────────────────
MDS_PREFIX_BASENAME="inversion_localpca"
MDS_PREFIX="${MDS_DIR}/${MDS_PREFIX_BASENAME}"

# ── Candidate tables ─────────────────────────────────────────────────────────
# CONTROL candidates (500 kb fixed-distance merge from STEP10_v2)
CAND_FILE="${MDS_PREFIX}.candidate_regions.tsv.gz"
# SNAKE candidates (primary method — from STEP10e_v3 merge_A regions)
SNAKE_CAND_FILE="${MDS_DIR}/snake_regions_multiscale/snake_candidate_regions.tsv.gz"

# ── R config (for R scripts that source a config.R) ─────────────────────────
R_CONFIG="${FOLLOWUPDIR}/config_inversion_followup.R"

# ── SV prior directory (STEP_C00 output) ────────────────────────────────────
# Shared between 2c_precomp (STEP_C00 writes, STEP_C01a reads) and any
# downstream consumer that wants sample × inv_id genotypes.
export SV_PRIOR_DIR="${SV_PRIOR_DIR:-${INVDIR}/03_sv_prior}"

# ── GHSL v6 directory (chat 14) ──────────────────────────────────────────────
# STEP_C04_snake3_ghsl_v6.R (heavy engine) + STEP_C04b_snake3_ghsl_classify.R
# emit three per-chrom RDS plus genome-wide TSVs here. Consumed by
# phase_4/4b/lib_ghsl_confirmation.R and phase_2/2d/run_all.R.
export GHSL_DIR="${GHSL_DIR:-${MDS_DIR}/snake_regions_multiscale/ghsl_v6}"
export GHSL_PANEL_LIB="${GHSL_PANEL_LIB:-${UTILS_DIR}/lib_ghsl_panel.R}"

# ── Defaults ─────────────────────────────────────────────────────────────────

# Theta / tP bridge defaults
THETA_WINDOW_BP_DEFAULT=50000
THETA_STEP_BP_DEFAULT=10000
THETA_SCALE_TAG_DEFAULT="win${THETA_WINDOW_BP_DEFAULT}.step${THETA_STEP_BP_DEFAULT}"

# Core local-PCA defaults
NPC=2
WINSIZE=100
WINSTEP=20

# Main MDS / discovery defaults
MDS_MODE_DEFAULT="chunked_2x"
MDS_SEED_DEFAULT=42
MDS_DIMS_DEFAULT=20
Z_THRESH_DEFAULT=3
BACKGROUND_MULTIPLIER_DEFAULT=2

# Positional control merge defaults
GAP_BP_CONTROL_DEFAULT=500000
MIN_WINDOWS_CONTROL_DEFAULT=3

# Snake output directories
SNAKE1_DIR="${MDS_DIR}/snake_regions_multiscale"
SNAKE1_DIAG_DIR="${MDS_DIR}/snake_diagnostics_v3"
SNAKE2_DIR="${MDS_DIR}/snake2_community"
SNAKE3_DIR="${MDS_DIR}/snake3_ghsl"
SNAKE4_DIR="${MDS_DIR}/snake4_threeband"
CONSENSUS_DIR="${MDS_DIR}/three_snake_consensus"
HMM_DIR="${MDS_DIR}/regime_hmm"
RESCUE_DIR="${MDS_DIR}/rescue_pass2"

# Precomp directories (v9.1)
PRECOMP_DIR="${SNAKE1_DIR}/precomp"
SIM_MATS_DIR="${PRECOMP_DIR}/sim_mats"

# Backward-compat aliases
SNAKE_PRIMARY_DIR="${SNAKE1_DIR}"
SNAKE_SECONDARY_DIR="${SNAKE1_DIR}"
SNAKE_DIAG_DIR="${SNAKE1_DIAG_DIR}"

# v8.5 output directories
LANDSCAPE_DIR="${MDS_DIR}/snake_regions_multiscale/landscape"
TRIANGLE_DIR="${MDS_DIR}/snake_regions_multiscale/triangles"
SCORING_DIR="${MDS_DIR}/snake_regions_multiscale/scoring"
DECOMPOSITION_DIR="${MDS_DIR}/snake_regions_multiscale/decomposition"
REGIME_DIR="${MDS_DIR}/snake_regions_multiscale/regime_engine"

# NN-smoothed similarity scales (v8.5.3; expanded 2026-04-17 FIX 21)
export NN_SIM_SCALES="20,40,80,120,160,200,240,320"
TRIANGLE_MULTISCALE_DIR="${MDS_DIR}/snake_regions_multiscale/triangles_multiscale"

# Current method philosophy
PRIMARY_CANDIDATE_METHOD="snake"
CONTROL_CANDIDATE_METHOD="posmerge_500kb"

# ---------------------------------------------------------------------------
# Original Li & Ralph 2019 / Huang et al. 2020 / Mérot et al. comparison
# ---------------------------------------------------------------------------
CLASSIC_LOCALPCA_OUTLIER_Z_DEFAULT=2
CLASSIC_MAX_GAP_WINDOWS_DEFAULT=20
CLASSIC_WINDOW_MODE_DEFAULT="fixed_bp"
CLASSIC_WINDOW_BP_SET_DEFAULT="1000,10000,100000"
CLASSIC_LD_R2_THRESHOLD_DEFAULT=0.4
CLASSIC_MIN_SNPS_FOR_GENOTYPING_DEFAULT=10000
CLASSIC_OUTLIER_SD_DESCRIPTION="MDS score > 2 SD above mean across windows"

# Regional PCA / genotype interpretation defaults
REGIONAL_PCA_DEFAULT_K_MIN=2
REGIONAL_PCA_DEFAULT_K_MAX=5
LEGACY_K3_REFERENCE_DEFAULT=3

# Reproducibility / comparison labels
METHOD_BRANCH_DEFAULT="snake_primary"
METHOD_BRANCH_CLASSIC="lostruct_classic"
METHOD_BRANCH_CONTROL="posmerge_control"

# ── Helper functions ─────────────────────────────────────────────────────────
inv_timestamp() { date '+%F %T'; }
inv_log() { echo "[$(inv_timestamp)] [INV] $*"; }
inv_err() { echo "[$(inv_timestamp)] [INV] [ERROR] $*" >&2; }
inv_die() { inv_err "$@"; exit 1; }

inv_check_file() {
  local f="$1" label="${2:-file}"
  [[ -s "$f" ]] || inv_die "Missing or empty ${label}: ${f}"
}

inv_init_dirs() {
  mkdir -p \
    "${INVDIR}" \
    "${DOSAGE_DIR}" "${LPCA_DIR}" "${MDS_DIR}" "${BRIDGE_DIR}" \
    "${REGIONAL_DIR}" "${COMBINED_DIR}" \
    "${FOLLOWUP_RESULTS}" "${PLOTS_RESULTS}" \
    "${LANDSCAPE_DIR}" "${TRIANGLE_DIR}" "${SCORING_DIR}" \
    "${DECOMPOSITION_DIR}" "${REGIME_DIR}" \
    "${PRECOMP_DIR}" "${SIM_MATS_DIR}" \
    "${TRIANGLE_MULTISCALE_DIR}" \
    "${REGISTRY_DIR}" \
    "${LOCAL_Q_DIR}" "${STATS_CACHE_DIR}" \
    "${SV_PRIOR_DIR}" "${GHSL_DIR}" \
    "${INVDIR}/logs"
}
