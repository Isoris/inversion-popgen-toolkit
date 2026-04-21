#!/usr/bin/env bash
# =============================================================================
# pipeline_bridge.sh — Universal bash loader for cross-module wiring
#
# Source this at the top of any SLURM launcher or shell script.
# It exports all env vars that registry_bridge.R and R scripts need.
#
# Usage (in SLURM launchers):
#   set -a; source pipeline_bridge.sh; set +a
#
# What it does:
#   1. Sources inversion config (if present)
#   2. Sources ancestry config (if present)
#   3. Sources module2b config (if present)
#   4. Exports REGISTRY_BRIDGE (+ legacy LOAD_BRIDGE alias), REGISTRY_DIR, sample paths
#   5. Ensures registry directory exists
#   6. Compiles C binaries if missing (instant_q, region_popstats)
#
# It will NOT conflict between configs because:
#   - All three configs use the same BASE
#   - Module-specific vars are prefixed (INVDIR, MODULE2B_*, etc.)
#   - Shared vars (REF, BAMLIST, SAMPLES_IND) resolve identically
# =============================================================================

# ── Project root ─────────────────────────────────────────────────────────────
export BASE="${BASE:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"

# ── Self-location ────────────────────────────────────────────────────────────
_PB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Source inversion config ──────────────────────────────────────────────────
# Prefer flattened layout; fall back to legacy v8.5.
for _try in \
    "${BASE}/00_inversion_config.sh" \
    "${BASE}/inversion_modules/00_inversion_config.sh" \
    "${BASE}/inversion_codebase_v8.5/00_inversion_config.sh"; do
  if [[ -f "$_try" ]]; then
    _INV_CFG="$_try"; break
  fi
done
if [[ -n "${_INV_CFG:-}" && -f "$_INV_CFG" ]]; then
  source "$_INV_CFG"
fi
unset _try

# ── Source ancestry config ───────────────────────────────────────────────────
_ANC_CFG="${BASE}/unified_ancestry/00_ancestry_config.sh"
if [[ -f "$_ANC_CFG" ]]; then
  source "$_ANC_CFG"
fi

# ── Source module2b config ───────────────────────────────────────────────────
_M2B_CFG="${BASE}/MODULE_2B_Structure/00_module2b_config.sh"
if [[ -f "$_M2B_CFG" ]]; then
  source "$_M2B_CFG"
fi

# ── Cross-module paths ───────────────────────────────────────────────────────
export SAMPLES_IND="${SAMPLES_IND:-${BASE}/het_roh/01_inputs_check/samples.ind}"
export SAMPLE_LIST="${SAMPLE_LIST:-${BASE}/popstruct_thin/list_of_samples_one_per_line_same_bamfile_list.tsv}"
export PRUNED_LIST="${PRUNED_LIST:-${BASE}/popstruct_thin/06_relatedness/pruned_samples.txt}"

# ── Sample registry ──────────────────────────────────────────────────────────
export REGISTRY_DIR="${REGISTRY_DIR:-${BASE}/sample_registry}"
mkdir -p "${REGISTRY_DIR}/groups" "${REGISTRY_DIR}/backups" 2>/dev/null || true

# ── Registry bridge (for R scripts to source) ────────────────────────────────
# Chat-18: the bridge file was renamed utils/load_bridge.R → utils/registry_bridge.R
# to make it obvious that this is the registry entry point. The old name is
# preserved as a shim so in-flight scripts keep working. Search order below
# prefers the new name but accepts the old one if the rename hasn't been
# deployed yet.
#
# Exports two env vars so scripts and SLURM jobs can find the bridge:
#   REGISTRY_BRIDGE (preferred, post-chat-18)
#   LOAD_BRIDGE     (legacy alias, same value, kept for compatibility)
if [[ -z "${REGISTRY_BRIDGE:-}" ]]; then
  for _try in \
    "${_PB_DIR}/registry_bridge.R" \
    "${BASE}/utils/registry_bridge.R" \
    "${BASE}/inversion_modules/utils/registry_bridge.R" \
    "${_PB_DIR}/../utils/registry_bridge.R" \
    "${_PB_DIR}/load_bridge.R" \
    "${BASE}/utils/load_bridge.R" \
    "${BASE}/inversion_modules/utils/load_bridge.R" \
    "${_PB_DIR}/../utils/load_bridge.R" \
    "${BASE}/inversion_codebase_v8.5/utils/load_bridge.R"; do
    if [[ -f "$_try" ]]; then
      export REGISTRY_BRIDGE="$_try"
      break
    fi
  done
  unset _try
  : "${REGISTRY_BRIDGE:=${BASE}/utils/registry_bridge.R}"
fi
# Legacy alias — keeps old SLURM launchers and R scripts working unchanged.
export LOAD_BRIDGE="${LOAD_BRIDGE:-${REGISTRY_BRIDGE}}"

# ── Unified ancestry module paths ────────────────────────────────────────────
export ANCESTRY_CONFIG="${_ANC_CFG}"
export INSTANT_Q_R="${BASE}/unified_ancestry/wrappers/instant_q.R"
export DISPATCHER_R="${BASE}/unified_ancestry/dispatchers/region_stats_dispatcher.R"
export INSTANT_Q_BIN="${BASE}/unified_ancestry/src/instant_q"
export POPSTATS_BIN="${BASE}/unified_ancestry/engines/fst_dxy/region_popstats"
export HOBS_WINDOWER_BIN="${BASE}/unified_ancestry/engines/hobs_hwe/scripts/hobs_windower"

# 2026-04-17: data dirs kept OUT of the code tree.
# LOCAL_Q_DIR defaults to $BASE/ancestry_cache (ancestry config also sets
# this; the `:=` default below respects whatever the ancestry config picked).
# RESULTS_REGISTRY_DIR is the chat-16 rewrite of chat-15 stats_cache — a
# first-class fourth registry. See registries/DATABASE_DESIGN.md.
: "${LOCAL_Q_DIR:=${BASE}/ancestry_cache}"
: "${RESULTS_REGISTRY_DIR:=${BASE}/registries/data/results_registry}"
: "${STATS_CACHE_DIR:=$RESULTS_REGISTRY_DIR}"   # deprecated alias, same dir
: "${SAMPLE_REGISTRY:=${BASE}/registries/data/sample_registry}"
: "${SAMPLE_GROUP:=all_226}"
: "${CANONICAL_K:=${DEFAULT_K:-8}}"
export LOCAL_Q_DIR RESULTS_REGISTRY_DIR STATS_CACHE_DIR \
       SAMPLE_REGISTRY SAMPLE_GROUP CANONICAL_K

# ── BEAGLE directory ─────────────────────────────────────────────────────────
export BEAGLE_DIR="${BASE}/popstruct_thin/04_beagle_byRF_majmin"

# ── Theme file ───────────────────────────────────────────────────────────────
export THEME_FILE="${BASE}/unified_ancestry/utils/theme_systems_plate.R"

# ── Rscript binary ──────────────────────────────────────────────────────────
export RSCRIPT_BIN="${RSCRIPT_BIN:-/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript}"

# ── Compile C binaries if missing ────────────────────────────────────────────
_compile_if_missing() {
  local bin="$1" src="$2" flags="${3:-}"
  if [[ ! -x "$bin" && -f "$src" ]]; then
    echo "[pipeline_bridge] Compiling $(basename "$bin")..."
    local dir="$(dirname "$src")"
    if [[ -f "${dir}/Makefile" ]]; then
      make -C "$dir" -j4 2>/dev/null || \
        gcc -O2 -o "$bin" "$src" $flags -lm -lz -lpthread 2>/dev/null || \
        echo "[pipeline_bridge] WARN: compile failed for $(basename "$bin")"
    else
      gcc -O2 -o "$bin" "$src" $flags -lm -lz -lpthread 2>/dev/null || \
        echo "[pipeline_bridge] WARN: compile failed for $(basename "$bin")"
    fi
  fi
}

_compile_if_missing "${INSTANT_Q_BIN}" "${BASE}/unified_ancestry/src/instant_q.cpp" "-lstdc++"
_compile_if_missing "${POPSTATS_BIN}" "${BASE}/unified_ancestry/engines/fst_dxy/region_popstats.c" ""
_compile_if_missing "${HOBS_WINDOWER_BIN}" "${BASE}/unified_ancestry/engines/hobs_hwe/scripts/hobs_windower.c" ""

# ── Helper: initialize sample registry from R ────────────────────────────────
init_registry_if_needed() {
  local master="${REGISTRY_DIR}/sample_master.tsv"
  if [[ ! -s "$master" ]]; then
    echo "[pipeline_bridge] Initializing sample registry..."
    "${RSCRIPT_BIN}" -e "
      source('${REGISTRY_BRIDGE}')
      # Registry auto-initializes from samples.ind on first load
      cat('[registry] Master:', nrow(reg\$get_master()), 'samples\n')

      # Register canonical groups
      all_ids <- get_sample_ids()
      register_group('all_226', all_ids, description = 'All 226 QC-pass samples')

      pruned <- get_pruned_ids()
      if (length(pruned) > 0) {
        register_group('unrelated_81', pruned,
                        description = 'NAToRA pruned unrelated set')
      }
    "
  fi
}

# ── Convenience: run init on first use ───────────────────────────────────────
init_registry_if_needed

# ── Summary ──────────────────────────────────────────────────────────────────
_pb_log() { echo "[$(date '+%F %T')] [pipeline_bridge] $*"; }
_pb_log "BASE=$BASE"
_pb_log "REGISTRY_BRIDGE=$REGISTRY_BRIDGE"
_pb_log "REGISTRY_DIR=$REGISTRY_DIR"
_pb_log "ANCESTRY_CONFIG=$ANCESTRY_CONFIG"
[[ -x "$INSTANT_Q_BIN" ]] && _pb_log "instant_q: compiled" || _pb_log "instant_q: NOT compiled"
[[ -x "$POPSTATS_BIN" ]]  && _pb_log "popstats: compiled"  || _pb_log "popstats: NOT compiled"

# Cleanup
unset _PB_DIR _INV_CFG _ANC_CFG _M2B_CFG _pb_log _compile_if_missing
