#!/usr/bin/env bash
# =============================================================================
# registry_loader.sh — minimal bash helpers for pipeline launchers
# =============================================================================
# Provides shell access to the registry. This is NOT a full API mirror —
# R and Python do the heavy lifting. Bash scripts use this to:
#
#   - resolve registry paths from env
#   - check "does candidate X have validation level >= Y?" as a gate
#   - list candidates at a given tier for per-candidate SLURM arrays
#
# Usage:
#   source registries/api/bash/registry_loader.sh
#   registry_resolve_paths                        # sets REGISTRIES, EVIDENCE_REGISTRY, etc.
#   registry_list_candidates_by_tier 2            # echoes one cid per line
#   registry_check_validation LG12_17 SUPPORTED   # exit 0 if >=, 1 if <
# =============================================================================

# ── Resolve paths (idempotent) ───────────────────────────────────────────────
registry_resolve_paths() {
  : "${BASE:=/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"
  : "${REGISTRIES:=${BASE}/inversion-popgen-toolkit/registries}"
  : "${SAMPLE_REGISTRY:=${REGISTRIES}/data/sample_registry}"
  : "${EVIDENCE_REGISTRY:=${REGISTRIES}/data/evidence_registry}"
  : "${INTERVAL_REGISTRY:=${REGISTRIES}/data/interval_registry}"
  : "${SCHEMA_DIR:=${REGISTRIES}/schemas}"

  export BASE REGISTRIES SAMPLE_REGISTRY EVIDENCE_REGISTRY INTERVAL_REGISTRY SCHEMA_DIR
}

# ── List candidates at a given tier ──────────────────────────────────────────
# Reads every per_candidate/<cid>/keys.tsv and filters on q7_tier.
# Outputs one candidate_id per line.
registry_list_candidates_by_tier() {
  local max_tier="${1:-2}"
  registry_resolve_paths
  local percand="${EVIDENCE_REGISTRY}/per_candidate"
  [[ -d "${percand}" ]] || { echo "[registry] no per_candidate dir: ${percand}" >&2; return 1; }

  for d in "${percand}"/*/; do
    [[ -d "$d" ]] || continue
    local cid
    cid=$(basename "$d")
    local kf="$d/keys.tsv"
    [[ -f "$kf" ]] || continue
    # Grab current tier from keys.tsv (last row for q7_tier wins)
    local tier
    tier=$(awk -F'\t' '$1 == "q7_tier" { v = $2 } END { print v }' "$kf")
    [[ -z "$tier" ]] && continue
    if (( tier <= max_tier )); then
      echo "$cid"
    fi
  done
}

# ── Check validation level gate ──────────────────────────────────────────────
# registry_check_validation CID MIN_LEVEL
#   MIN_LEVEL ∈ {NONE, UNCERTAIN, SUPPORTED, VALIDATED}
# Returns 0 if current level >= MIN_LEVEL, 1 otherwise (SUSPECT always fails).
registry_check_validation() {
  local cid="$1"
  local min="$2"
  registry_resolve_paths

  declare -A rank=( [NONE]=0 [UNCERTAIN]=1 [SUPPORTED]=2 [VALIDATED]=3 [SUSPECT]=-1 )
  local min_rank="${rank[$min]:-0}"

  local kf="${EVIDENCE_REGISTRY}/per_candidate/${cid}/keys.tsv"
  [[ -f "$kf" ]] || return 1
  local current
  current=$(awk -F'\t' '$1 == "q6_group_validation" { v = $2 } END { print v }' "$kf")
  [[ -z "$current" ]] && current="NONE"
  local cur_rank="${rank[$current]:-0}"

  if (( cur_rank >= min_rank )); then
    return 0
  else
    return 1
  fi
}

# ── Convenience: echo validation level of a candidate ─────────────────────────
registry_get_validation() {
  local cid="$1"
  registry_resolve_paths
  local kf="${EVIDENCE_REGISTRY}/per_candidate/${cid}/keys.tsv"
  [[ -f "$kf" ]] || { echo "NONE"; return; }
  awk -F'\t' '$1 == "q6_group_validation" { v = $2 } END { if (v == "") v = "NONE"; print v }' "$kf"
}

# ── Convenience: check that a sample group exists ─────────────────────────────
registry_has_group() {
  local gid="$1"
  registry_resolve_paths
  [[ -f "${SAMPLE_REGISTRY}/groups/${gid}.txt" ]]
}

# If sourced, do nothing further. If executed directly, run self-test.
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  registry_resolve_paths
  echo "REGISTRIES=${REGISTRIES}"
  echo "SAMPLE_REGISTRY=${SAMPLE_REGISTRY}"
  echo "EVIDENCE_REGISTRY=${EVIDENCE_REGISTRY}"
  echo "INTERVAL_REGISTRY=${INTERVAL_REGISTRY}"
  echo "SCHEMA_DIR=${SCHEMA_DIR}"
  echo ""
  echo "# usage examples:"
  echo "#   registry_list_candidates_by_tier 2"
  echo "#   registry_check_validation LG12_17 SUPPORTED && echo ok || echo not yet"
  echo "#   registry_get_validation LG12_17"
  echo "#   registry_has_group inv_LG12_17_RECOMBINANT && echo yes"
fi
