#!/usr/bin/env bash
# =============================================================================
# LAUNCH_group_cheats.sh — phase 4d launcher
# =============================================================================
# Runs cheat30 (GDS), cheat6 (ancestry jackknife standalone), and burden
# regression, BUT gates each on the corresponding group validation level.
# This is the concrete demonstration of "downstream cheats check validation
# before running" that the architecture doc describes.
#
# Called by run_phase4.sh after 4c_hypothesis completes.
#
# Usage:
#   sbatch LAUNCH_group_cheats.sh --chroms LG01,LG02,...
# =============================================================================

#SBATCH -A lt200308
#SBATCH -J 4d_cheats
#SBATCH -t 04:00:00
#SBATCH -c 4
#SBATCH --mem=16G

set -euo pipefail

# ── Parse ────────────────────────────────────────────────────────────────────
CHROMS_CSV=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --chroms) CHROMS_CSV="$2"; shift 2 ;;
    *) shift ;;
  esac
done

# ── Config + registry paths ──────────────────────────────────────────────────
BASE="${BASE:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"
MODULE_DIR="${BASE}/inversion-popgen-toolkit"
# shellcheck disable=SC1091
source "${MODULE_DIR}/00_inversion_config.sh" 2>/dev/null || true
# shellcheck disable=SC1091
source "${MODULE_DIR}/registries/api/bash/registry_loader.sh"
registry_resolve_paths

echo "[LAUNCH_group_cheats] REGISTRIES=${REGISTRIES}"
echo "[LAUNCH_group_cheats] chroms=${CHROMS_CSV:-ALL}"

# ── Get candidate list (tier ≤ 3) ────────────────────────────────────────────
CIDS=($(registry_list_candidates_by_tier 3))
echo "[LAUNCH_group_cheats] ${#CIDS[@]} candidates at tier<=3"

# ── Per-candidate dispatch with validation gates ─────────────────────────────
for cid in "${CIDS[@]}"; do
  level=$(registry_get_validation "$cid")
  echo "[LAUNCH_group_cheats] ${cid}: validation=${level}"

  # cheat30 GDS — requires NONE (always runs)
  echo "  cheat30 GDS → running"
  Rscript "${MODULE_DIR}/cheats/cheat30_gds_by_genotype.R" \
    --candidate "$cid" \
    --outdir "${EVIDENCE_REGISTRY}/per_candidate/${cid}/raw/" \
    || echo "    cheat30 failed for $cid" >&2

  # cheat6 ancestry jackknife — requires UNCERTAIN
  if registry_check_validation "$cid" UNCERTAIN; then
    echo "  cheat6 jackknife → running"
    Rscript "${MODULE_DIR}/cheats/cheat6_ancestry_jackknife_v934.R" \
      --candidate "$cid" \
      --outdir "${EVIDENCE_REGISTRY}/per_candidate/${cid}/raw/" \
      || echo "    cheat6 failed for $cid" >&2
  else
    echo "  cheat6 jackknife → SKIP (validation below UNCERTAIN)"
  fi

  # Fst sub-block of age_evidence — requires SUPPORTED
  if registry_check_validation "$cid" SUPPORTED; then
    echo "  age_evidence Fst sub-block → running"
    Rscript "${MODULE_DIR}/scripts/compute_age_fst_subblock.R" \
      --candidate "$cid" \
      --outdir "${EVIDENCE_REGISTRY}/per_candidate/${cid}/raw/" \
      || echo "    age_fst failed for $cid" >&2
  else
    echo "  age_evidence Fst sub-block → SKIP (below SUPPORTED)"
  fi

  # Burden regression — requires VALIDATED
  if registry_check_validation "$cid" VALIDATED; then
    echo "  burden regression → running"
    Rscript "${MODULE_DIR}/burden/STEP_C01f_c_burden_regression.R" \
      --candidate "$cid" \
      --outdir "${EVIDENCE_REGISTRY}/per_candidate/${cid}/raw/" \
      || echo "    burden failed for $cid" >&2
  else
    echo "  burden regression → SKIP (below VALIDATED)"
  fi
done

echo "[LAUNCH_group_cheats] done"
