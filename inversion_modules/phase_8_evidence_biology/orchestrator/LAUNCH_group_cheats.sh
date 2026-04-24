#!/usr/bin/env bash
# =============================================================================
# LAUNCH_group_cheats.sh — phase 4d production launcher
# =============================================================================
# Runs cheat30 (GDS), cheat6 (ancestry jackknife), cheat27/28/29 (mechanism),
# and conditionally the Fst age sub-block + burden regression — gating each
# on the corresponding group validation level.
#
# Called by run_phase4.sh after 4c_hypothesis completes.
#
# Usage:
#   sbatch LAUNCH_group_cheats.sh --chroms LG01,LG02,...
#
# FIX 50 (chat 10, 2026-04-17): promoted from LAUNCH_group_cheats_example.sh.
# Pass 15 (2026-04-24): scripts regrouped by question under phase_8_evidence_biology/.
# Path corrections:
#   ORIGINAL EXAMPLE                  PRE-PASS-15                           POST-PASS-15 (current)
#   --------------------------------  ------------------------------------  -------------------------------------
#   cheats/cheat30_gds_by_genotype.R  4f_group_dependent/cheat30_*.R        q5_age_and_origin/cheat30_*.R
#   cheats/cheat6_ancestry_*.R        4f_group_dependent/cheat6_*.R         q7_existence_audit/cheat6_*.R
#   (cheat27 / 28 / 29 — added in FIX 50)                                    q4_mechanism/cheat{27,28,29}_*.R
#   scripts/compute_age_fst_subblock  (ASPIRATIONAL — not in tarball)       still aspirational
#   burden/STEP_C01f_c_burden_*       (ASPIRATIONAL — not in tarball)       still aspirational; see SPEC_VS_REALITY.md
#
# Also added: cheat27 + cheat28 + cheat29 calls (these 3 exist in q4_mechanism/
# but the example launcher never called them). They require NONE validation
# level (mechanism analysis is group-independent), so no gate is needed.
#
# NOTE: The two aspirational scripts (compute_age_fst_subblock.R,
# STEP_C01f_c_burden_regression.R) are left in place as commented-out calls
# so that when the scripts arrive in the tarball, uncommenting is the only
# action needed. Their validation gates (SUPPORTED / VALIDATED) remain
# active if-guarded via registry_check_validation.
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
PHASE8_DIR="${MODULE_DIR}/inversion_modules/phase_8_evidence_biology"

# shellcheck disable=SC1091
source "${MODULE_DIR}/00_inversion_config.sh" 2>/dev/null || true
# shellcheck disable=SC1091
source "${MODULE_DIR}/registries/api/bash/registry_loader.sh"
registry_resolve_paths

echo "[LAUNCH_group_cheats] REGISTRIES=${REGISTRIES}"
echo "[LAUNCH_group_cheats] PHASE8_DIR=${PHASE8_DIR}"
echo "[LAUNCH_group_cheats] chroms=${CHROMS_CSV:-ALL}"

# ── Get candidate list (tier ≤ 3) ────────────────────────────────────────────
CIDS=($(registry_list_candidates_by_tier 3))
echo "[LAUNCH_group_cheats] ${#CIDS[@]} candidates at tier<=3"

# ── Per-candidate dispatch with validation gates ─────────────────────────────
for cid in "${CIDS[@]}"; do
  level=$(registry_get_validation "$cid")
  echo "[LAUNCH_group_cheats] ${cid}: validation=${level}"
  RAW_OUT="${EVIDENCE_REGISTRY}/per_candidate/${cid}/raw/"
  mkdir -p "${RAW_OUT}"

  # ── cheat30 GDS — requires NONE (always runs) ────────────────────────────
  echo "  cheat30 GDS → running"
  Rscript "${PHASE8_DIR}/q5_age_and_origin/cheat30_gds_by_genotype.R" \
    --candidate "$cid" \
    --outdir "${RAW_OUT}" \
    || echo "    cheat30 failed for $cid" >&2

  # ── cheat27 SD/NAHR substrate — requires NONE (junction sequence) ────────
  echo "  cheat27 SD/NAHR → running"
  Rscript "${PHASE8_DIR}/q4_mechanism/cheat27_sd_nahr_substrate.R" \
    --candidate "$cid" \
    --outdir "${RAW_OUT}" \
    || echo "    cheat27 failed for $cid" >&2

  # ── cheat28 tandem-repeat context — requires NONE ────────────────────────
  echo "  cheat28 tandem-repeat → running"
  Rscript "${PHASE8_DIR}/q4_mechanism/cheat28_tandem_repeat_context.R" \
    --candidate "$cid" \
    --outdir "${RAW_OUT}" \
    || echo "    cheat28 failed for $cid" >&2

  # ── cheat29 junction forensics — requires NONE ───────────────────────────
  echo "  cheat29 junction forensics → running"
  Rscript "${PHASE8_DIR}/q4_mechanism/cheat29_junction_forensics.R" \
    --candidate "$cid" \
    --outdir "${RAW_OUT}" \
    || echo "    cheat29 failed for $cid" >&2

  # ── cheat6 ancestry jackknife — requires UNCERTAIN ───────────────────────
  if registry_check_validation "$cid" UNCERTAIN; then
    echo "  cheat6 jackknife → running"
    Rscript "${PHASE8_DIR}/q7_existence_audit/cheat6_ancestry_jackknife_v934.R" \
      --candidate "$cid" \
      --outdir "${RAW_OUT}" \
      || echo "    cheat6 failed for $cid" >&2
  else
    echo "  cheat6 jackknife → SKIP (validation below UNCERTAIN)"
  fi

  # ── Fst sub-block of age_evidence — requires SUPPORTED ───────────────────
  # ASPIRATIONAL: compute_age_fst_subblock.R not yet in tarball (chat-9
  # Finding W). When it lands, uncomment the call below.
  if registry_check_validation "$cid" SUPPORTED; then
    echo "  age_evidence Fst sub-block → PENDING SCRIPT (SUPPORTED reached)"
    # Rscript "${MODULE_DIR}/scripts/compute_age_fst_subblock.R" \
    #   --candidate "$cid" \
    #   --outdir "${RAW_OUT}" \
    #   || echo "    age_fst failed for $cid" >&2
  else
    echo "  age_evidence Fst sub-block → SKIP (below SUPPORTED)"
  fi

  # ── Burden regression — requires VALIDATED ───────────────────────────────
  # ASPIRATIONAL: STEP_C01f_c_burden_regression.R not yet in tarball.
  if registry_check_validation "$cid" VALIDATED; then
    echo "  burden regression → PENDING SCRIPT (VALIDATED reached)"
    # Rscript "${MODULE_DIR}/burden/STEP_C01f_c_burden_regression.R" \
    #   --candidate "$cid" \
    #   --outdir "${RAW_OUT}" \
    #   || echo "    burden failed for $cid" >&2
  else
    echo "  burden regression → SKIP (below VALIDATED)"
  fi
done

echo "[LAUNCH_group_cheats] done"
