#!/usr/bin/env bash
# =============================================================================
# submit_precomp_chain.sh
#
# One-shot submitter for the full precomp chain, with SLURM dependencies so
# each step only starts if the previous one completed OK:
#
#   1. LAUNCH_instant_q_precompute.slurm   (ancestry K=8 — 28 array tasks)
#   2. LAUNCH_STEP_C00_sv_prior.slurm      (SV prior — single job, all chroms)
#   3. LAUNCH_STEP_C01a_precompute.slurm   (precompute with localQ + SV merge)
#
# Usage:
#   ./submit_precomp_chain.sh                 # default: K_SWEEP=8, all chroms
#   K_SWEEP="8 12" ./submit_precomp_chain.sh  # K=8 and K=12 precompute
#   SKIP_ANCESTRY=1 ./submit_precomp_chain.sh # skip step 1 (already cached)
#   SKIP_SV=1 ./submit_precomp_chain.sh       # skip step 2 (SV prior already built)
#
# Prereqs (pre-flight checks from chat-17 handoff):
#   - SAMPLE_GROUP=all_226 registered in sample_registry
#   - registry_loader.R sources clean on LANTA R
#   - load_bridge.R STEP 6.5 auto-registers all_226 on first init
# =============================================================================

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="${BASE:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"

ANC_LAUNCHER="${BASE}/unified_ancestry/launchers/LAUNCH_instant_q_precompute.slurm"
C00_LAUNCHER="${HERE}/LAUNCH_STEP_C00_sv_prior.slurm"
C01A_LAUNCHER="${HERE}/LAUNCH_STEP_C01a_precompute.slurm"

for f in "${C00_LAUNCHER}" "${C01A_LAUNCHER}"; do
  [[ -f "$f" ]] || { echo "Missing: $f" >&2; exit 1; }
done

K_SWEEP="${K_SWEEP:-8}"
SKIP_ANCESTRY="${SKIP_ANCESTRY:-0}"
SKIP_SV="${SKIP_SV:-0}"

DEPS=""

# ── 1. Ancestry precompute ───────────────────────────────────────────────
if [[ "${SKIP_ANCESTRY}" == "1" ]]; then
  echo "[chain] SKIP_ANCESTRY=1 — skipping step 1"
else
  [[ -f "${ANC_LAUNCHER}" ]] || { echo "Missing ancestry launcher: ${ANC_LAUNCHER}" >&2; exit 1; }
  echo "[chain] Submitting ancestry precompute (K_SWEEP=${K_SWEEP})"
  JID_ANC=$(K_SWEEP="${K_SWEEP}" sbatch --parsable "${ANC_LAUNCHER}")
  echo "[chain]   ancestry job: ${JID_ANC}"
  DEPS="afterok:${JID_ANC}"
fi

# ── 2. SV prior ──────────────────────────────────────────────────────────
if [[ "${SKIP_SV}" == "1" ]]; then
  echo "[chain] SKIP_SV=1 — skipping step 2"
else
  echo "[chain] Submitting SV prior${DEPS:+ (depends: ${DEPS})}"
  # SV prior doesn't actually need ancestry output, but we serialize it here
  # so the whole chain either proceeds or blocks together. If you prefer
  # parallel, change to `--dependency=singleton` or drop the dep here.
  if [[ -n "${DEPS}" ]]; then
    JID_SV=$(sbatch --parsable --dependency="${DEPS}" "${C00_LAUNCHER}")
  else
    JID_SV=$(sbatch --parsable "${C00_LAUNCHER}")
  fi
  echo "[chain]   SV prior job: ${JID_SV}"
  DEPS="afterok:${JID_SV}"
fi

# ── 3. C01a precompute — depends on BOTH ancestry and SV prior ───────────
# If either was skipped, its job id is not in DEPS; the chain relies on
# whoever ran last. Common case: ancestry + SV prior both submitted, DEPS
# is afterok:${JID_SV} — and JID_SV itself already depended on JID_ANC,
# so transitivity covers it.
echo "[chain] Submitting C01a precompute${DEPS:+ (depends: ${DEPS})}"
if [[ -n "${DEPS}" ]]; then
  JID_C01A=$(sbatch --parsable --dependency="${DEPS}" "${C01A_LAUNCHER}")
else
  JID_C01A=$(sbatch --parsable "${C01A_LAUNCHER}")
fi
echo "[chain]   C01a job: ${JID_C01A}"

echo ""
echo "[chain] All jobs submitted. Monitor with:"
echo "  squeue -u \$USER"
echo "  tail -f logs/inv_C01a_${JID_C01A}.out"
