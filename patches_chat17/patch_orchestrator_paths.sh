#!/usr/bin/env bash
# =============================================================================
# patch_orchestrator_paths.sh
#
# Chat-17 patch. The phase-4 orchestrators (run_phase4.sh, run_phase4b.sh)
# carry stale paths from an earlier directory layout. This patch rewrites
# them in-place to match the current tree:
#
#   OLD                                          NEW
#   ${BASE}/inversion-popgen-toolkit         →   ${BASE}/inversion_modules
#   ${MODULE_DIR}/launchers/LAUNCH_X.sh      →   ${MODULE_DIR}/phase_4_postprocessing/4?/LAUNCH_X.sh
#   ${MODULE_DIR}/phase4b_rewrite/R          →   ${MODULE_DIR}/phase_4_postprocessing/4b_group_proposal
#   ${MODULE_DIR}/phase4b_rewrite/python     →   ${MODULE_DIR}/phase_4_postprocessing/4b_group_proposal
#
# Backup suffix: .bk_chat17
#
# Usage:  bash patch_orchestrator_paths.sh [--dry-run]
# =============================================================================

set -euo pipefail

BASE="${BASE:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"
MODULE_DIR="${MODULE_DIR:-${BASE}/inversion_modules}"
ORCH_DIR="${MODULE_DIR}/phase_4_postprocessing/orchestrator"
DRY_RUN=0
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=1

RUN_PHASE4="${ORCH_DIR}/run_phase4.sh"
RUN_PHASE4B="${ORCH_DIR}/run_phase4b.sh"

for f in "${RUN_PHASE4}" "${RUN_PHASE4B}"; do
  [[ -f "$f" ]] || { echo "[patch] missing: $f" >&2; exit 1; }
done

patch_file() {
  local file="$1"; shift
  echo "[patch] ${file}"
  if [[ ${DRY_RUN} -eq 1 ]]; then
    echo "  (dry-run — showing what would change)"
    diff <(cat "$file") <(eval "$*" < "$file") | head -40 || true
    return
  fi
  cp "$file" "${file}.bk_chat17"
  eval "$*" < "${file}.bk_chat17" > "$file"
  echo "  backup: ${file}.bk_chat17"
}

# ── run_phase4.sh ────────────────────────────────────────────────────────
# 1. Base module path
# 2. Launcher paths: ${MODULE_DIR}/launchers/LAUNCH_X.sh →
#    per-phase subdirs. The 4a/4c/4e launchers live under
#    phase_4_postprocessing/4?/LAUNCH_X.sh.
patch_file "${RUN_PHASE4}" \
  "sed -E \
    -e 's|\\\$\\{BASE\\}/inversion-popgen-toolkit|\\\$\\{BASE\\}/inversion_modules|g' \
    -e 's|\\\$\\{MODULE_DIR\\}/launchers/LAUNCH_C01d_scoring_pass1\\.sh|\\\$\\{MODULE_DIR\\}/phase_4_postprocessing/4a_existence_layers/LAUNCH_C01d_scoring_pass1.sh|g' \
    -e 's|\\\$\\{MODULE_DIR\\}/launchers/LAUNCH_C01g_boundary\\.sh|\\\$\\{MODULE_DIR\\}/phase_4_postprocessing/4a_existence_layers/LAUNCH_C01g_boundary.sh|g' \
    -e 's|\\\$\\{MODULE_DIR\\}/launchers/LAUNCH_C01f_hypothesis\\.sh|\\\$\\{MODULE_DIR\\}/phase_4_postprocessing/4c_group_validation/LAUNCH_C01f_hypothesis.sh|g' \
    -e 's|\\\$\\{MODULE_DIR\\}/launchers/LAUNCH_C01d_scoring_pass2\\.sh|\\\$\\{MODULE_DIR\\}/phase_4_postprocessing/4e_final_classification/LAUNCH_C01d_scoring_pass2.sh|g' \
    -e 's|\\\$\\{MODULE_DIR\\}/launchers/LAUNCH_characterize_classify\\.sh|\\\$\\{MODULE_DIR\\}/phase_4_postprocessing/4e_final_classification/LAUNCH_characterize_classify.sh|g' \
    -e 's|\\\$\\{MODULE_DIR\\}/launchers/LAUNCH_group_cheats\\.sh|\\\$\\{MODULE_DIR\\}/phase_4_postprocessing/orchestrator/LAUNCH_group_cheats.sh|g'"

# ── run_phase4b.sh ───────────────────────────────────────────────────────
# MODULE_DIR base + R / python script paths.
patch_file "${RUN_PHASE4B}" \
  "sed -E \
    -e 's|\\\$\\{BASE\\}/inversion-popgen-toolkit|\\\$\\{BASE\\}/inversion_modules|g' \
    -e 's|\\\$\\{MODULE_DIR\\}/phase4b_rewrite/R|\\\$\\{MODULE_DIR\\}/phase_4_postprocessing/4b_group_proposal|g' \
    -e 's|\\\$\\{MODULE_DIR\\}/phase4b_rewrite/python|\\\$\\{MODULE_DIR\\}/phase_4_postprocessing/4b_group_proposal|g'"

echo ""
echo "[patch] done."
if [[ ${DRY_RUN} -eq 0 ]]; then
  echo "[patch] verify with:"
  echo "         bash -n ${RUN_PHASE4}"
  echo "         bash -n ${RUN_PHASE4B}"
  echo "         diff ${RUN_PHASE4}.bk_chat17 ${RUN_PHASE4}"
fi
