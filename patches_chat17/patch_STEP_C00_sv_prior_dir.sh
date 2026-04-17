#!/usr/bin/env bash
# =============================================================================
# patch_STEP_C00_sv_prior_dir.sh
#
# Chat-17 patch. The `SV_PRIOR_DIR` env var in 00_inversion_config.sh points
# to `${INVDIR}/03_sv_prior`, but STEP_C00_build_sv_prior.R hardcodes its
# output dir to `${INVDIR}/06_mds_candidates/snake_regions_multiscale/sv_prior`.
#
# This patch makes STEP_C00 respect the env var when set, falling back to
# the old hardcoded path for back-compat. Canonical resolution:
#
#   OLD (line 120):
#     SV_PRIOR_DIR <- file.path(INVDIR, "06_mds_candidates/snake_regions_multiscale/sv_prior")
#
#   NEW:
#     SV_PRIOR_DIR <- Sys.getenv("SV_PRIOR_DIR",
#       file.path(INVDIR, "06_mds_candidates/snake_regions_multiscale/sv_prior"))
#
# After applying, either:
#   (a) Export SV_PRIOR_DIR in 00_inversion_config.sh to whatever path is
#       already used on disk (the MDS_DIR one), or
#   (b) Change 00_inversion_config.sh so SV_PRIOR_DIR points at the
#       MDS_DIR location. The chat-17 C00 launcher already uses the script's
#       hardcoded path, so either direction works — the patch just makes the
#       script honour whatever the env says.
#
# Backup suffix: .bk_chat17
#
# Usage:  bash patch_STEP_C00_sv_prior_dir.sh [--dry-run]
# =============================================================================

set -euo pipefail

BASE="${BASE:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"
SCRIPT="${BASE}/inversion_modules/phase_2_discovery/2c_precomp/STEP_C00_build_sv_prior.R"
DRY_RUN=0
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=1

[[ -f "${SCRIPT}" ]] || { echo "[patch] missing: ${SCRIPT}" >&2; exit 1; }

# Use Python for literal string replace — no regex pitfalls with embedded parens / quotes.
python3 - "${SCRIPT}" "${DRY_RUN}" <<'PYEOF'
import sys, shutil
script_path, dry_run = sys.argv[1], sys.argv[2] == "1"

OLD = 'SV_PRIOR_DIR <- file.path(INVDIR, "06_mds_candidates/snake_regions_multiscale/sv_prior")'
NEW = ('SV_PRIOR_DIR <- Sys.getenv("SV_PRIOR_DIR",\n'
       '  file.path(INVDIR, "06_mds_candidates/snake_regions_multiscale/sv_prior"))')

with open(script_path, "r") as fh:
    text = fh.read()

if 'Sys.getenv("SV_PRIOR_DIR"' in text:
    print(f"[patch] already applied — SV_PRIOR_DIR env lookup present in {script_path}")
    sys.exit(0)

if OLD not in text:
    print(f"[patch] target line not found in {script_path}", file=sys.stderr)
    print(f"        expected: {OLD}", file=sys.stderr)
    sys.exit(2)

count = text.count(OLD)
if count != 1:
    print(f"[patch] expected exactly 1 occurrence, found {count}", file=sys.stderr)
    sys.exit(2)

patched = text.replace(OLD, NEW)

print(f"[patch] {script_path}")
if dry_run:
    print("  (dry-run — would replace the hardcoded SV_PRIOR_DIR line)")
    print(f"  OLD: {OLD}")
    print(f"  NEW: {NEW.splitlines()[0]}")
    print(f"       {NEW.splitlines()[1]}")
    sys.exit(0)

shutil.copyfile(script_path, script_path + ".bk_chat17")
with open(script_path, "w") as fh:
    fh.write(patched)
print(f"  backup: {script_path}.bk_chat17")
PYEOF

if [[ ${DRY_RUN} -eq 0 ]]; then
  echo ""
  echo "[patch] verify:"
  echo "         grep -n 'SV_PRIOR_DIR <-' ${SCRIPT}"
  echo "         diff ${SCRIPT}.bk_chat17 ${SCRIPT}"
fi
