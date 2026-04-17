#!/usr/bin/env bash
# =============================================================================
# patch_rename_window_dt.sh
#
# Chat-17 patch. Renames the misleadingly-named per-window diagnostic table
# variable and its on-disk file:
#
#   inv_like_dt                      →  window_dt
#   window_inv_likeness.tsv.gz       →  window_dt.tsv.gz
#   snake_inv_likeness.tsv.gz        →  window_dt.tsv.gz  (stale ref in diag_common.R)
#
# Rationale: the table and file were named after their flagship column
# (`inv_likeness`, a het + dip + discreteness composite) as if that were
# their only content. In reality the same table carries ~40 columns
# including MDS z-scores, dip tests, dosage stats, band stamps,
# family_likeness, localQ_*, SV overlap — the name systematically
# misleads anyone auditing the code. "window_dt" is honest: it's the
# per-window diagnostic table. Column `inv_likeness` is kept (that name
# IS accurate for the column itself).
#
# Files touched (verified against the chat-16 tarball; run against the
# HPC tree to confirm nothing new has appeared):
#
#   Variable rename (inv_like_dt → window_dt):
#     inversion_modules/phase_2_discovery/2c_precomp/STEP_C01a_precompute.R
#     inversion_modules/phase_2_discovery/2c_precomp/STEP_C01b_1_seeded_regions.R
#     inversion_modules/phase_2_discovery/2c_precomp/diags/STEP_C01a_diag_common.R
#     inversion_modules/phase_2_discovery/2c_precomp/diags/STEP_C01a_diag_matrices.R
#     inversion_modules/phase_2_discovery/2c_precomp/diags/STEP_C01a_diag_mds.R
#     inversion_modules/phase_2_discovery/2c_precomp/diags/STEP_C01a_diag_profiles.R
#     inversion_modules/phase_2_discovery/2c_precomp/diags/STEP_C01a_diag_summaries.R
#     inversion_modules/phase_2_discovery/2c_precomp/patches/patch_C01a_precompute_flashlight.R
#     inversion_modules/phase_4_postprocessing/4d_group_dependent/cheat28_tandem_repeat_context.R
#     unified_ancestry/wrappers/instant_q.R
#
#   Filename rename (window_inv_likeness.tsv.gz → window_dt.tsv.gz):
#     inversion_modules/phase_2_discovery/2c_precomp/STEP_C01a_precompute.R       (writer)
#     inversion_modules/phase_2_discovery/2c_precomp/STEP_C01b_1_seeded_regions.R (reader)
#     inversion_modules/phase_2_discovery/2c_precomp/diags/STEP_C01a_diag_common.R (stale ref)
#
# Unchanged:
#   - Column name `inv_likeness` (accurate for the column itself)
#   - Function `compute_inv_likeness_all()` (still returns the table that
#     bears the inv_likeness column, among ~40 others)
#   - The precomp RDS is not touched — its internal list element is `dt`,
#     not inv_like_dt, so it's already OK.
#
# Backup suffix: .bk_chat17
# Usage:  bash patch_rename_window_dt.sh [--dry-run]
# =============================================================================

set -euo pipefail

BASE="${BASE:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"
DRY_RUN=0
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=1

# List of files (relative to $BASE)
FILES=(
  "inversion_modules/phase_2_discovery/2c_precomp/STEP_C01a_precompute.R"
  "inversion_modules/phase_2_discovery/2c_precomp/STEP_C01b_1_seeded_regions.R"
  "inversion_modules/phase_2_discovery/2c_precomp/diags/STEP_C01a_diag_common.R"
  "inversion_modules/phase_2_discovery/2c_precomp/diags/STEP_C01a_diag_matrices.R"
  "inversion_modules/phase_2_discovery/2c_precomp/diags/STEP_C01a_diag_mds.R"
  "inversion_modules/phase_2_discovery/2c_precomp/diags/STEP_C01a_diag_profiles.R"
  "inversion_modules/phase_2_discovery/2c_precomp/diags/STEP_C01a_diag_summaries.R"
  "inversion_modules/phase_2_discovery/2c_precomp/patches/patch_C01a_precompute_flashlight.R"
  "inversion_modules/phase_4_postprocessing/4d_group_dependent/cheat28_tandem_repeat_context.R"
  "unified_ancestry/wrappers/instant_q.R"
)

N_PATCHED=0
N_SKIPPED=0
N_MISSING=0

for rel in "${FILES[@]}"; do
  path="${BASE}/${rel}"
  if [[ ! -f "${path}" ]]; then
    echo "[patch] MISSING: ${rel}" >&2
    N_MISSING=$((N_MISSING + 1))
    continue
  fi

  # Python does the work — literal replace, no regex escaping, handles all
  # three renames + idempotency check atomically per file.
  # Capture python exit code safely under `set -e`. The `|| rc=$?` pattern
  # lets us read the non-zero exit without the whole script dying.
  rc=0
  python3 - "${path}" "${DRY_RUN}" <<'PYEOF' || rc=$?
import sys, shutil
path, dry_run = sys.argv[1], sys.argv[2] == "1"

with open(path, "r") as fh:
    text = fh.read()

# Check if already patched (window_dt is the telltale — unique post-rename)
# Only declare "already patched" if inv_like_dt is GONE (otherwise partial).
already = ("window_dt" in text) and ("inv_like_dt" not in text)
if already:
    print(f"[skip] already patched: {path}")
    sys.exit(10)   # distinguish from 0 so caller can count skips

# Apply the three renames. Order matters:
#   1. filename renames (most specific first)
#   2. variable rename (broader)
replacements = [
    ("snake_inv_likeness.tsv.gz",   "window_dt.tsv.gz"),
    ("window_inv_likeness.tsv.gz",  "window_dt.tsv.gz"),
    ("inv_like_dt",                 "window_dt"),
]

patched = text
n_hits = 0
for old, new in replacements:
    count = patched.count(old)
    if count > 0:
        patched = patched.replace(old, new)
        n_hits += count

if n_hits == 0:
    print(f"[skip] no target strings found: {path}")
    sys.exit(10)

print(f"[patch] {path}  ({n_hits} substitutions)")
if dry_run:
    sys.exit(0)

shutil.copyfile(path, path + ".bk_chat17")
with open(path, "w") as fh:
    fh.write(patched)
PYEOF
  case $rc in
    0)  N_PATCHED=$((N_PATCHED + 1)) ;;
    10) N_SKIPPED=$((N_SKIPPED + 1)) ;;
    *)  echo "[patch] FAIL: ${rel} (python rc=$rc)" >&2; exit $rc ;;
  esac
done

echo ""
echo "[patch] summary: ${N_PATCHED} patched, ${N_SKIPPED} skipped, ${N_MISSING} missing"
if [[ ${DRY_RUN} -eq 0 && ${N_PATCHED} -gt 0 ]]; then
  echo ""
  echo "[patch] verify:"
  echo "  - All R files still parse: find \$BASE -name '*.R.bk_chat17' | while read bk; do"
  echo "      orig=\${bk%.bk_chat17}; python3 \$BASE/_rcheck.py \$orig; done"
  echo "  - Sanity grep — no inv_like_dt references outside archive:"
  echo "      grep -rln 'inv_like_dt' --include='*.R' \$BASE 2>/dev/null | grep -v _archive | grep -v bk_chat17 | grep -v deprecated"
  echo "      (expect: empty output)"
  echo ""
  echo "[patch] IMPORTANT: existing precomp output files on disk still have"
  echo "        the OLD name (window_inv_likeness.tsv.gz). Either:"
  echo "          (a) rename them to match:"
  echo "              find \$INVDIR -name 'window_inv_likeness.tsv.gz' -exec \\"
  echo "                rename 's/window_inv_likeness/window_dt/' {} +"
  echo "          (b) symlink for back-compat:"
  echo "              for f in \$(find \$INVDIR -name 'window_inv_likeness.tsv.gz'); do"
  echo "                ln -sf \$(basename \$f) \$(dirname \$f)/window_dt.tsv.gz"
  echo "              done"
  echo "          (c) just re-run precomp — which chat 18 will do anyway."
fi
