#!/usr/bin/env bash
# =============================================================================
# patch_rename_scientific_names.sh
#
# Chat-17 patch. Applies Tier A of the RENAMING.md plan that was only
# partially executed — internal identifiers still carry informal nicknames
# from the pipeline's early development even though the flagship filenames
# were already renamed. See NAMING_AUDIT.md (in this patches dir) for the
# full triage and rationale for what IS and IS NOT in this patch.
#
# RENAMES (literal, whole-word-ish via surrounding context):
#
#   flashlight          → sv_prior              (concept term / var / CLI)
#   FLASHLIGHT_LOADER   → SV_PRIOR_LOADER       (env var)
#   --flashlight_dir    → --sv_prior_dir        (CLI flag)
#   --flashlight        → --sv_prior            (CLI flag, short form)
#   snake_id            → region_id             (column / field / param)
#   snake_phase         → extension_phase       (column in decision log)
#   core_family         → scale_tier            (column / field)
#
# NOT RENAMED HERE (see NAMING_AUDIT.md Tier B/C):
#   - Filenames (cheat*.R, *_snake3_*) — require coordinated launcher updates
#   - Directory `snake_regions_multiscale/` — touches many launchers
#   - Column `inv_likeness` (accurate name for the score column itself)
#   - `peel` / `peeling` (documented scientifically in-place)
#   - `STAIR_*`, `NN_*`, `CONSENSUS_*`, `LANDSCAPE_*` (function-describing)
#
# BACK-COMPAT PRESERVED:
#   The C01g boundary script reads `sv_flashlight_<chr>.rds` as a fallback
#   for legacy files. That fallback READER is kept as-is. This patch only
#   renames the variable names, CLI flags, and log prefixes — not the
#   fallback-filename string literal.
#
# NOTE on chat-17 launchers:
#   LAUNCH_C01g_boundary.sh currently passes `--flashlight <dir>` to the
#   R script. If you apply this patch, the launcher must also be updated
#   to pass `--sv_prior` instead. The patch handles the R-side; launcher
#   fix is one `sed` line noted in the post-patch verification block.
#
# Backup suffix: .bk_chat17_scientific
# Usage:  bash patch_rename_scientific_names.sh [--dry-run]
# =============================================================================

set -euo pipefail

BASE="${BASE:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"
DRY_RUN=0
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=1

# List of files (relative to $BASE) — audited from the chat-16 tarball.
# Re-audit on HPC before applying if the tree has drifted.
FILES=(
  "inversion_modules/phase_2_discovery/2c_precomp/STEP_C00_build_sv_prior.R"
  "inversion_modules/phase_2_discovery/2c_precomp/STEP_C01a_precompute.R"
  "inversion_modules/phase_2_discovery/2c_precomp/patches/patch_C01a_precompute_flashlight.R"
  "inversion_modules/phase_2_discovery/2c_precomp/patches/patch_C01b2_merge_flashlight.R"
  "inversion_modules/phase_2_discovery/2c_precomp/patches/patch_C01d_scoring_flashlight.R"
  "inversion_modules/phase_2_discovery/2c_precomp/patches/patch_C01f_hypothesis_flashlight.R"
  "inversion_modules/phase_2_discovery/2c_precomp/patches/patch_C01i_decomposition_flashlight.R"
  "inversion_modules/phase_2_discovery/2c_precomp/patches/patch_C04_ghsl_flashlight.R"
  "inversion_modules/phase_4_postprocessing/4a_existence_layers/STEP_C01d_candidate_scoring_wired_25_v934_registry.R"
  "inversion_modules/phase_4_postprocessing/4a_existence_layers/STEP_C01g_boundary_catalog_wired_4_8_10_11_17_21_v934_registry.R"
  "inversion_modules/phase_4_postprocessing/4a_existence_layers/STEP_C01l_local_structure_segments.R"
  "inversion_modules/phase_4_postprocessing/4b_group_proposal/STEP_C01i_b_multi_recomb.R"
  "inversion_modules/phase_4_postprocessing/4b_group_proposal/STEP_C01i_d_seal.R"
  "inversion_modules/phase_4_postprocessing/4b_group_proposal/STEP_C01i_decompose.R"
  "inversion_modules/phase_4_postprocessing/4b_group_proposal/lib_decompose_helpers.R"
  "inversion_modules/phase_4_postprocessing/4b_group_proposal/lib_recomb_combination.R"
  "inversion_modules/phase_4_postprocessing/4b_group_proposal/lib_step03_seed_loader.R"
  "inversion_modules/phase_4_postprocessing/4e_final_classification/compute_candidate_status.R"
  "inversion_modules/phase_5_followup/codebase_8.5_experimental/current_followup/STEP40_candidate_internal_coherence.R"
  "inversion_modules/utils/debug_core_pcas.R"
)

N_PATCHED=0
N_SKIPPED=0
N_MISSING=0

for rel in "${FILES[@]}"; do
  path="${BASE}/${rel}"
  if [[ ! -f "${path}" ]]; then
    echo "[rename] MISSING: ${rel}" >&2
    N_MISSING=$((N_MISSING + 1))
    continue
  fi

  # Capture python exit code safely under `set -e`. The `|| rc=$?` pattern
  # lets us read the non-zero exit without the whole script dying.
  rc=0
  python3 - "${path}" "${DRY_RUN}" <<'PYEOF' || rc=$?
import sys, shutil, re
path, dry_run = sys.argv[1], sys.argv[2] == "1"

with open(path, "r") as fh:
    text = fh.read()

# Ordered list of renames. CLI flags first (most specific), then uppercase
# env vars, then variable/concept names. Using regex with word-boundary
# anchors to avoid clobbering substrings in comments or paths.
#
# Order matters: --flashlight_dir must be replaced BEFORE --flashlight
# so the short form doesn't clobber the long-form match.
RENAMES = [
    # CLI flags (longest first)
    (r"--flashlight_dir", "--sv_prior_dir"),
    (r"--flashlight\b",   "--sv_prior"),
    # Uppercase env vars
    (r"\bFLASHLIGHT_LOADER\b", "SV_PRIOR_LOADER"),
    (r"\bFLASH_DIR\b",         "SV_PRIOR_DIR"),
    # Concept term (case-sensitive lowercase)
    (r"\bflashlight_dir\b",  "sv_prior_dir"),
    (r"\bflashlight\b",      "sv_prior"),
    # Snake/core renames
    (r"\bsnake_id\b",        "region_id"),
    (r"\bsnake_phase\b",     "extension_phase"),
    (r"\bcore_family\b",     "scale_tier"),
    # Log prefix
    (r"\[flashlight\]", "[sv_prior]"),
]

# Idempotency check: if none of the OLD tokens match anymore, skip.
any_hit = False
n_hits_per_pattern = []
for pat, new in RENAMES:
    n = len(re.findall(pat, text))
    n_hits_per_pattern.append((pat, n))
    if n > 0:
        any_hit = True

if not any_hit:
    print(f"[skip] already clean: {path}")
    sys.exit(10)

patched = text
total_hits = 0
for pat, new in RENAMES:
    patched, n = re.subn(pat, new, patched)
    total_hits += n

# Preserve one known-safe legacy reference: the fallback reader for
# sv_flashlight_<chr>.rds in C01g should NOT be renamed because that file
# may still exist on disk from prior runs. Check that we didn't clobber it.
if "sv_flashlight_" in text and "sv_flashlight_" not in patched:
    print(f"[patch] ERROR: clobbered back-compat 'sv_flashlight_' in {path}", file=sys.stderr)
    sys.exit(3)

print(f"[rename] {path}  ({total_hits} substitutions)")
if dry_run:
    for pat, n in n_hits_per_pattern:
        if n > 0:
            print(f"           {pat}: {n}")
    sys.exit(0)

shutil.copyfile(path, path + ".bk_chat17_scientific")
with open(path, "w") as fh:
    fh.write(patched)
PYEOF
  case $rc in
    0)  N_PATCHED=$((N_PATCHED + 1)) ;;
    10) N_SKIPPED=$((N_SKIPPED + 1)) ;;
    *)  echo "[rename] FAIL: ${rel} (python rc=$rc)" >&2; exit $rc ;;
  esac
done

echo ""
echo "[rename] summary: ${N_PATCHED} patched, ${N_SKIPPED} already clean, ${N_MISSING} missing"

if [[ ${DRY_RUN} -eq 0 && ${N_PATCHED} -gt 0 ]]; then
  echo ""
  echo "[rename] Also patch the chat-17 LAUNCH_C01g_boundary.sh (it passes --flashlight):"
  echo "  sed -i.bk_chat17_scientific \\"
  echo "    -e 's|--flashlight  \"|--sv_prior  \"|g' \\"
  echo "    \${BASE}/inversion_modules/phase_4_postprocessing/4a_existence_layers/LAUNCH_C01g_boundary.sh"
  echo ""
  echo "[rename] Verify:"
  echo "  # All R files still parse"
  echo "  for bk in \$(find \${BASE} -name '*.bk_chat17_scientific'); do"
  echo "    orig=\${bk%.bk_chat17_scientific}"
  echo "    [[ \$orig == *.R ]] && python3 \${BASE}/_rcheck.py \"\$orig\""
  echo "  done"
  echo ""
  echo "  # No 'flashlight' references outside archive/doc"
  echo "  grep -rln 'flashlight' --include='*.R' --include='*.py' --include='*.sh' \\"
  echo "    \${BASE} 2>/dev/null | grep -v _archive | grep -v deprecated | \\"
  echo "    grep -v RENAMING.md | grep -v NAMING_AUDIT.md | grep -v bk_chat17"
fi
