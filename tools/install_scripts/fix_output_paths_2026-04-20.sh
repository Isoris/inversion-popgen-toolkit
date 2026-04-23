#!/usr/bin/env bash
# =============================================================================
# fix_output_paths_2026-04-20.sh
# -----------------------------------------------------------------------------
# Relocates phase_qc_shelf outputs from the misplaced top-level folder
# ($BASE/inversion_modules/phase_4_postprocessing/4b_qc_triage/results/) to the correct location
# ($BASE/inversion_localpca_v7/phase_qc_shelf_results/) and fixes the default
# in 00_config.sh so future runs land in the right place.
#
# Does THREE things:
#   1. Patch 00_config.sh (live module + bundle copy) so the QC_OUT default
#      points to inversion_localpca_v7/phase_qc_shelf_results
#   2. Move the existing results/ tree to the new location
#   3. Remove the now-empty $BASE/inversion_modules/ mistake folder
#
# Safe mode by default: preview only. Pass --apply to really do it.
#
# Usage:
#   bash fix_output_paths_2026-04-20.sh                # preview only
#   bash fix_output_paths_2026-04-20.sh --apply        # actually apply
# =============================================================================
set -euo pipefail

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
LIVE_CONFIG="${BASE}/inversion-popgen-toolkit/inversion_modules/phase_4_postprocessing/4b_qc_triage/00_config.sh"
BUNDLE_CONFIG="${BASE}/toolkit_update_2026-04-20/modules/phase_qc_shelf/00_config.sh"

OLD_RESULTS="${BASE}/inversion_modules/phase_4_postprocessing/4b_qc_triage/results"
NEW_RESULTS="${BASE}/inversion_localpca_v7/phase_qc_shelf_results"
OLD_MISTAKE_DIR="${BASE}/inversion_modules"

APPLY=false
[[ "${1:-}" == "--apply" ]] && APPLY=true

ts="$(date +%Y%m%d_%H%M%S)"
LOG="${BASE}/fix_output_paths.${ts}.log"
BACKUP_DIR="${BASE}/.fix_output_paths_backup_${ts}"

say() { printf '%s\n' "$*" | tee -a "$LOG"; }
show() {
  say "  $ $*"
  if [[ "$APPLY" == true ]]; then
    eval "$@"
  fi
}

say "================================================================"
say "  fix_output_paths_2026-04-20"
say "  Mode: $([[ "$APPLY" == true ]] && echo "APPLY" || echo "DRY RUN")"
say "  Log:  $LOG"
say "================================================================"

# =============================================================================
# Step 1: Patch 00_config.sh files (live + bundle)
# =============================================================================
say ""
say "--- Step 1: Patch QC_OUT default in 00_config.sh ---"

OLD_LINE=': "${QC_OUT:=${BASE}/inversion_modules/phase_4_postprocessing/4b_qc_triage/results}"'
NEW_LINE=': "${QC_OUT:=${BASE}/inversion_localpca_v7/phase_qc_shelf_results}"'

for target in "$LIVE_CONFIG" "$BUNDLE_CONFIG"; do
  if [[ ! -f "$target" ]]; then
    say "  ⚠️  missing: $target"
    continue
  fi

  if grep -qF "$OLD_LINE" "$target"; then
    say "  WILL PATCH  $target"
    say "    old: $OLD_LINE"
    say "    new: $NEW_LINE"
    if [[ "$APPLY" == true ]]; then
      mkdir -p "$BACKUP_DIR"
      cp "$target" "$BACKUP_DIR/$(basename "$target").$(echo "$target" | tr '/' '_').backup"
    fi
    # Escape for sed: | delimiter avoids the / in paths
    show "sed -i 's|${OLD_LINE}|${NEW_LINE}|' '$target'"
  elif grep -qF "$NEW_LINE" "$target"; then
    say "  = already patched  $target"
  else
    say "  ⚠️  pattern not found in $target — the file may differ from expected"
    say "    look for a QC_OUT default and update it manually to point at"
    say "    \${BASE}/inversion_localpca_v7/phase_qc_shelf_results"
  fi
done

# =============================================================================
# Step 2: Move existing results to the new location
# =============================================================================
say ""
say "--- Step 2: Move existing results/ to new location ---"

if [[ ! -d "$OLD_RESULTS" ]]; then
  say "  ✓ no existing results to move (already done?)"
else
  old_size=$(du -sh "$OLD_RESULTS" 2>/dev/null | awk '{print $1}')
  say "  source:   $OLD_RESULTS  ($old_size)"
  say "  target:   $NEW_RESULTS"

  if [[ -d "$NEW_RESULTS" ]]; then
    say "  ⚠️  target exists already"
    say "  files currently in target:"
    ls -la "$NEW_RESULTS" 2>/dev/null | head -10 | sed 's/^/      /' | tee -a "$LOG"
    say ""
    say "  Strategy: move each subdir (figures/ logs/ tracks/ popstats_groups/)"
    say "  individually; if a subdir name clashes, add a .YYYYMMDD suffix."

    for sub in figures logs tracks popstats_groups; do
      src="$OLD_RESULTS/$sub"
      dst="$NEW_RESULTS/$sub"
      if [[ ! -d "$src" ]]; then
        say "    - skip  $sub/  (not present in old location)"
        continue
      fi
      if [[ -e "$dst" ]]; then
        say "    ⚠️  CONFLICT  $dst already exists"
        say "       preview: ls of src:"
        ls "$src" | head -5 | sed 's/^/         /' | tee -a "$LOG"
        say "       preview: ls of dst:"
        ls "$dst" | head -5 | sed 's/^/         /' | tee -a "$LOG"
        say "       would rename src -> ${dst}.preserved_${ts}/"
        show "mv '$src' '${dst}.preserved_${ts}'"
      else
        show "mkdir -p '$NEW_RESULTS'"
        show "mv '$src' '$dst'"
      fi
    done
  else
    show "mkdir -p '$(dirname "$NEW_RESULTS")'"
    show "mv '$OLD_RESULTS' '$NEW_RESULTS'"
  fi
fi

# =============================================================================
# Step 3: Remove the mistake folder if it's now empty
# =============================================================================
say ""
say "--- Step 3: Clean up mistake folder ---"

if [[ ! -d "$OLD_MISTAKE_DIR" ]]; then
  say "  ✓ already gone"
else
  # Count remaining non-empty contents
  remaining=$(find "$OLD_MISTAKE_DIR" -mindepth 1 ! -empty 2>/dev/null | wc -l)
  if [[ "$remaining" == "0" ]]; then
    say "  WILL REMOVE  $OLD_MISTAKE_DIR  (empty — safe)"
    show "rm -rf '$OLD_MISTAKE_DIR'"
  else
    say "  ⚠️  NOT EMPTY — preserved for manual review"
    say "  Contents:"
    find "$OLD_MISTAKE_DIR" -mindepth 1 2>/dev/null | head -15 | sed 's/^/      /' | tee -a "$LOG"
    say "  (Not removing. Inspect manually and decide.)"
  fi
fi

# =============================================================================
# Summary
# =============================================================================
say ""
say "================================================================"
if [[ "$APPLY" == true ]]; then
  say "  ✅ APPLIED."
  say "  Backup of config files: $BACKUP_DIR/"
  say "  Log: $LOG"
  say ""
  say "  Verify:"
  say "    ls $NEW_RESULTS"
  say "    grep QC_OUT $LIVE_CONFIG"
  say "    ls $OLD_MISTAKE_DIR 2>/dev/null && echo 'mistake dir still there' || echo 'mistake dir gone'"
else
  say "  🔎 DRY RUN only. Nothing changed on disk."
  say "  Re-run with --apply to actually do it."
fi
say "================================================================"
