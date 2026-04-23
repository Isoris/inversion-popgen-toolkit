#!/usr/bin/env bash
# =============================================================================
# install_update_2026-04-20.sh
# -----------------------------------------------------------------------------
# Applies the toolkit_update_2026-04-20 bundle to your live toolkit.
#
# Designed for safety: dry-run by default, --apply to really do it. Backups
# of every overwritten file. No destructive action until you say --apply.
#
# Run from: wherever you extracted the bundle (i.e. alongside
#           toolkit_update_2026-04-20/, inversion-popgen-toolkit/, and
#           inversion_modules/).
#
# Usage:
#   bash install_update_2026-04-20.sh                # preview (dry run)
#   bash install_update_2026-04-20.sh --apply        # actually apply
#   bash install_update_2026-04-20.sh --apply --skip-v21    # skip v2.1 patcher
#                                                           # (if already done)
# =============================================================================
set -euo pipefail

# --- Paths (edit these if you're in a different layout) ---------------------
BUNDLE="${PWD}/toolkit_update_2026-04-20"
LIVE_TOOLKIT="${PWD}/inversion-popgen-toolkit"
LIVE_MODULE="${PWD}/inversion-popgen-toolkit/inversion_modules/phase_4_postprocessing/4b_qc_triage"

# --- Flags ------------------------------------------------------------------
APPLY=false
SKIP_V21=false
for arg in "$@"; do
  case "$arg" in
    --apply)    APPLY=true ;;
    --skip-v21) SKIP_V21=true ;;
    *) echo "unknown arg: $arg"; exit 1 ;;
  esac
done

# --- Sanity --------------------------------------------------------------
[[ -d "$BUNDLE" ]]         || { echo "ERROR: $BUNDLE not found. Run from the dir that contains it."; exit 1; }
[[ -d "$LIVE_TOOLKIT" ]]   || { echo "ERROR: $LIVE_TOOLKIT not found."; exit 1; }
[[ -d "$LIVE_MODULE" ]]    || { echo "ERROR: $LIVE_MODULE not found."; exit 1; }

ts="$(date +%Y%m%d_%H%M%S)"
BACKUP_DIR="${PWD}/.install_backup_${ts}"
LOG_FILE="${PWD}/install_update_2026-04-20.${ts}.log"

say() { printf '%s\n' "$*" | tee -a "$LOG_FILE"; }
do_cmd() {
  # Run a command; echo it; only actually run if --apply.
  say "  $ $*"
  if [[ "$APPLY" == true ]]; then
    eval "$@"
  fi
}

say "================================================================"
say "  Install update 2026-04-20"
say "  Mode: $([[ "$APPLY" == true ]] && echo "APPLY" || echo "DRY RUN (preview only)")"
say "  Bundle: $BUNDLE"
say "  Live toolkit: $LIVE_TOOLKIT"
say "  Live module:  $LIVE_MODULE"
say "  Backups:      $BACKUP_DIR"
say "  Log:          $LOG_FILE"
say "================================================================"

# =============================================================================
# Step 0: Check if v2.1 is already applied (gives us info for step 1)
# =============================================================================
say ""
say "--- Step 0: Check whether v2.1 is already applied ---"

V21_APPLIED=false
if [[ -f "${LIVE_MODULE}/R/q09_gap_characterization.R" ]] && \
   [[ -f "${LIVE_MODULE}/STEP_Q09_gap_characterization.sh" ]] && \
   [[ -f "${LIVE_MODULE}/run_chrom.sh" ]]; then
  V21_APPLIED=true
  say "  ✓ v2.1 signature files present — v2.1 appears installed"
else
  say "  ✗ v2.1 signature files missing — v2.1 NOT yet installed"
  [[ -f "${LIVE_MODULE}/R/q09_gap_characterization.R" ]] || say "    missing: R/q09_gap_characterization.R"
  [[ -f "${LIVE_MODULE}/STEP_Q09_gap_characterization.sh" ]] || say "    missing: STEP_Q09_gap_characterization.sh"
  [[ -f "${LIVE_MODULE}/run_chrom.sh" ]] || say "    missing: run_chrom.sh"
fi

# =============================================================================
# Step 1: Run the v2.1 upgrade patcher (if needed)
# =============================================================================
say ""
say "--- Step 1: v2.1 patcher ---"

if [[ "$SKIP_V21" == true ]]; then
  say "  SKIPPED (--skip-v21 flag set)"
elif [[ "$V21_APPLIED" == true ]]; then
  say "  SKIPPED (v2.1 already applied per step 0)"
else
  say "  Need to apply UPGRADE_v2.1.sh"
  if [[ "$APPLY" == true ]]; then
    say "  Running UPGRADE_v2.1.sh --dry-run first:"
    (cd "$LIVE_MODULE" && bash "${BUNDLE}/UPGRADE_v2.1.sh" --dry-run 2>&1 | sed 's/^/    /' | tee -a "$LOG_FILE")
    say ""
    say "  Now applying UPGRADE_v2.1.sh (real run):"
    (cd "$LIVE_MODULE" && bash "${BUNDLE}/UPGRADE_v2.1.sh" 2>&1 | sed 's/^/    /' | tee -a "$LOG_FILE")
  else
    say "  Would run: cd $LIVE_MODULE && bash ${BUNDLE}/UPGRADE_v2.1.sh"
  fi
fi

# =============================================================================
# Step 2: Overlay the 2 standardized files (00_config.sh + Q01) with backups
# =============================================================================
say ""
say "--- Step 2: Overlay standardized 00_config.sh + STEP_Q01 ---"

if [[ "$APPLY" == true ]]; then
  mkdir -p "$BACKUP_DIR"
fi

for name in 00_config.sh STEP_Q01_snp_density.sh; do
  src="${BUNDLE}/modules/phase_qc_shelf/${name}"
  dst="${LIVE_MODULE}/${name}"
  if [[ ! -f "$src" ]]; then
    say "  ❌ bundle missing $src"
    continue
  fi
  if [[ -f "$dst" ]]; then
    size_old=$(stat -c%s "$dst" 2>/dev/null || echo "?")
    size_new=$(stat -c%s "$src" 2>/dev/null || echo "?")
    md5_old=$(md5sum "$dst" 2>/dev/null | awk '{print $1}')
    md5_new=$(md5sum "$src" 2>/dev/null | awk '{print $1}')
    if [[ "$md5_old" == "$md5_new" ]]; then
      say "  = IDENTICAL  $name  (no change needed)"
    else
      say "  ≠ WILL UPDATE  $name  (old=${size_old}b → new=${size_new}b)"
      do_cmd "cp '$dst' '$BACKUP_DIR/$name.backup'"
      do_cmd "cp '$src' '$dst'"
      do_cmd "bash -n '$dst' && echo '    syntax OK' || echo '    ⚠️ syntax FAIL'"
    fi
  else
    say "  + WILL CREATE  $name  (not currently in live module)"
    do_cmd "cp '$src' '$dst'"
  fi
done

# =============================================================================
# Step 3: Install tools/ and docs/ at toolkit root (new dirs, no overlay risk)
# =============================================================================
say ""
say "--- Step 3: Install tools/ and docs/ at toolkit root ---"

for sub in tools docs; do
  src_dir="${BUNDLE}/${sub}"
  dst_dir="${LIVE_TOOLKIT}/${sub}"
  [[ -d "$src_dir" ]] || { say "  ❌ bundle missing $src_dir"; continue; }

  if [[ ! -d "$dst_dir" ]]; then
    say "  + WILL CREATE  $dst_dir/"
    do_cmd "mkdir -p '$dst_dir'"
  else
    say "  = EXISTS       $dst_dir/  (will overlay files)"
  fi

  for f in "$src_dir"/*; do
    name=$(basename "$f")
    target="$dst_dir/$name"
    if [[ -f "$target" ]]; then
      md5_old=$(md5sum "$target" 2>/dev/null | awk '{print $1}')
      md5_new=$(md5sum "$f" 2>/dev/null | awk '{print $1}')
      if [[ "$md5_old" == "$md5_new" ]]; then
        say "    = identical    $sub/$name"
      else
        say "    ≠ will update  $sub/$name"
        do_cmd "cp '$target' '$BACKUP_DIR/$sub.$name.backup' 2>/dev/null || true"
        do_cmd "cp '$f' '$target'"
      fi
    else
      say "    + will create  $sub/$name"
      do_cmd "cp '$f' '$target'"
    fi
  done

  # chmod python tools
  if [[ "$sub" == "tools" ]]; then
    do_cmd "chmod +x '$dst_dir'/*.py 2>/dev/null || true"
  fi
done

# =============================================================================
# Step 4: Summary
# =============================================================================
say ""
say "================================================================"
if [[ "$APPLY" == true ]]; then
  say "  ✅ Applied. Backups in: $BACKUP_DIR"
  say "     Log in: $LOG_FILE"
  say ""
  say "  Next steps:"
  say "    1. Test Q01:"
  say "         cd $LIVE_MODULE"
  say "         bash STEP_Q01_snp_density.sh C_gar_LG28 2>&1 | head -40"
  say "       Look for: 🧾 config snapshot, ▶ banner, 📥 📤 file previews, 📊 summary, ✅ done"
  say ""
  say "    2. Git commit the changes:"
  say "         cd $LIVE_TOOLKIT"
  say "         git status"
  say "         git add -A"
  say "         git commit -m 'Update 2026-04-20: phase_qc_shelf v2.1 + logging helpers + audit docs'"
else
  say "  🔎 DRY RUN complete. Nothing changed on disk."
  say "     Review the actions above. If you agree, re-run with --apply."
fi
say "================================================================"
