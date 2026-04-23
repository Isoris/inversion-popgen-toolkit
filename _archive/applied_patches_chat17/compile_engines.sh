#!/usr/bin/env bash
# =============================================================================
# compile_engines.sh
#
# Chat-17 patch — build all C/C++ engine binaries the pipeline depends on.
#
# pipeline_bridge.sh already auto-compiles instant_q + region_popstats +
# hobs_windower on first source, but (a) it silences make errors with
# `2>/dev/null`, making real compile failures hard to diagnose, and
# (b) it does NOT cover the phase-5 engines (rare_sfs_pairwise,
# export_q_residual_dosage). This helper compiles all three makefiles
# explicitly, with verbose output, so a first-time setup gets a clean
# diagnostic if anything fails.
#
# Binaries produced:
#   unified_ancestry/src/instant_q
#     — used by LAUNCH_instant_q_precompute.slurm (ancestry track; chat 16 rewrite)
#
#   unified_ancestry/engines/fst_dxy/region_popstats
#     — used by LAUNCH_region_popstats.slurm (FST/dXY/theta tracks)
#
#   MODULE_5B_inversion_followup/engines/rare_sfs_pairwise
#     — used by LAUNCH_rare_sfs_pairwise.slurm (phase 5 followup)
#   MODULE_5B_inversion_followup/engines/export_q_residual_dosage
#     — used by LAUNCH_q_residual_dosage.slurm (phase 5 followup)
#
# Not covered here:
#   — hobs_windower (already in pipeline_bridge auto-compile, different dir)
#   — any binary not shipped with source (e.g., NGSadmix, ANGSD, etc.)
#
# Usage:
#   bash compile_engines.sh                  # compile all
#   bash compile_engines.sh --clean          # clean all first
#   bash compile_engines.sh --test           # compile + run make test
#   bash compile_engines.sh --only instant_q # compile one
# =============================================================================

set -euo pipefail

BASE="${BASE:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"

DO_CLEAN=0
DO_TEST=0
ONLY=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --clean) DO_CLEAN=1; shift ;;
    --test)  DO_TEST=1;  shift ;;
    --only)  ONLY="$2"; shift 2 ;;
    -h|--help)
      sed -n '2,/^# ===/p' "$0" | sed '/^# ===/d' | sed 's/^# \{0,1\}//'
      exit 0 ;;
    *) echo "unknown arg: $1" >&2; exit 2 ;;
  esac
done

# ── Engine registry ──────────────────────────────────────────────────────
# Each entry: name | makefile_dir | expected_binary_path
ENGINES=(
  "instant_q|${BASE}/unified_ancestry/src|${BASE}/unified_ancestry/src/instant_q"
  "region_popstats|${BASE}/unified_ancestry/engines/fst_dxy|${BASE}/unified_ancestry/engines/fst_dxy/region_popstats"
  "phase5_engines|${BASE}/inversion_modules/phase_5_followup/engines|${BASE}/inversion_modules/phase_5_followup/engines/rare_sfs_pairwise"
)
# Phase-5 Makefile builds TWO binaries; tracking one is enough for existence
# check since make all builds both.

# hobs_windower has a .c file but no Makefile — compile directly.
# Matches the gcc line in utils/pipeline_bridge.sh _compile_if_missing.
HOBS_SRC="${BASE}/unified_ancestry/engines/hobs_hwe/scripts/hobs_windower.c"
HOBS_BIN="${BASE}/unified_ancestry/engines/hobs_hwe/scripts/hobs_windower"

ok_count=0
fail_count=0
skipped=0

for entry in "${ENGINES[@]}"; do
  IFS='|' read -r name mkdir binpath <<< "$entry"

  if [[ -n "$ONLY" && "$ONLY" != "$name" ]]; then
    continue
  fi

  echo "================================================================"
  echo "  $name"
  echo "  Makefile dir: $mkdir"
  echo "  Binary:       $binpath"
  echo "================================================================"

  if [[ ! -d "$mkdir" ]]; then
    echo "  SKIP: makefile dir missing"
    skipped=$((skipped + 1))
    continue
  fi
  if [[ ! -f "${mkdir}/Makefile" ]]; then
    echo "  SKIP: no Makefile in $mkdir"
    skipped=$((skipped + 1))
    continue
  fi

  pushd "$mkdir" > /dev/null

  if [[ $DO_CLEAN -eq 1 ]]; then
    echo "  make clean..."
    make clean || true
  fi

  echo "  make..."
  if make; then
    if [[ $DO_TEST -eq 1 ]]; then
      echo "  make test..."
      make test || {
        echo "  FAIL: make test"
        fail_count=$((fail_count + 1))
        popd > /dev/null
        continue
      }
    fi
    ok_count=$((ok_count + 1))
  else
    echo "  FAIL: make returned non-zero"
    fail_count=$((fail_count + 1))
  fi

  popd > /dev/null
  echo ""
done

# ── hobs_windower — no Makefile, raw gcc ────────────────────────────────
if [[ -z "$ONLY" || "$ONLY" == "hobs_windower" ]]; then
  echo "================================================================"
  echo "  hobs_windower (no Makefile — direct gcc)"
  echo "  Source: $HOBS_SRC"
  echo "  Binary: $HOBS_BIN"
  echo "================================================================"
  if [[ ! -f "$HOBS_SRC" ]]; then
    echo "  SKIP: source missing"
    skipped=$((skipped + 1))
  else
    if [[ $DO_CLEAN -eq 1 ]]; then
      rm -f "$HOBS_BIN"
    fi
    echo "  gcc -O2..."
    if gcc -O2 -o "$HOBS_BIN" "$HOBS_SRC" -lm -lz -lpthread; then
      ok_count=$((ok_count + 1))
    else
      echo "  FAIL: gcc returned non-zero"
      fail_count=$((fail_count + 1))
    fi
  fi
  echo ""
fi

# ── Dead-wood warning ────────────────────────────────────────────────────
# Stale copies of rare_sfs_pairwise.c and export_q_residual_dosage.c exist
# under unified_ancestry/engines/fst_dxy/ even though the fst_dxy Makefile
# explicitly no longer builds them (comment at top of that Makefile). The
# canonical copies are in inversion_modules/phase_5_followup/engines/.
#
# We don't delete them (not our call), but flag if both sets exist so the
# user doesn't accidentally build the wrong copy later.
STALE_FST_SRCS=(
  "${BASE}/unified_ancestry/engines/fst_dxy/rare_sfs_pairwise.c"
  "${BASE}/unified_ancestry/engines/fst_dxy/export_q_residual_dosage.c"
)
stale_present=0
for f in "${STALE_FST_SRCS[@]}"; do
  [[ -f "$f" ]] && stale_present=1
done
if [[ $stale_present -eq 1 ]]; then
  echo "[compile_engines] NOTE: stale source copies exist under"
  echo "                   unified_ancestry/engines/fst_dxy/ that the Makefile"
  echo "                   explicitly does NOT build. Canonical copies are in"
  echo "                   inversion_modules/phase_5_followup/engines/."
  echo "                   Consider cleaning up after HPC validation."
fi

# ── Verify binaries exist + are executable ───────────────────────────────
echo "================================================================"
echo "  Binary status check"
echo "================================================================"
for entry in "${ENGINES[@]}"; do
  IFS='|' read -r name mkdir binpath <<< "$entry"
  [[ -n "$ONLY" && "$ONLY" != "$name" ]] && continue
  if [[ -x "$binpath" ]]; then
    echo "  OK   $binpath"
  else
    echo "  MISS $binpath"
  fi
done

# Also check the phase-5 second binary
if [[ -z "$ONLY" || "$ONLY" == "phase5_engines" ]]; then
  bin2="${BASE}/inversion_modules/phase_5_followup/engines/export_q_residual_dosage"
  if [[ -x "$bin2" ]]; then
    echo "  OK   $bin2"
  else
    echo "  MISS $bin2"
  fi
fi

if [[ -z "$ONLY" || "$ONLY" == "hobs_windower" ]]; then
  if [[ -x "$HOBS_BIN" ]]; then
    echo "  OK   $HOBS_BIN"
  else
    echo "  MISS $HOBS_BIN"
  fi
fi

echo ""
echo "[compile_engines] summary: ${ok_count} built, ${fail_count} failed, ${skipped} skipped"
[[ $fail_count -eq 0 ]] || exit 1
