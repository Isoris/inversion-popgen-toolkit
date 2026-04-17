#!/usr/bin/env bash
# =============================================================================
# patch_rename_scientific_names_v2.sh
#
# Chat-17 follow-on to patch_rename_scientific_names.sh. The v1 patch used
# whole-word regex `\bflashlight\b` which (correctly) did NOT match compound
# identifiers like `flashlight_seeds`, `try_load_flashlight`, `.has_flashlight`,
# or status strings like `"flashlight_only"`. This v2 catches those.
#
# RENAMES applied here:
#
#   flashlight_(\w+)    → sv_prior_\1           (function/var prefix form)
#   (\w+)_flashlight\b  → \1_sv_prior           (function/var suffix form)
#   "flashlight_only"   → "sv_prior_only"       (status string literals)
#   "flashlight_sparse" → "sv_prior_sparse"
#   "flashlight+step03" → "sv_prior+step03"
#   "flashlight"        → "sv_prior"            (bare status string)
#
# Examples of what this will rename:
#   flashlight_seeds     → sv_prior_seeds
#   flashlight_seeded    → sv_prior_seeded
#   flashlight_available → sv_prior_available
#   flashlight_path      → sv_prior_path
#   try_load_flashlight  → try_load_sv_prior
#   load_flashlight      → load_sv_prior
#   .has_flashlight      → .has_sv_prior
#   s3_flashlight_hemi   → s3_sv_prior_hemi
#
# EXPLICITLY PRESERVED (hardcoded path strings to a stale directory that
# isn't present in the tree; renaming won't fix that they point at
# nothing; that's a separate issue for a future chat):
#
#   "flashlight_v2/cache"             → keep as-is (broken path, not ours)
#   "flashlight_v2", "cheats"         → keep
#   "flashlight_v2", "utils"          → keep
#   "flashlight_loader_v2.R"          → keep (source file that doesn't exist)
#   "sv_flashlight_<chr>.rds"         → keep (legacy on-disk filename fallback)
#
# The regex uses a negative lookbehind + lookahead to avoid matching
# the `flashlight_v2` directory name and the `sv_flashlight_` filename
# prefix.
#
# Backup suffix: .bk_chat17_scientific_v2
# Usage:  bash patch_rename_scientific_names_v2.sh [--dry-run]
# =============================================================================

set -euo pipefail

BASE="${BASE:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"
DRY_RUN=0
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=1

# Audit all *.R and *.py files in the active tree and touch only ones that
# contain one of the compound patterns. Let Python do the detection to avoid
# shell regex-escape hell.
mapfile -t FILES < <(
  grep -rln 'flashlight_[a-zA-Z]\|_flashlight\b\|"flashlight"\|"flashlight_only"\|"flashlight_sparse"\|"flashlight+step03"' \
    --include="*.R" --include="*.py" --include="*.sh" "${BASE}" 2>/dev/null \
    | grep -v _archive | grep -v deprecated | grep -v bk_chat17 \
    | grep -v RENAMING.md | grep -v _bk_rename.py || true
)

if [[ ${#FILES[@]} -eq 0 ]]; then
  echo "[rename-v2] no files contain v2-target patterns; nothing to do"
  exit 0
fi

echo "[rename-v2] will scan ${#FILES[@]} file(s):"
for f in "${FILES[@]}"; do echo "  - ${f#$BASE/}"; done
echo ""

N_PATCHED=0
N_SKIPPED=0

for path in "${FILES[@]}"; do
  rc=0
  python3 - "${path}" "${DRY_RUN}" <<'PYEOF' || rc=$?
import sys, shutil, re
path, dry_run = sys.argv[1], sys.argv[2] == "1"

with open(path, "r") as fh:
    text = fh.read()

original = text

# Apply renames in order. Each uses a regex that is CAREFUL to not
# clobber the known-preserved strings (flashlight_v2, sv_flashlight_,
# flashlight_loader_v2).

# Rule 1: flashlight_X where X is an identifier continuation, but NOT
# flashlight_v2 (the stale dir), NOT flashlight_loader_v2 (stale file),
# NOT flashlight_only/sparse (handled as status strings below).
#
# Negative lookahead: (?!v2\b|loader_v2|only\b|sparse\b|seeds\b)
# Actually we WANT to rename flashlight_seeds — that's a real identifier.
# So the skip list is just: v2 and loader_v2.
def rename_prefix(m):
    tail = m.group(1)
    if tail.startswith("v2") or tail.startswith("loader_v2"):
        return m.group(0)  # keep as-is
    return f"sv_prior_{tail}"

# Match flashlight_ followed by one or more word chars
text = re.sub(r"\bflashlight_(\w+)", rename_prefix, text)

# Rule 2: X_flashlight where X is an identifier prefix, but NOT
# sv_flashlight (the legacy filename prefix).
def rename_suffix(m):
    head = m.group(1)
    if head == "sv":
        return m.group(0)  # keep sv_flashlight_ as-is
    return f"{head}_sv_prior"

# Match word chars followed by _flashlight at a word boundary
# Need to be careful: _flashlight could be part of a larger identifier
# like try_load_flashlight. We want the ENTIRE preceding identifier.
# Use a lazy match with \w+ and require it's at a word boundary start.
text = re.sub(r"\b(\w+?)_flashlight\b", rename_suffix, text)

# Rule 2b: flashlight in the MIDDLE of a compound identifier.
# e.g. s3_flashlight_hemi — neither prefix nor suffix rule catches this
# because flashlight is surrounded by underscores (word chars on both
# sides, no \b). Match _flashlight_ anywhere and rename to _sv_prior_,
# but preserve the bounded cases flashlight_v2 and loader_v2 suffix.
def rename_middle(m):
    prefix = m.group(1)  # what came before _flashlight_
    suffix = m.group(2)  # what comes after _flashlight_
    # Skip the stale-directory patterns
    if suffix.startswith("v2") or suffix.startswith("loader_v2"):
        return m.group(0)
    return f"{prefix}_sv_prior_{suffix}"

# Match (identifier)_flashlight_(identifier)
text = re.sub(r"(\w+)_flashlight_(\w+)", rename_middle, text)

# Rule 3: status string literals. These have already been caught by Rule 1
# if they were bare identifiers (flashlight_only → sv_prior_only), but in
# quoted string form the regex also applies because \w doesn't distinguish
# quotes. Let's verify:
# "flashlight_only" — matches flashlight_only via rule 1 (tail = "only\"...")
# Actually the tail would include the closing quote? No: \w+ only matches
# word chars, so flashlight_only in "flashlight_only" matches just _only.
# So rule 1 handles them already. Good.

# Rule 4: bare "flashlight" status string (from lib_step03_seed_loader.R
# line 178: source = "flashlight"). This is a literal value, not an
# identifier. Rule 1 doesn't touch it (no underscore follows). Rule 2
# doesn't touch it (no underscore precedes). Handle explicitly.
# Match "flashlight" as a complete quoted string (flanked by ").
text = re.sub(r'"flashlight"', '"sv_prior"', text)

# Rule 5: "flashlight+step03" — a composite status string. Not caught by
# rules 1/2 because + isn't a word char and there's no _ boundary.
text = text.replace('"flashlight+step03"', '"sv_prior+step03"')

# Idempotency + safety checks
if text == original:
    print(f"[skip] no v2 patterns matched: {path}")
    sys.exit(10)

# Safety: we should NOT have touched these known-preserved strings.
preserved = [
    ("flashlight_v2",        "hardcoded stale directory name"),
    ("sv_flashlight_",       "legacy RDS filename fallback (C01g)"),
    ("_flashlight.rds",      "legacy RDS filename fallback (C01a)"),
    ("flashlight_loader_v2", "stale source file name"),
]
for needle, label in preserved:
    orig_count = original.count(needle)
    new_count  = text.count(needle)
    if orig_count != new_count:
        print(f"[rename-v2] ERROR: clobbered preserved string '{needle}' ({label}) "
              f"in {path} (was {orig_count}, now {new_count})", file=sys.stderr)
        sys.exit(3)

# Count substitutions (approximate — by counting flashlight references)
before_fl = len(re.findall(r'flashlight', original))
after_fl  = len(re.findall(r'flashlight', text))
n_hits = before_fl - after_fl

print(f"[rename-v2] {path}  ({n_hits} flashlight refs renamed)")
if dry_run:
    sys.exit(0)

shutil.copyfile(path, path + ".bk_chat17_scientific_v2")
with open(path, "w") as fh:
    fh.write(text)
PYEOF
  case $rc in
    0)  N_PATCHED=$((N_PATCHED + 1)) ;;
    10) N_SKIPPED=$((N_SKIPPED + 1)) ;;
    *)  echo "[rename-v2] FAIL: ${path} (python rc=$rc)" >&2; exit $rc ;;
  esac
done

echo ""
echo "[rename-v2] summary: ${N_PATCHED} patched, ${N_SKIPPED} already clean"

if [[ ${DRY_RUN} -eq 0 && ${N_PATCHED} -gt 0 ]]; then
  echo ""
  echo "[rename-v2] Verify:"
  echo "  # All R files still parse"
  echo "  for bk in \$(find \${BASE} -name '*.bk_chat17_scientific_v2'); do"
  echo "    orig=\${bk%.bk_chat17_scientific_v2}"
  echo "    [[ \$orig == *.R ]] && python3 \${BASE}/_rcheck.py \"\$orig\""
  echo "  done"
  echo ""
  echo "  # Only remaining flashlight refs should be:"
  echo "  # (a) 'sv_flashlight_' (legacy filename fallback)"
  echo "  # (b) 'flashlight_v2/...' (stale dir, not our scope)"
  echo "  # (c) comments/docstrings mentioning flashlight historically"
  echo "  grep -rn 'flashlight' --include='*.R' --include='*.py' --include='*.sh' \\"
  echo "    \${BASE} 2>/dev/null | grep -v _archive | grep -v bk_chat17 | \\"
  echo "    grep -v RENAMING | grep -v 'sv_flashlight_' | grep -v 'flashlight_v2'"
fi
