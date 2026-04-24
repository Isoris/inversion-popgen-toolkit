#!/usr/bin/env python3
"""
scan_script_contracts.py — walk scripts + schemas, cross-check contracts.

See docs/SCRIPT_CONTRACT.md for the contract header spec.

USAGE
    cd <repo root>
    python3 tools/scan_script_contracts.py
    python3 tools/scan_script_contracts.py --strict   # exit nonzero on drift

WHAT IT CHECKS
    1. Every script with a REGISTRY_CONTRACT block parses cleanly.
    2. Each BLOCKS_WRITTEN entry's schema path exists.
    3. Each keys: list is a subset of its schema's keys_extracted keys
       (for status=WIRED only; other statuses are allowed to have empty
       or partial keys: lists by design).
    4. Each schema in registries/schemas/structured_block_schemas/ has at
       least one script claiming BLOCKS_WRITTEN on its block_type.
    5. No key is claimed by more than one WIRED contract (double-writers
       are usually a mistake, unless both writers agree on semantics —
       the scanner emits a warning, not a hard failure).

EXIT CODES
    0   All checks passed (or --strict not set and only warnings).
    1   Hard drift: contract names a key that isn't in its schema's
        keys_extracted, OR BLOCKS_WRITTEN entry points to a missing
        schema file. With --strict, any warning also triggers exit 1.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable

REPO_ROOT = Path(__file__).resolve().parent.parent

# Where to look
SCRIPT_DIRS = [
    REPO_ROOT / "inversion_modules",
    REPO_ROOT / "unified_ancestry",
    REPO_ROOT / "utils",
]
SCHEMA_DIR = REPO_ROOT / "registries" / "schemas" / "structured_block_schemas"

# Skip archived / superseded content
SKIP_DIR_FRAGMENTS = ("/_archive", "/_archive_superseded", "/superseded",
                      "/.git/", "/node_modules/")

SCRIPT_SUFFIXES = {".R", ".py", ".sh"}


# =============================================================================
# Contract block parsing
# =============================================================================

CONTRACT_RE = re.compile(r"REGISTRY_CONTRACT\s*$", re.MULTILINE)
BLOCK_ENTRY_RE = re.compile(
    r"^\s*[#*]*\s*-\s+(?P<block>[\w]+)\s*:\s*(?P<schema>[^\s]+)\s*$"
)
# Continuation lines under a BLOCKS_WRITTEN entry
FIELD_RE = re.compile(r"^\s*[#*]*\s*(keys|keys_na|status|note):\s*(?P<val>.*)$")

# Accept these status markers
STATUS_OK = {"WIRED", "PRODUCES_BUT_NOT_WIRED"}


def iter_script_files() -> Iterable[Path]:
    for root in SCRIPT_DIRS:
        if not root.exists():
            continue
        for p in root.rglob("*"):
            if not p.is_file():
                continue
            if p.suffix not in SCRIPT_SUFFIXES:
                continue
            s = str(p)
            if any(frag in s for frag in SKIP_DIR_FRAGMENTS):
                continue
            yield p


def read_head(path: Path, n_lines: int = 300) -> list[str]:
    try:
        with path.open("r", encoding="utf-8", errors="replace") as f:
            return [next(f) for _ in range(n_lines)]
    except StopIteration:
        try:
            return path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)
        except Exception:
            return []
    except Exception:
        return []


def find_contract_start(lines: list[str]) -> int | None:
    for i, line in enumerate(lines):
        if CONTRACT_RE.search(line):
            return i
    return None


def parse_contract(lines: list[str], start: int, path: Path) -> dict | None:
    """Parse a REGISTRY_CONTRACT block starting at line index `start`.

    Returns a dict:
      {
        "path": path,
        "blocks": [ { "block_type", "schema_path", "keys": [...],
                       "keys_na": [...], "status", "note" }, ... ],
        "keys_in": ["none"] or list,
        "parse_errors": [ ... ]
      }

    Stops at the first blank line that follows a non-blank in-contract line
    (so the block ends wherever the script's header block naturally does).
    Also stops at divider lines of ==== or ---- of length >= 30.
    """
    errors: list[str] = []
    blocks: list[dict] = []
    keys_in: list[str] = []
    mode = "HEADER"         # initial state
    current: dict | None = None

    # Strip leading comment markers from each line
    def strip_comment(s: str) -> str:
        s = s.rstrip("\n")
        # common comment prefixes: '#', '# ', '//', '/*', '*', '"""'
        s = re.sub(r"^\s*(#|//|/\*|\*)\s?", "", s)
        return s

    # Track the indent of the REGISTRY_CONTRACT keyword for sub-field parsing
    i = start + 1
    while i < len(lines):
        raw = lines[i]
        # Divider line ends the block
        if re.match(r"^\s*[#*/]*\s*[=\-]{20,}\s*$", raw):
            break
        stripped = strip_comment(raw)
        s = stripped.strip()
        if s == "":
            # blank line within the contract — tolerated (for readability),
            # but two in a row ends the block
            if mode == "DONE_EMPTY":
                break
            mode = "DONE_EMPTY"
            i += 1
            continue
        mode = "ACTIVE"

        # Major section markers
        if s.startswith("BLOCKS_WRITTEN"):
            # Support "BLOCKS_WRITTEN: none" as an explicit declaration
            # that this is a compute engine or utility with no registry
            # writes — the file is allowed to have a contract block so
            # readers can confirm it was deliberately scoped.
            rest = s[len("BLOCKS_WRITTEN"):].lstrip(":").strip()
            if rest.lower() in {"none", "(none)", "()", "— compute engine"}:
                # Synthesize a stub entry so the scanner doesn't report
                # "no BLOCKS_WRITTEN entries parsed"
                current = {
                    "block_type": "_NONE",
                    "schema_path": "",
                    "keys": [],
                    "keys_na": [],
                    "status": "NONE_DECLARED",
                    "note": "explicitly declared no registry writes",
                }
                blocks.append(current)
                current = None
            i += 1
            continue
        if s.startswith("KEYS_IN"):
            # Inline or next-line value
            rest = s[len("KEYS_IN"):].lstrip(":").strip()
            if rest:
                keys_in = [k.strip() for k in rest.split(",") if k.strip()]
            i += 1
            # swallow continuation lines starting with "-"
            while i < len(lines):
                nxt_raw = lines[i]
                if re.match(r"^\s*[#*/]*\s*[=\-]{20,}\s*$", nxt_raw):
                    break
                nxt = strip_comment(nxt_raw).strip()
                if not nxt:
                    break
                m_dash = re.match(r"^-\s+(.+)$", nxt)
                if m_dash:
                    keys_in.append(m_dash.group(1).strip())
                    i += 1
                    continue
                # new section — stop swallowing
                break
            continue

        # Entry start line: "- <block>: <schema>"
        m = BLOCK_ENTRY_RE.match(stripped)
        if m:
            if current is not None:
                blocks.append(current)
            current = {
                "block_type": m.group("block"),
                "schema_path": m.group("schema"),
                "keys": [],
                "keys_na": [],
                "status": None,
                "note": None,
            }
            i += 1
            continue

        # Sub-field lines (keys / keys_na / status / note)
        mf = FIELD_RE.match(stripped)
        if mf:
            if current is None:
                errors.append(
                    f"{path}:{i+1}: sub-field '{mf.group(1)}' with no open entry"
                )
                i += 1
                continue
            field = mf.group(1)
            val = mf.group(2).strip()
            # Handle multi-line keys: / keys_na: wrapped values
            if field in ("keys", "keys_na"):
                # Collect continuation lines that start with whitespace only
                collected = [val]
                j = i + 1
                while j < len(lines):
                    cont_raw = lines[j]
                    if re.match(r"^\s*[#*/]*\s*[=\-]{20,}\s*$", cont_raw):
                        break
                    cont = strip_comment(cont_raw).rstrip()
                    if cont == "" or cont.strip() == "":
                        break
                    # if this line looks like a new field/entry, stop
                    if BLOCK_ENTRY_RE.match(cont) or FIELD_RE.match(cont):
                        break
                    if cont.strip().startswith(("BLOCKS_WRITTEN", "KEYS_IN")):
                        break
                    collected.append(cont.strip())
                    j += 1
                combined = " ".join(collected).strip()
                # Strip surrounding parens like "(none)"
                if combined.lower() in {"(none)", "none", "()"}:
                    keys_list = []
                else:
                    keys_list = [k.strip() for k in combined.split(",") if k.strip()]
                current[field] = keys_list
                i = j
                continue
            elif field == "status":
                current["status"] = val
                i += 1
                continue
            elif field == "note":
                # capture note — can continue on subsequent lines with
                # consistent indent, but we'll settle for the single-line
                # value for now (the scanner doesn't validate note content)
                current["note"] = val
                i += 1
                continue

        # Line that isn't recognized but doesn't end the block: skip
        i += 1

    if current is not None:
        blocks.append(current)

    if not blocks:
        errors.append(f"{path}: REGISTRY_CONTRACT block found but no BLOCKS_WRITTEN entries parsed")

    return {
        "path": path,
        "blocks": blocks,
        "keys_in": keys_in,
        "parse_errors": errors,
    }


# =============================================================================
# Schema loading
# =============================================================================

def load_schema(schema_path: Path) -> tuple[dict | None, str | None]:
    if not schema_path.exists():
        return None, f"schema file not found: {schema_path}"
    try:
        return json.loads(schema_path.read_text()), None
    except json.JSONDecodeError as e:
        return None, f"schema JSON invalid: {e}"


def schema_extracted_keys(schema: dict) -> list[str]:
    """Return declared extracted keys, with {side} templates expanded.

    Templated keys like 'q3_{side}_bp' are expanded into the set of
    concrete keys the writer will emit: 5prime, 3prime (canonical) plus
    left, right (legacy dual-write). A wired contract may name either
    the canonical forms, the legacy forms, or both.

    Non-templated keys pass through unchanged.
    """
    ke = schema.get("keys_extracted", []) or []
    out: list[str] = []
    for entry in ke:
        if not isinstance(entry, dict):
            continue
        k = entry.get("key")
        if not isinstance(k, str):
            continue
        if "{side}" in k:
            for sub in ("5prime", "3prime", "left", "right"):
                out.append(k.replace("{side}", sub))
        else:
            out.append(k)
    return out


def enumerate_schemas() -> dict[str, Path]:
    """block_type -> schema path, for the canonical schema dir.

    Prefers explicit `block_type`; falls back to `title`, then to the
    filename without `.schema.json`. The fallbacks catch schemas that
    were authored before `block_type` was standardized (a handful of
    chat-13 regime-family schemas still use `title` instead).
    """
    m: dict[str, Path] = {}
    if not SCHEMA_DIR.exists():
        return m
    for p in sorted(SCHEMA_DIR.glob("*.schema.json")):
        try:
            d = json.loads(p.read_text())
        except Exception:
            continue
        bt = d.get("block_type") or d.get("title")
        if not isinstance(bt, str):
            # Derive from filename: "regime_segments.schema.json" -> "regime_segments"
            bt = p.name[:-len(".schema.json")]
        m.setdefault(bt, p)
    return m


# =============================================================================
# Main
# =============================================================================

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1] if __doc__ else "")
    ap.add_argument("--strict", action="store_true",
                    help="exit 1 even on warnings (missing coverage, "
                         "duplicated keys)")
    ap.add_argument("--verbose", "-v", action="store_true")
    args = ap.parse_args()

    schemas = enumerate_schemas()

    parsed: list[dict] = []
    parse_errors: list[str] = []
    for path in iter_script_files():
        lines = read_head(path)
        if not lines:
            continue
        start = find_contract_start(lines)
        if start is None:
            continue
        res = parse_contract(lines, start, path)
        parsed.append(res)
        parse_errors.extend(res["parse_errors"])

    # -------------------------------------------------------------------------
    # Check 1: every BLOCKS_WRITTEN entry's schema_path exists + keys subset
    # -------------------------------------------------------------------------
    hard_errors: list[str] = []
    warnings: list[str] = []

    # block_type -> list of (path, entry, schema_keys_set)
    by_block: dict[str, list[tuple[Path, dict, set]]] = defaultdict(list)

    for res in parsed:
        path = res["path"]
        for entry in res["blocks"]:
            # Explicitly declared "no blocks written" — engines / utilities
            if entry.get("block_type") == "_NONE":
                continue
            schema_path = REPO_ROOT / entry["schema_path"]
            schema, err = load_schema(schema_path)
            if err:
                hard_errors.append(f"{path}: block {entry['block_type']}: {err}")
                continue
            declared_keys = set(schema_extracted_keys(schema))
            entry_keys = set(entry.get("keys") or [])
            status = (entry.get("status") or "").upper()

            # Drift check: contract names a key not in schema
            if status == "WIRED":
                extra = entry_keys - declared_keys
                if extra:
                    hard_errors.append(
                        f"{path}: block {entry['block_type']} (WIRED): "
                        f"contract names keys not in schema.keys_extracted: "
                        f"{sorted(extra)}"
                    )

            # Schema block_type mismatch
            schema_bt = schema.get("block_type")
            if schema_bt and schema_bt != entry["block_type"]:
                warnings.append(
                    f"{path}: block {entry['block_type']} points to schema "
                    f"{schema_path} whose block_type is '{schema_bt}'"
                )

            # Status check
            if status and status != "WIRED" and status != "PRODUCES_BUT_NOT_WIRED" \
                    and not status.startswith("BLOCKED_ON_"):
                warnings.append(
                    f"{path}: block {entry['block_type']}: unknown status "
                    f"'{entry['status']}' (expected WIRED, PRODUCES_BUT_NOT_WIRED, "
                    f"or BLOCKED_ON_*)"
                )

            by_block[entry["block_type"]].append((path, entry, declared_keys))

    # -------------------------------------------------------------------------
    # Check 2: every schema in SCHEMA_DIR has a claimant
    # -------------------------------------------------------------------------
    uncovered: list[str] = []
    for bt, sp in sorted(schemas.items()):
        if bt not in by_block:
            uncovered.append(bt)

    # -------------------------------------------------------------------------
    # Check 3: a WIRED contract "owns" every key it declares; check for
    #          duplicate WIRED writers claiming the same key (warn).
    # -------------------------------------------------------------------------
    key_claimants: dict[tuple[str, str], list[Path]] = defaultdict(list)
    for bt, entries in by_block.items():
        for path, entry, _schema_keys in entries:
            status = (entry.get("status") or "").upper()
            if status != "WIRED":
                continue
            for k in entry.get("keys") or []:
                key_claimants[(bt, k)].append(path)
    for (bt, k), paths in sorted(key_claimants.items()):
        if len(paths) > 1:
            warnings.append(
                f"key {k} in block {bt} is claimed WIRED by multiple "
                f"scripts: {[str(p) for p in paths]}"
            )

    # -------------------------------------------------------------------------
    # Report
    # -------------------------------------------------------------------------
    print(f"[scan] scripts with REGISTRY_CONTRACT: {len(parsed)}")
    print(f"[scan] schemas in registry: {len(schemas)}")

    # group by status
    by_status: dict[str, int] = defaultdict(int)
    for res in parsed:
        for entry in res["blocks"]:
            by_status[(entry.get("status") or "UNSPECIFIED").upper()] += 1
    print(f"[scan] BLOCKS_WRITTEN entries by status:")
    for s, n in sorted(by_status.items()):
        print(f"         {s:30s} {n:4d}")

    if args.verbose:
        print(f"\n[scan] schemas covered by ≥1 contract:")
        for bt in sorted(schemas):
            if bt in by_block:
                who = ", ".join(sorted({str(p.name) for p, _, _ in by_block[bt]}))
                print(f"         {bt:40s} ← {who}")

    if parse_errors:
        print(f"\n[scan] PARSE ERRORS ({len(parse_errors)}):")
        for m in parse_errors:
            print(f"         {m}")

    if hard_errors:
        print(f"\n[scan] HARD DRIFT ({len(hard_errors)}):")
        for m in hard_errors:
            print(f"         ERR {m}")

    if warnings:
        print(f"\n[scan] WARNINGS ({len(warnings)}):")
        for m in warnings:
            print(f"         WRN {m}")

    if uncovered:
        print(f"\n[scan] UNCOVERED SCHEMAS ({len(uncovered)}) — no script "
              f"declares BLOCKS_WRITTEN for these block_types:")
        for bt in sorted(uncovered):
            print(f"         {bt}  ({schemas[bt].name})")

    # Exit status
    if hard_errors:
        return 1
    if args.strict and (warnings or uncovered or parse_errors):
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
