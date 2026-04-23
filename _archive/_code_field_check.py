#!/usr/bin/env python3
"""
_code_field_check.py — validate that every keys_extracted `from` path
can be satisfied by a field that the source code actually writes.

For each of the three BK schemas, we read the corresponding source
script and extract the list of top-level fields written in block_data.
Then we check that every `from` path's head segment is in that set.

This is static analysis (regex-over-source), not runtime verification.
It will miss conditional fields written in some branches but not others
— those are reported as warnings, not errors, because the schema allows
them to be absent.
"""
from __future__ import annotations
import json
import re
import sys
from pathlib import Path
from collections import defaultdict


SCHEMA_DIR = Path("registries/schemas/structured_block_schemas")

# Each schema -> (source file, function/context hint for grepping the right block)
TARGETS = {
    "internal_dynamics.schema.json": {
        "sources": ["inversion_modules/phase_4_postprocessing/4b_group_proposal/STEP_C01i_decompose.R"],
        "block_type": "internal_dynamics",
    },
    "recombinant_map.schema.json": {
        "sources": ["inversion_modules/phase_4_postprocessing/4b_group_proposal/STEP_C01i_b_multi_recomb.R"],
        "block_type": "recombinant_map",
    },
    "internal_ancestry_composition.schema.json": {
        "sources": ["inversion_modules/phase_4_postprocessing/4b_group_proposal/STEP_C01i_c_nested_composition.py",
                    "inversion_modules/phase_4_postprocessing/4b_group_proposal/nested_composition_core.py"],
        "block_type": "internal_ancestry_composition",
    },
}


# Patterns to find field names being written inside a block_data / dict construction.
# - R:  block_data <- list( candidate_id = cid, chrom = chr, ... )
#       Also: write_block(... data = list(parent_id = ..., ...))
# - Py: data = { "parent_id": cid, ... }
#       and  return { "parent_id": pid, ... }
# We intentionally over-approximate by collecting any "symbol =" or "string_key:"
# pattern in regions near block_data / data dict starts.
R_KEY_PAT = re.compile(r"^\s*([a-zA-Z_][a-zA-Z_0-9]*)\s*=", re.MULTILINE)
PY_KEY_PAT = re.compile(r"""^\s*["']([a-zA-Z_][a-zA-Z_0-9]*)["']\s*:""", re.MULTILINE)


def extract_fields_r(text: str) -> set[str]:
    """R: find all 'name =' lines inside block_data <- list(...) bodies.

    KNOWN LIMITATION: the comma-splitter under-reports fields when a
    value is a multi-line if/else expression, because a top-level comma
    inside the body CAN appear inside such an expression and prematurely
    close a segment. On STEP_C01i_decompose.R this returns ~15/27 fields.
    As a second pass we use a fallback heuristic: any line matching
    ^<whitespace><ident><whitespace>=<non-=> inside the block_data
    body is also admitted. This is permissive (could false-positive on
    local-variable assignments that happen to sit inside the list() body),
    but internal_dynamics block_data contains no such locals — everything
    between `block_data <- list(` and its matching `)` is a field
    assignment. For the other R script (multi_recomb.R) both paths give
    the same answer. Cross-check the flagged set with the schema diff
    before ruling a field missing.
    """
    fields = set()
    for m in re.finditer(r"block_data\s*<-\s*list\s*\(", text):
        start = m.end()
        depth = 1
        i = start
        in_str = None  # None, '"', or "'"
        while i < len(text) and depth > 0:
            c = text[i]
            if in_str:
                if c == "\\" and i + 1 < len(text):
                    i += 2
                    continue
                if c == in_str:
                    in_str = None
            else:
                if c in ('"', "'"):
                    in_str = c
                elif c == "(":
                    depth += 1
                elif c == ")":
                    depth -= 1
            i += 1
        body = text[start:i - 1]

        # Split body into top-level comma-separated segments
        segs = []
        depth2 = 0
        in_str2 = None
        buf = []
        for c in body:
            if in_str2:
                buf.append(c)
                if c == "\\":
                    continue
                if c == in_str2:
                    in_str2 = None
                continue
            if c in ('"', "'"):
                in_str2 = c
                buf.append(c)
                continue
            if c == "(":
                depth2 += 1
                buf.append(c)
            elif c == ")":
                depth2 -= 1
                buf.append(c)
            elif c == "," and depth2 == 0:
                segs.append("".join(buf))
                buf = []
            else:
                buf.append(c)
        if buf:
            segs.append("".join(buf))

        for seg in segs:
            # First non-blank line of the segment should start with `name =`
            stripped = seg.lstrip("\n\r ")
            mm = re.match(r"^([a-zA-Z_][a-zA-Z_0-9.]*)\s*=[^=]", stripped)
            if mm:
                fields.add(mm.group(1))

        # Permissive fallback: admit any `^<ws><ident> =<not-=>` line inside
        # the body. Handles multi-line if/else values where the strict
        # comma-splitter misfires.
        for line in body.splitlines():
            mm = re.match(r"^\s+([a-zA-Z_][a-zA-Z_0-9.]*)\s*=[^=]", line)
            if mm:
                fields.add(mm.group(1))
    return fields


def extract_fields_py(text: str) -> set[str]:
    """Python: find string-keyed dict bodies immediately following 'return {' or 'data = {'."""
    fields = set()
    for pattern in [r"return\s*\{", r"data\s*=\s*\{"]:
        for m in re.finditer(pattern, text):
            start = m.end()
            depth = 1
            i = start
            while i < len(text) and depth > 0:
                if text[i] == "{":
                    depth += 1
                elif text[i] == "}":
                    depth -= 1
                i += 1
            body = text[start:i - 1]
            # Top-level string keys only (depth 0 inside body)
            depth2 = 0
            buf = []
            segs = []
            for ch in body:
                if ch in "{[(":
                    depth2 += 1
                    buf.append(ch)
                elif ch in "}])":
                    depth2 -= 1
                    buf.append(ch)
                elif ch == "," and depth2 == 0:
                    segs.append("".join(buf))
                    buf = []
                else:
                    buf.append(ch)
            if buf:
                segs.append("".join(buf))
            for seg in segs:
                mm = re.match(r'^\s*["\']([a-zA-Z_][a-zA-Z_0-9]*)["\']\s*:', seg)
                if mm:
                    fields.add(mm.group(1))
    return fields


def gather_fields(sources: list[str]) -> set[str]:
    fields: set[str] = set()
    for src in sources:
        p = Path(src)
        if not p.exists():
            print(f"  [warn] source not found: {src}", file=sys.stderr)
            continue
        text = p.read_text()
        if src.endswith(".R"):
            fields |= extract_fields_r(text)
        elif src.endswith(".py"):
            fields |= extract_fields_py(text)
    return fields


def main() -> int:
    errs: list[str] = []
    warns: list[str] = []
    for fname, cfg in TARGETS.items():
        print(f"\n=== {fname} ===")
        schema_path = SCHEMA_DIR / fname
        schema = json.loads(schema_path.read_text())
        kex = schema.get("keys_extracted", [])

        written = gather_fields(cfg["sources"])
        print(f"  fields written by {cfg['sources'][0].split('/')[-1]}"
              f"{'+' + str(len(cfg['sources']) - 1) + ' more' if len(cfg['sources']) > 1 else ''}: "
              f"{sorted(written)}")
        print()

        for ke in kex:
            head = ke["from"].split(".")[0]
            key = ke["key"]
            status = "OK   " if head in written else "MISS "
            print(f"  {status} {key:40s} <- {ke['from']}")
            if head not in written:
                errs.append(f"[{fname}] '{key}' from '{ke['from']}' — "
                            f"head '{head}' NOT written by source")

    print()
    print("=" * 60)
    if errs:
        print(f"ERRORS ({len(errs)}):")
        for e in errs:
            print(f"  - {e}")
        return 1
    print("ALL OK — every keys_extracted head segment is a field the source code writes")
    return 0


if __name__ == "__main__":
    sys.exit(main())
