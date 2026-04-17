#!/usr/bin/env python3
"""
_schema_check.py — validate the three BK schemas.

Checks:
1. Parse as JSON.
2. Validate against JSON Schema draft-07 meta-schema.
3. For each keys_extracted entry, verify:
   - `key` is present and unique within the schema.
   - `from` is a dotted path rooted at a property declared in `properties`.
4. Report all findings; non-zero exit if any hard failure.

Run from work root.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path


SCHEMA_DIR = Path("registries/schemas/structured_block_schemas")

SCHEMAS_TO_CHECK = [
    "internal_dynamics.schema.json",
    "recombinant_map.schema.json",
    "internal_ancestry_composition.schema.json",
]


def resolve_path_in_properties(props: dict, path: str) -> tuple[bool, str]:
    """Walk a dotted path through a draft-07 'properties' tree.
    Returns (ok, note). ok=True means the first segment is declared in
    properties (we only validate the first segment; deeper segments are
    often nested-object keys whose sub-schema 'additionalProperties=true'
    or whose full shape is captured in description rather than properties).
    """
    parts = path.split(".")
    head = parts[0]
    if head not in props:
        return False, f"'{head}' not in top-level properties"
    # For single-segment paths we're done.
    if len(parts) == 1:
        return True, "ok"
    # Dotted path: try to walk into nested 'properties' where available.
    node = props[head]
    for p in parts[1:]:
        sub = node.get("properties") if isinstance(node, dict) else None
        if not isinstance(sub, dict):
            # Nested leaf without declared sub-properties — allowed (matches
            # e.g. gate_params.min_dco_bp when gate_params is declared with
            # properties; or structure_counts.* when additionalProperties
            # is used). Accept.
            return True, "ok (deep-leaf, not fully declared)"
        if p not in sub:
            return False, f"'{p}' not in properties of '{'.'.join(parts[:parts.index(p)])}'"
        node = sub[p]
    return True, "ok"


def main() -> int:
    errs: list[str] = []
    warns: list[str] = []

    # Optional: draft-07 meta-schema validation if jsonschema is available.
    try:
        import jsonschema  # type: ignore
        meta_ok = True
    except ImportError:
        meta_ok = False
        warns.append("jsonschema not installed; skipping draft-07 meta-schema validation")

    for fname in SCHEMAS_TO_CHECK:
        p = SCHEMA_DIR / fname
        print(f"\n=== {fname} ===")
        if not p.exists():
            errs.append(f"[{fname}] missing")
            continue

        # 1. JSON parse
        try:
            with open(p) as f:
                schema = json.load(f)
        except json.JSONDecodeError as e:
            errs.append(f"[{fname}] JSON parse error: {e}")
            continue
        print(f"  JSON OK ({p.stat().st_size} bytes)")

        # 2. draft-07 meta-schema
        if meta_ok:
            try:
                jsonschema.Draft7Validator.check_schema(schema)
                print("  draft-07 meta-schema OK")
            except jsonschema.SchemaError as e:
                errs.append(f"[{fname}] draft-07 SchemaError: {e.message}")

        # 3. keys_extracted checks
        kex = schema.get("keys_extracted")
        if kex is None:
            errs.append(f"[{fname}] no keys_extracted block")
            continue

        props = schema.get("properties", {})
        seen_keys: set[str] = set()
        seen_froms: dict[str, str] = {}
        for ke in kex:
            k = ke.get("key")
            frm = ke.get("from")
            if not k or not frm:
                errs.append(f"[{fname}] malformed ke: {ke}")
                continue
            if k in seen_keys:
                errs.append(f"[{fname}] duplicate key '{k}'")
            seen_keys.add(k)
            if frm in seen_froms:
                warns.append(f"[{fname}] two keys resolve from same path '{frm}': "
                             f"'{seen_froms[frm]}' and '{k}'")
            seen_froms[frm] = k
            ok, note = resolve_path_in_properties(props, frm)
            status = "OK" if ok else "MISSING"
            print(f"  {status:7s} {k:40s} <- {frm}  ({note})")
            if not ok:
                errs.append(f"[{fname}] key '{k}' from '{frm}' — {note}")

        print(f"  {len(kex)} keys_extracted entries, {len(seen_keys)} unique keys")

    print()
    print("=" * 60)
    if warns:
        print(f"WARNINGS ({len(warns)}):")
        for w in warns:
            print(f"  - {w}")
    if errs:
        print(f"ERRORS ({len(errs)}):")
        for e in errs:
            print(f"  - {e}")
        return 1
    print(f"ALL OK — {len(SCHEMAS_TO_CHECK)} schemas validated")
    return 0


if __name__ == "__main__":
    sys.exit(main())
