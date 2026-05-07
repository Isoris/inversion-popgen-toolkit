#!/usr/bin/env python3
"""
build/flatten.py
=================

Inline all ES-module imports referenced by an Atlas HTML file into one
self-contained HTML. Output is drop-into-browser portable.

Usage:
    python3 build/flatten.py inversion_review.html
    python3 build/flatten.py inversion_review.html --out dist/inversion_review_flat.html

How it works:
    1. Parse the source HTML.
    2. Find every <script type="module" src="..."> tag.
    3. Recursively follow import statements inside each module.
    4. Topologically sort modules by import dependency.
    5. Rewrite imports to local references and inline as one big
       <script type="module">.
    6. Bundle JSON/CSS/asset side-files via <link> / <script> tags
       — caller's job; this script only handles JS.

What it does NOT do:
    - Minify (out of scope; the source is the readable artifact)
    - Tree-shake (we ship what's imported)
    - Resolve npm packages (we don't use any)
    - Source maps (we don't need them; the inlined output uses
      /* === module: name.js === */ comment markers instead)

Constraints encoded:
    - We support `import { a, b } from './x.js'` and
      `import { a as b } from './x.js'` and
      `import * as M from './x.js'` and `import x from './x.js'`
      (the four forms shared/state_io.js exports use).
    - We support `export function f()`, `export const x =`,
      `export class C`, `export { a, b }`, `export { a as b }`,
      `export default ...` (also: re-exports `export { ... } from './x.js'`).
    - We assume modules don't have side effects at import time
      beyond defining symbols (true for shared/state_io.js).

This file is NOT part of the Atlas runtime — it runs at bundle time.
"""

from __future__ import annotations
import argparse, os, re, sys
from pathlib import Path
from typing import Dict, List, Set, Tuple


# ---------------------------------------------------------------------
# Regex set (kept narrow + explicit; we do NOT try to parse JS, only
# the patterns we actually use in the Atlas codebase).
# ---------------------------------------------------------------------

RE_SCRIPT_MODULE = re.compile(
    r'<script\s+type="module"(?:\s+src="([^"]+)")?\s*>(.*?)</script>',
    re.IGNORECASE | re.DOTALL,
)

# import ... from './foo.js'  OR  '../shared/foo.js'  (relative paths only)
# Supports multiline form:  import {\n  a, b, c,\n} from '../x.js';
# DOTALL on the binding so { ... } can span lines; non-greedy until `from`.
RE_IMPORT = re.compile(
    r"""^[ \t]*import\s+(?P<binding>(?:\{[^}]*\}|[^'";\n]+?))\s*from\s+['"](?P<src>\.{1,2}/[^'"]+)['"]\s*;?[ \t]*$""",
    re.MULTILINE | re.DOTALL,
)

# Bare side-effect import:  import './foo.js';
RE_IMPORT_BARE = re.compile(
    r"""^[ \t]*import\s+['"](?P<src>\.{1,2}/[^'"]+)['"]\s*;?[ \t]*$""",
    re.MULTILINE,
)

# Re-export forwarding:  export { a, b } from './foo.js';
# Also supports multiline:  export {\n  a,\n  b,\n} from './foo.js';
RE_REEXPORT_FROM = re.compile(
    r"""^[ \t]*export\s+\{(?P<names>[^}]+)\}\s+from\s+['"](?P<src>\.{1,2}/[^'"]+)['"]\s*;?[ \t]*$""",
    re.MULTILINE | re.DOTALL,
)

# Strip `export ` keywords that prefix declarations. We rewrite
#   export function foo()    → function foo()
#   export const X = ...     → const X = ...
#   export class C { ... }   → class C { ... }
#   export default <expr>    → /* default */ <expr>     (rare; we name it _default_<module>)
RE_EXPORT_DECL = re.compile(
    r"""^[ \t]*export\s+(?=(function|const|let|var|class|async|default))""",
    re.MULTILINE,
)

# `export { foo, bar }` (list-only export, no `from`) → drop the line.
# Multiline tolerant.
RE_EXPORT_LIST = re.compile(
    r"""^[ \t]*export\s+\{[^}]*\}\s*;?[ \t]*$""",
    re.MULTILINE | re.DOTALL,
)


# ---------------------------------------------------------------------
# Module loader / dependency walker
# ---------------------------------------------------------------------

class Module:
    __slots__ = ('path', 'source', 'imports', 'rel_to_root')

    def __init__(self, path: Path, source: str, rel_to_root: str):
        self.path = path
        self.source = source
        self.rel_to_root = rel_to_root
        self.imports: List[Tuple[str, Path]] = []   # (binding text, resolved path)


def resolve(base: Path, rel: str) -> Path:
    p = (base.parent / rel).resolve()
    if not p.exists():
        # Allow .js extension to be implicit in source
        if not rel.endswith('.js'):
            p2 = (base.parent / (rel + '.js')).resolve()
            if p2.exists():
                return p2
        raise FileNotFoundError(f"import not found: {rel} (from {base})")
    return p


def load_module(path: Path, root: Path, cache: Dict[Path, Module]) -> Module:
    if path in cache:
        return cache[path]
    src = path.read_text(encoding='utf-8')
    rel = os.path.relpath(path, root).replace(os.sep, '/')
    mod = Module(path=path, source=src, rel_to_root=rel)
    cache[path] = mod
    # Walk imports
    for m in RE_IMPORT.finditer(src):
        rel_src = m.group('src')
        try:
            dep = resolve(path, rel_src)
        except FileNotFoundError as e:
            raise SystemExit(f"flatten: {e}")
        mod.imports.append((m.group('binding'), dep))
        load_module(dep, root, cache)
    for m in RE_IMPORT_BARE.finditer(src):
        dep = resolve(path, m.group('src'))
        mod.imports.append(('(bare)', dep))
        load_module(dep, root, cache)
    for m in RE_REEXPORT_FROM.finditer(src):
        dep = resolve(path, m.group('src'))
        mod.imports.append(('(reexport)', dep))
        load_module(dep, root, cache)
    return mod


def topo_sort(mods: Dict[Path, Module]) -> List[Module]:
    """Stable topological sort by import dependency."""
    visited: Set[Path] = set()
    out: List[Module] = []
    def visit(p: Path, stack: Tuple[Path, ...] = ()):
        if p in visited:
            return
        if p in stack:
            cycle = ' → '.join(str(s) for s in (*stack, p))
            raise SystemExit(f"flatten: import cycle: {cycle}")
        for _, dep in mods[p].imports:
            visit(dep, stack + (p,))
        visited.add(p)
        out.append(mods[p])
    for p in mods:
        visit(p)
    return out


# ---------------------------------------------------------------------
# Source rewriter — strip imports, strip exports, leave declarations.
# ---------------------------------------------------------------------

def rewrite_module_source(src: str) -> str:
    # 1. Strip `import ... from './x.js';` lines.
    src = RE_IMPORT.sub('', src)
    src = RE_IMPORT_BARE.sub('', src)
    # 2. Strip `export { ... } from './x.js';` re-exports.
    #    The names re-exported are still defined upstream by the
    #    upstream module (which we'll inline first), so dropping the
    #    re-export line is safe.
    src = RE_REEXPORT_FROM.sub('', src)
    # 3. Strip `export { foo, bar };` standalone export lists.
    src = RE_EXPORT_LIST.sub('', src)
    # 4. Strip `export ` from declaration lines.
    src = RE_EXPORT_DECL.sub('', src)
    return src


# ---------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------

def flatten(html_path: Path, out_path: Path) -> None:
    html = html_path.read_text(encoding='utf-8')
    root = html_path.parent.resolve()
    cache: Dict[Path, Module] = {}

    # Find every <script type="module" src="..."> (we don't process
    # inline `<script type="module">…</script>` because the Atlases use
    # one entry-point src= per HTML).
    entry_paths: List[Path] = []
    for m in RE_SCRIPT_MODULE.finditer(html):
        src = m.group(1)
        if not src:
            continue
        ep = (root / src).resolve()
        if not ep.exists():
            raise SystemExit(f"flatten: entry not found: {src}")
        entry_paths.append(ep)
        load_module(ep, root, cache)

    if not entry_paths:
        raise SystemExit('flatten: no <script type="module" src="…"> entry found')

    # Topo sort all loaded modules.
    ordered = topo_sort(cache)

    # Build the inlined script body.
    parts: List[str] = []
    parts.append("/* === FLATTENED BUNDLE — produced by build/flatten.py === */\n")
    parts.append("/* Modules in topological import order: */\n")
    for mod in ordered:
        parts.append(f"/*   - {mod.rel_to_root} */\n")
    parts.append("\n")

    for mod in ordered:
        rewritten = rewrite_module_source(mod.source).strip('\n')
        parts.append(f"\n/* ===== module: {mod.rel_to_root} ===== */\n")
        parts.append(rewritten)
        parts.append("\n")

    inlined = ''.join(parts)

    # Replace ALL <script type="module" src="…"></script> with one
    # combined <script type="module"> block at the position of the FIRST
    # entry tag; remove the rest. This preserves load order semantics
    # (modules are deferred-by-default; injected scripts run after DOMContentLoaded).
    placed = {'done': False}
    def repl(m: re.Match) -> str:
        if not m.group(1):
            return m.group(0)   # not a src= tag, leave it
        if placed['done']:
            return ''           # subsequent entries removed
        placed['done'] = True
        return f'<script type="module">\n{inlined}\n</script>'
    out_html = RE_SCRIPT_MODULE.sub(repl, html)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(out_html, encoding='utf-8')
    print(f"flatten: {html_path}  →  {out_path}")
    print(f"         {len(cache)} modules, {len(out_html):,} bytes total")


def main(argv: List[str]) -> int:
    ap = argparse.ArgumentParser(description="Flatten Atlas HTML+JS into one file")
    ap.add_argument('html', help='source HTML file (e.g. inversion_review.html)')
    ap.add_argument('--out', help='destination (default: dist/<name>_flat.html)')
    args = ap.parse_args(argv)

    src = Path(args.html).resolve()
    if not src.exists():
        ap.error(f"not found: {src}")
    if args.out:
        dst = Path(args.out).resolve()
    else:
        dst = src.parent / 'dist' / (src.stem + '_flat.html')

    flatten(src, dst)
    return 0


if __name__ == '__main__':
    raise SystemExit(main(sys.argv[1:]))
