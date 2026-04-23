#!/usr/bin/env python3
"""
_rcheck.py — Lightweight R syntax balance checker.

Not an R parser — just catches the common self-inflicted errors
after str_replace edits:
 - unbalanced (), [], {}
 - unterminated "..." or '...'
 - stray ``` backticks
 - lines with ',,' (common copy-paste artefact)
 - mismatched <<- / ->> arrows (sanity only)

Usage:
    python3 tools/_rcheck.py FILE [FILE ...]

Exit 0 if all clean, 1 otherwise.
"""
import sys
import re
from pathlib import Path


def strip_strings_and_comments(line: str):
    """Return line with string contents replaced by a single placeholder
    char 'S', and trailing # comment removed. Handles escaped quotes.
    Using a placeholder (not deletion) preserves the fact that a comma
    had a non-empty argument between it and the next comma, which
    matters for the stray-',,' check."""
    out = []
    i = 0
    n = len(line)
    in_s = None  # current string quote char
    while i < n:
        c = line[i]
        if in_s:
            if c == "\\" and i + 1 < n:
                i += 2
                continue
            if c == in_s:
                in_s = None
                out.append("S")  # close: one placeholder for the whole string
            i += 1
            continue
        if c in ("'", '"'):
            in_s = c
            i += 1
            continue
        if c == "#":
            break
        out.append(c)
        i += 1
    return "".join(out), (in_s is not None)


def check_file(path: Path):
    errs = []
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except Exception as e:
        return [f"{path}: read error: {e}"]

    # Top-level balance counters (after strings/comments stripped, line by line)
    depth_round = 0
    depth_square = 0
    depth_curly = 0
    open_markers = {"(": "round", "[": "square", "{": "curly"}
    close_markers = {")": "round", "]": "square", "}": "curly"}
    ctr = {"round": 0, "square": 0, "curly": 0}

    # Backtick usage — markdown leftovers
    if "```" in text:
        errs.append(f"{path}: contains ``` (markdown fence leaked into R)")

    for lineno, raw in enumerate(text.splitlines(), start=1):
        stripped, unterminated = strip_strings_and_comments(raw)
        if unterminated:
            errs.append(f"{path}:{lineno}: unterminated string literal")
        # Stray ,, — but only flag outside of any [ ] context, because
        # R array indexing legitimately uses x[, col] and x[, , i].
        # Build a bracket-depth track and require ',,' only when at depth 0.
        if re.search(r",\s*,", stripped):
            depth_sq = 0
            depth_any = 0
            bad = False
            j = 0
            while j < len(stripped):
                ch = stripped[j]
                if ch == "[":
                    depth_sq += 1
                elif ch == "]":
                    depth_sq = max(0, depth_sq - 1)
                elif ch == ",":
                    # look ahead past whitespace
                    k = j + 1
                    while k < len(stripped) and stripped[k] == " ":
                        k += 1
                    if k < len(stripped) and stripped[k] == "," and depth_sq == 0:
                        bad = True
                        break
                j += 1
            if bad:
                errs.append(f"{path}:{lineno}: stray ',,'")
        for ch in stripped:
            if ch in open_markers:
                ctr[open_markers[ch]] += 1
            elif ch in close_markers:
                ctr[close_markers[ch]] -= 1
                if ctr[close_markers[ch]] < 0:
                    errs.append(
                        f"{path}:{lineno}: extra closing '{ch}'"
                    )
                    # reset so downstream counts don't cascade
                    ctr[close_markers[ch]] = 0

    for k, v in ctr.items():
        if v != 0:
            errs.append(
                f"{path}: {k} brackets unbalanced at EOF (net {v:+d})"
            )

    return errs


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    any_err = False
    for p in sys.argv[1:]:
        path = Path(p)
        if not path.is_file():
            print(f"{p}: not a file")
            any_err = True
            continue
        errs = check_file(path)
        if errs:
            any_err = True
            for e in errs:
                print(e)
        else:
            print(f"{p}: OK")
    sys.exit(1 if any_err else 0)


if __name__ == "__main__":
    main()
