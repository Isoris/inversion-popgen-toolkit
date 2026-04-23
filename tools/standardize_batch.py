#!/usr/bin/env python3
"""
standardize_batch.py — apply the standardized header/echo format to many files

THIS IS THE SECOND PHASE. Do NOT run it until you've tested the super prompt
on ~5 sample files by hand and are happy with the output quality.

Workflow:

  # Phase 1: manual iteration (done separately in a chat interface)
  # Phase 2 (this script): once the prompt is frozen, apply it in bulk

  python3 standardize_batch.py \\
    --file_list files_to_standardize.txt \\
    --prompt STANDARDIZE_SCRIPT.md \\
    --api_key $ANTHROPIC_API_KEY \\
    --model claude-sonnet-4-6 \\
    --staging_dir ./staged_updates \\
    --dry_run

  # If dry run looks good:
  python3 standardize_batch.py ... --apply

Safety features:
  - Every file gets committed to staging_dir first
  - Diff is computed against original; any non-comment/non-echo changes
    are flagged and NOT auto-applied
  - --apply is required to actually modify files in place
  - All changes are logged to standardize.log

Dependencies:
  pip install anthropic
"""

import argparse
import difflib
import os
import subprocess
import sys
from pathlib import Path
from datetime import datetime

try:
    from anthropic import Anthropic
except ImportError:
    print("ERROR: pip install anthropic", file=sys.stderr)
    sys.exit(1)


COMMENT_PREFIXES = {
    ".sh": ("#",),
    ".R": ("#",),
    ".r": ("#",),
    ".py": ("#",),
    ".pl": ("#",),
    ".awk": ("#",),
    ".c": ("//", "/*", " *", "*/"),
    ".cpp": ("//", "/*", " *", "*/"),
}


def is_comment_or_blank(line: str, ext: str) -> bool:
    s = line.strip()
    if not s:
        return True
    prefixes = COMMENT_PREFIXES.get(ext, ("#",))
    return any(s.startswith(p) for p in prefixes)


def diff_is_safe(original: str, updated: str, ext: str) -> tuple[bool, list[str]]:
    """
    Return (safe, reasons). Safe means every changed line is either a
    comment, blank, or a recognized echo/print statement.
    """
    orig_lines = original.splitlines()
    new_lines = updated.splitlines()
    diff = list(difflib.unified_diff(orig_lines, new_lines, lineterm=""))

    flags = []
    echo_markers = (
        "qc_log", "printf", "echo", "message(", "print(", "fprintf(stderr"
    )
    in_hunk = False
    for line in diff:
        if line.startswith("@@"):
            in_hunk = True
            continue
        if not in_hunk:
            continue
        if line.startswith("+") and not line.startswith("+++"):
            content = line[1:]
            if is_comment_or_blank(content, ext):
                continue
            if any(m in content for m in echo_markers):
                continue
            flags.append(f"added non-comment non-echo: {content[:100]}")
        elif line.startswith("-") and not line.startswith("---"):
            content = line[1:]
            if is_comment_or_blank(content, ext):
                continue
            # Removed logic lines are always a red flag
            flags.append(f"removed content: {content[:100]}")
    return (len(flags) == 0, flags)


def standardize_one(client: Anthropic, prompt: str, file_path: Path,
                    model: str, max_tokens: int = 16000) -> tuple[str, str, str]:
    """
    Return (diff_summary, risk_flags, updated_content).
    """
    content = file_path.read_text(errors="replace")
    msg = client.messages.create(
        model=model,
        max_tokens=max_tokens,
        messages=[{
            "role": "user",
            "content": (
                prompt + "\n\n---\n\n"
                + f"FILE PATH: {file_path}\n"
                + f"FILE EXTENSION: {file_path.suffix}\n\n"
                + f"FILE CONTENTS:\n```{file_path.suffix.lstrip('.')}\n"
                + content + "\n```\n"
            )
        }]
    )
    out = msg.content[0].text
    # Parse the three sections
    diff_summary = ""
    risk_flags = ""
    updated = ""
    if "**DIFF SUMMARY**" in out:
        after = out.split("**DIFF SUMMARY**", 1)[1]
        if "**RISK FLAGS**" in after:
            diff_summary, rest = after.split("**RISK FLAGS**", 1)
            if "**UPDATED FILE**" in rest:
                risk_flags, file_block = rest.split("**UPDATED FILE**", 1)
                # Extract the first code block from file_block
                if "```" in file_block:
                    parts = file_block.split("```")
                    # parts[0] is preamble, parts[1] is the code block,
                    # the language tag is the first line of parts[1]
                    if len(parts) >= 2:
                        block = parts[1]
                        lines = block.split("\n", 1)
                        if len(lines) > 1:
                            updated = lines[1]
    return (diff_summary.strip(), risk_flags.strip(), updated)


def main():
    ap = argparse.ArgumentParser(
        description="Batch apply standardized headers/echoes to scripts.")
    ap.add_argument("--file_list", required=True,
                    help="Text file, one script path per line")
    ap.add_argument("--prompt", required=True,
                    help="Path to the frozen super prompt (e.g. STANDARDIZE_SCRIPT.md)")
    ap.add_argument("--model", default="claude-opus-4-7")
    ap.add_argument("--api_key", default=os.environ.get("ANTHROPIC_API_KEY"))
    ap.add_argument("--staging_dir", default="./staged_updates")
    ap.add_argument("--log", default="standardize.log")
    ap.add_argument("--apply", action="store_true",
                    help="Actually modify files in place (default: staging only)")
    ap.add_argument("--dry_run", action="store_true",
                    help="Show plan, don't call API")
    args = ap.parse_args()

    if not args.api_key and not args.dry_run:
        ap.error("--api_key required or set ANTHROPIC_API_KEY env")

    files = [Path(ln.strip()) for ln in Path(args.file_list).read_text().splitlines()
             if ln.strip() and not ln.startswith("#")]
    print(f"[standardize] {len(files)} files queued", file=sys.stderr)

    if args.dry_run:
        for f in files:
            print(f"  would process: {f}")
        return

    prompt = Path(args.prompt).read_text()
    staging = Path(args.staging_dir)
    staging.mkdir(parents=True, exist_ok=True)
    log_f = open(args.log, "a")
    log_f.write(f"\n=== {datetime.utcnow().isoformat()} run ===\n")

    client = Anthropic(api_key=args.api_key)

    stats = {"processed": 0, "safe": 0, "flagged": 0, "applied": 0, "errors": 0}
    for f in files:
        if not f.exists():
            log_f.write(f"MISSING: {f}\n")
            stats["errors"] += 1
            continue
        try:
            diff, risks, updated = standardize_one(client, prompt, f, args.model)
        except Exception as e:
            log_f.write(f"ERROR {f}: {e}\n")
            stats["errors"] += 1
            continue

        if not updated:
            log_f.write(f"EMPTY_RESPONSE: {f}\n")
            stats["errors"] += 1
            continue

        original = f.read_text(errors="replace")
        safe, flags = diff_is_safe(original, updated, f.suffix)

        # Always stage
        stage_path = staging / f.name
        stage_path.write_text(updated)
        log_f.write(f"\n--- {f} ---\n")
        log_f.write(f"SAFE: {safe}\n")
        log_f.write(f"DIFF SUMMARY: {diff}\n")
        log_f.write(f"RISK FLAGS: {risks}\n")
        if not safe:
            log_f.write(f"AUTO-SAFETY FLAGS: {flags}\n")
            stats["flagged"] += 1
        else:
            stats["safe"] += 1

        if args.apply and safe:
            f.write_text(updated)
            stats["applied"] += 1
            print(f"  applied: {f}", file=sys.stderr)
        else:
            print(f"  staged ({'SAFE' if safe else 'FLAGGED'}): {f}", file=sys.stderr)

        stats["processed"] += 1

    log_f.write(f"\nFINAL: {stats}\n")
    log_f.close()
    print(f"\n[standardize] done: {stats}", file=sys.stderr)
    print(f"  staging: {staging}", file=sys.stderr)
    print(f"  log: {args.log}", file=sys.stderr)


if __name__ == "__main__":
    main()
