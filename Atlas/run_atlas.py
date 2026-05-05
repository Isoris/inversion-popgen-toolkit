#!/usr/bin/env python3
"""
run_atlas.py — single-command launcher for the merged atlas server (turn 145).

Cross-platform equivalent of run_atlas.sh. Use this on Windows or wherever
bash isn't convenient.

    python3 run_atlas.py
    python3 run_atlas.py --port 8765
    python3 run_atlas.py --no-popstats
    python3 run_atlas.py --check               # dry run
    python3 run_atlas.py --config custom.yaml

Behavior:
  - The /file + /compute subsystem is always enabled with project_root set
    to this script's parent directory (the Atlas/ root).
  - The popstats subsystem is enabled iff
    server_turn1/popstats_server.config.yaml exists, or --config <path> is
    passed.
  - Default bind: 127.0.0.1:8765. Override via --host / --port.

See docs/SERVER_AUDIT_2026-05-05.md for the full architecture map.
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path

ATLAS_ROOT = Path(__file__).resolve().parent
SERVER_DIR = ATLAS_ROOT / "server_turn1"
SERVER_PY = SERVER_DIR / "popstats_server.py"
DEFAULT_CFG = SERVER_DIR / "popstats_server.config.yaml"
EXAMPLE_CFG = SERVER_DIR / "popstats_server.config.example.yaml"


def _err(msg: str, code: int = 2) -> None:
    print(f"error: {msg}", file=sys.stderr)
    sys.exit(code)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description="Launch the merged atlas server.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--no-popstats", action="store_true",
                   help="skip popstats even if a config file is present")
    p.add_argument("--config", type=Path, default=None,
                   help="non-default popstats config path")
    p.add_argument("--host", default="127.0.0.1", help="bind host (default 127.0.0.1)")
    p.add_argument("--port", default=8765, type=int, help="bind port (default 8765)")
    p.add_argument("--check", action="store_true",
                   help="dry run — print what would happen and exit 0")
    args = p.parse_args(argv)

    if not SERVER_PY.is_file():
        _err(f"merged server not found at {SERVER_PY}\n"
             f"  expected layout: <atlas_root>/server_turn1/popstats_server.py")

    cfg_to_use: "Path | None" = None
    if not args.no_popstats:
        if args.config is not None:
            if not args.config.is_file():
                _err(f"--config file not found: {args.config}")
            cfg_to_use = args.config.resolve()
        elif DEFAULT_CFG.is_file():
            cfg_to_use = DEFAULT_CFG

    # Preflight: importable deps. We check this in the SAME interpreter
    # that will exec the child server, since that's the one whose
    # site-packages matters.
    try:
        import fastapi  # noqa: F401
        import uvicorn  # noqa: F401
    except ImportError:
        _err("required Python packages missing.\n"
             f"  install: {sys.executable} -m pip install "
             f"fastapi uvicorn pyyaml numpy pandas pydantic\n"
             f"  (see {SERVER_DIR / 'requirements.txt'})", code=3)

    print("=" * 73)
    print(" atlas server launcher")
    print("=" * 73)
    print(f"  Atlas root:    {ATLAS_ROOT}")
    print(f"  Server:        {SERVER_PY}")
    print(f"  Bind:          http://{args.host}:{args.port}")
    print(f"  Project root:  {ATLAS_ROOT}        (/file + /compute subsystem)")
    if cfg_to_use is not None:
        print(f"  Popstats cfg:  {cfg_to_use}        (popstats subsystem)")
    else:
        print(f"  Popstats:      disabled (no config provided)")
        if not args.no_popstats and args.config is None:
            print(f"                 to enable:")
            print(f"                   cp {EXAMPLE_CFG} \\")
            print(f"                      {DEFAULT_CFG}")
            print(f"                   then edit paths and re-run this launcher.")
    print("=" * 73)

    if args.check:
        print("(--check) dry run only; not starting server")
        return 0

    # Build CLI for the merged server. We invoke as a subprocess so its
    # uvicorn "import string" mode resolves the right module by name.
    cmd = [
        sys.executable, str(SERVER_PY),
        "--project-root", str(ATLAS_ROOT),
        "--host", args.host,
        "--port", str(args.port),
    ]
    if cfg_to_use is not None:
        cmd.extend(["--config", str(cfg_to_use)])

    # Run from server_turn1/ so popstats_server's relative imports resolve.
    if hasattr(os, "execvpe"):
        os.chdir(SERVER_DIR)
        os.execvpe(cmd[0], cmd, os.environ)
    else:
        # Windows fallback: subprocess + propagate exit code.
        import subprocess
        return subprocess.call(cmd, cwd=str(SERVER_DIR))

    return 0   # unreachable on POSIX (execvpe replaces process)


if __name__ == "__main__":
    sys.exit(main())
