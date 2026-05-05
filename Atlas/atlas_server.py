#!/usr/bin/env python3
"""
atlas_server.py — compatibility shim (turn 145).

Pre-turn-145 this was a standalone http.server-based localhost server
serving /health, /file/<path>, and /compute/<name>. Those endpoints now
live in server_turn1/popstats_server.py as the "/file + /compute"
subsystem of the merged atlas server.

This shim keeps the legacy entry point working — any script or doc that
says "run python3 atlas_server.py <project_root>" continues to function
unchanged. Internally it just calls the merged server with
--project-root set to the positional argument.

Usage (unchanged from pre-turn-145):

    python3 atlas_server.py /path/to/project_folder [--port 8765] [--host 127.0.0.1]

For the recommended flow that also starts popstats, use:

    ./run_atlas.sh                       # POSIX shells
    python3 run_atlas.py                 # cross-platform

Both launchers wrap the merged server with sensible defaults and pick up
server_turn1/popstats_server.config.yaml automatically when present.

See docs/SERVER_AUDIT_2026-05-05.md for the full architecture.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

DEFAULT_PORT = 8765
SERVER_VERSION = "v2-shim-turn145"

ATLAS_ROOT = Path(__file__).resolve().parent
SERVER_DIR = ATLAS_ROOT / "server_turn1"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Atlas server shim — delegates to the merged server")
    parser.add_argument("project_root", type=Path,
                        help="Folder atlas will read/write under "
                             "(passed as --project-root to the merged server)")
    parser.add_argument("--port", type=int, default=DEFAULT_PORT)
    parser.add_argument("--host", default="127.0.0.1",
                        help="Bind address. Default 127.0.0.1 (localhost only)")
    parser.add_argument("--config", default=None, type=Path,
                        help="optional popstats config (turn 145 — also "
                             "spins up popstats subsystem)")
    args = parser.parse_args(argv)

    root = args.project_root.resolve()
    if not root.is_dir():
        print(f"error: {root} is not a directory", file=sys.stderr)
        return 2

    if not SERVER_DIR.is_dir():
        print(f"error: server_turn1/ not found at {SERVER_DIR}", file=sys.stderr)
        print("  this shim expects the merged server at "
              "<atlas_root>/server_turn1/popstats_server.py", file=sys.stderr)
        return 2

    # Import the merged server with its expected sys.path
    sys.path.insert(0, str(SERVER_DIR))
    try:
        import popstats_server as ps  # type: ignore
    except ImportError as e:
        print(f"error: failed to import merged server: {e}", file=sys.stderr)
        print("  install fastapi + uvicorn + pyyaml + numpy + pandas, "
              "then retry", file=sys.stderr)
        return 3

    # Bring up the /file + /compute subsystem with the supplied project root
    ps._bootstrap_file(root)

    # Optionally bring up popstats
    if args.config is not None:
        try:
            ps._bootstrap(args.config)
        except Exception as e:
            print(f"warning: --config bootstrap failed: {e}", file=sys.stderr)
            print("  popstats endpoints will return 503", file=sys.stderr)

    print(f"atlas_server.py shim {SERVER_VERSION}")
    print(f"  project_root: {root}")
    print(f"  listening on: http://{args.host}:{args.port}")
    print(f"  computes:     {sorted(ps.COMPUTE_REGISTRY.keys())}")
    if args.config is not None:
        print(f"  popstats:     {'enabled' if ps.ENGINES else 'failed'}")
    else:
        print(f"  popstats:     disabled (no --config)")
    print(f"  Ctrl+C to stop.")

    import uvicorn
    uvicorn.run(ps.app, host=args.host, port=args.port, log_level="info")
    return 0


if __name__ == "__main__":
    sys.exit(main())
