#!/usr/bin/env python3
"""
atlas_server.py — Reference local server for the Inversion Atlas.

Usage:
    python3 atlas_server.py /path/to/project_folder [--port 8765]

Provides a small HTTP API on http://localhost:PORT for atlas to use:

    GET  /health              → { "ok": true, "atlas_server": "v1" }
    GET  /file/<rel_path>     → file contents (text or json)
    POST /file/<rel_path>     → write body to file (creates dirs as needed)
    POST /compute/<name>      → run a named compute, return JSON

Designed for solo research use on localhost. NO authentication, NO HTTPS,
NO multi-user safety. Bind to 127.0.0.1 (default) so only the local
machine can reach it. If you need remote access, do it via SSH tunnel.

The compute names this reference server understands are stubs — they
return placeholder responses so atlas's adapter wiring can be tested
end-to-end. Real compute implementations (numpy contingency matrices,
etc.) are intentionally not included; you'll typically write project-
specific compute logic on top of this skeleton.

Why not Flask / FastAPI?
    Stdlib http.server is enough for one user on localhost, has zero
    install cost, and runs on any Python 3.7+ without a venv. If you
    later want async I/O or middleware, swap in FastAPI — the protocol
    is the same.
"""

import argparse
import datetime as dt
import json
import os
import sys
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import unquote


SERVER_VERSION = "v1"
DEFAULT_PORT = 8765


class AtlasHandler(BaseHTTPRequestHandler):
    project_root: Path = None  # set per-server-instance below

    # ---- CORS ---------------------------------------------------------------

    def _cors(self):
        # Allow atlas (which may be opened as file://...) to talk to us.
        # `*` is the laziest but safest choice for a localhost-only server.
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")

    def do_OPTIONS(self):
        self.send_response(204)
        self._cors()
        self.end_headers()

    # ---- helpers ------------------------------------------------------------

    def _send_json(self, status, payload):
        body = json.dumps(payload).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self._cors()
        self.end_headers()
        self.wfile.write(body)

    def _send_text(self, status, text, content_type="text/plain"):
        body = text.encode("utf-8") if isinstance(text, str) else text
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self._cors()
        self.end_headers()
        self.wfile.write(body)

    def _safe_path(self, rel: str) -> Path:
        """Resolve a request path inside the project root. Rejects
        traversal attempts (`..`) by checking the resolved path is under
        project_root."""
        p = (self.project_root / rel).resolve()
        # Path.is_relative_to is 3.9+; do the manual check for portability
        try:
            p.relative_to(self.project_root.resolve())
        except ValueError:
            raise PermissionError(f"path escapes project root: {rel}")
        return p

    # ---- routes -------------------------------------------------------------

    def do_GET(self):
        path = self.path.split("?", 1)[0]
        if path == "/health":
            self._send_json(200, {
                "ok": True,
                "atlas_server": SERVER_VERSION,
                "project_root": str(self.project_root),
                "now": dt.datetime.now().isoformat(timespec="seconds"),
            })
            return

        if path.startswith("/file/"):
            rel = unquote(path[len("/file/"):])
            try:
                target = self._safe_path(rel)
            except PermissionError as e:
                self._send_json(403, {"ok": False, "error": str(e)})
                return
            if not target.exists():
                self._send_json(404, {"ok": False, "error": f"not found: {rel}"})
                return
            if target.is_dir():
                # Directory listing as JSON
                entries = []
                for entry in sorted(target.iterdir()):
                    entries.append({
                        "name": entry.name,
                        "is_dir": entry.is_dir(),
                        "size": entry.stat().st_size if entry.is_file() else None,
                    })
                self._send_json(200, {"ok": True, "dir": rel, "entries": entries})
                return
            ct = "application/json" if target.suffix == ".json" else "text/plain"
            try:
                self._send_text(200, target.read_bytes(), content_type=ct)
            except Exception as e:
                self._send_json(500, {"ok": False, "error": str(e)})
            return

        self._send_json(404, {"ok": False, "error": f"unknown route: {path}"})

    def do_POST(self):
        path = self.path.split("?", 1)[0]
        length = int(self.headers.get("Content-Length") or 0)
        body = self.rfile.read(length) if length > 0 else b""

        if path.startswith("/file/"):
            rel = unquote(path[len("/file/"):])
            try:
                target = self._safe_path(rel)
            except PermissionError as e:
                self._send_json(403, {"ok": False, "error": str(e)})
                return
            try:
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(body)
                self._send_json(200, {"ok": True, "path": rel, "bytes": len(body)})
            except Exception as e:
                self._send_json(500, {"ok": False, "error": str(e)})
            return

        if path.startswith("/compute/"):
            name = unquote(path[len("/compute/"):])
            try:
                args = json.loads(body) if body else {}
            except Exception as e:
                self._send_json(400, {"ok": False, "error": f"invalid JSON: {e}"})
                return
            handler = COMPUTE_REGISTRY.get(name)
            if handler is None:
                self._send_json(404, {
                    "ok": False,
                    "error": f"unknown compute: {name}",
                    "registered": sorted(COMPUTE_REGISTRY.keys()),
                })
                return
            try:
                result = handler(args, project_root=self.project_root)
                self._send_json(200, {"ok": True, "result": result})
            except Exception as e:
                self._send_json(500, {"ok": False, "error": f"{type(e).__name__}: {e}"})
            return

        self._send_json(404, {"ok": False, "error": f"unknown route: {path}"})

    # Quiet logs — atlas hits /health every 30s
    def log_message(self, fmt, *args):
        if "/health" in fmt % args:
            return
        super().log_message(fmt, *args)


# -----------------------------------------------------------------------------
# Compute registry
#
# Add real implementations as you need them. Each function takes
# (args: dict, project_root: Path) and returns a JSON-serializable result.
# The atlas adapter calls these via POST /compute/<name>.
# -----------------------------------------------------------------------------

def _compute_echo(args, project_root):
    """Sanity check: echoes args + adds metadata. Useful for verifying
    the round-trip from atlas works end-to-end."""
    return {
        "echo": args,
        "project_root": str(project_root),
        "server_version": SERVER_VERSION,
    }


def _compute_list_files(args, project_root):
    """List files matching a glob under project_root. Used by the export
    UI to populate file-pickers."""
    pattern = args.get("pattern", "*")
    matches = sorted(str(p.relative_to(project_root))
                     for p in project_root.glob(pattern)
                     if p.is_file())
    return {"pattern": pattern, "files": matches}


COMPUTE_REGISTRY = {
    "echo": _compute_echo,
    "list_files": _compute_list_files,
    # Add: "contingency_matrix", "inheritance_full", "linkage_disequilibrium", ...
    # Each one a function (args, project_root) -> JSON-serializable.
}


# -----------------------------------------------------------------------------
# Entry point
# -----------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(description="Atlas local server")
    parser.add_argument("project_root", type=Path,
                        help="Folder atlas will read/write under")
    parser.add_argument("--port", type=int, default=DEFAULT_PORT)
    parser.add_argument("--host", default="127.0.0.1",
                        help="Bind address. Default 127.0.0.1 (localhost only)")
    args = parser.parse_args(argv)

    root = args.project_root.resolve()
    if not root.is_dir():
        print(f"error: {root} is not a directory", file=sys.stderr)
        return 2

    AtlasHandler.project_root = root
    server = ThreadingHTTPServer((args.host, args.port), AtlasHandler)
    print(f"atlas_server.py {SERVER_VERSION}")
    print(f"  project_root: {root}")
    print(f"  listening on: http://{args.host}:{args.port}")
    print(f"  computes:     {sorted(COMPUTE_REGISTRY.keys())}")
    print(f"  Ctrl+C to stop.")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nstopping...")
        server.shutdown()
    return 0


if __name__ == "__main__":
    sys.exit(main())
