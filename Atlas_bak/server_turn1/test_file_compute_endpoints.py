"""Tests for the /file + /compute subsystem (merged from atlas_server.py
into popstats_server.py during turn 145 server-unify).

Independent of the popstats subsystem — these run with no popstats config
loaded; the merged server treats the two subsystems as orthogonal.

Run: python3 server_turn1/test_file_compute_endpoints.py
"""
from __future__ import annotations

import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "server_turn1"))

from fastapi.testclient import TestClient

import popstats_server as ps


def _fresh_client(project_root: Path) -> TestClient:
    """Build a fresh TestClient with the file subsystem bootstrapped to the
    given project root and the popstats subsystem deliberately left unset.
    """
    # Reset module-level state between tests so we can swap project roots.
    ps.PROJECT_ROOT = None
    ps.ENGINES = None
    ps.CACHE = None
    ps.SAMPLES = None
    ps._bootstrap_file(project_root)
    return TestClient(ps.app)


class TestHealthEndpoint(unittest.TestCase):
    def test_health_reports_both_subsystems(self):
        with tempfile.TemporaryDirectory() as td:
            client = _fresh_client(Path(td))
            r = client.get("/health")
            self.assertEqual(r.status_code, 200)
            j = r.json()
            self.assertTrue(j["ok"])
            self.assertEqual(j["service"], "atlas_server (merged)")
            self.assertIn("subsystems", j)
            self.assertIn("file", j["subsystems"])
            self.assertIn("popstats", j["subsystems"])
            # File subsystem is up
            self.assertTrue(j["subsystems"]["file"]["ready"])
            self.assertEqual(
                Path(j["subsystems"]["file"]["project_root"]).resolve(),
                Path(td).resolve())
            # Popstats subsystem is down (no --config)
            self.assertFalse(j["subsystems"]["popstats"]["ready"])
            self.assertIn("reason", j["subsystems"]["popstats"])

    def test_health_alias(self):
        """Both /health and /api/health return identical payloads."""
        with tempfile.TemporaryDirectory() as td:
            client = _fresh_client(Path(td))
            a = client.get("/health").json()
            b = client.get("/api/health").json()
            # Ignore `now` field which differs by milliseconds
            a.pop("now", None); b.pop("now", None)
            self.assertEqual(a, b)


class TestFileEndpointGET(unittest.TestCase):
    def test_get_file_text(self):
        with tempfile.TemporaryDirectory() as td:
            (Path(td) / "hello.txt").write_text("hi there")
            client = _fresh_client(Path(td))
            r = client.get("/file/hello.txt")
            self.assertEqual(r.status_code, 200)
            self.assertEqual(r.text, "hi there")
            self.assertTrue(r.headers["content-type"].startswith("text/plain"))

    def test_get_file_json(self):
        with tempfile.TemporaryDirectory() as td:
            (Path(td) / "data.json").write_text('{"k": 1}')
            client = _fresh_client(Path(td))
            r = client.get("/file/data.json")
            self.assertEqual(r.status_code, 200)
            self.assertEqual(r.headers["content-type"], "application/json")
            self.assertEqual(r.json(), {"k": 1})

    def test_get_file_nested_path(self):
        with tempfile.TemporaryDirectory() as td:
            nested = Path(td) / "a" / "b" / "c.json"
            nested.parent.mkdir(parents=True)
            nested.write_text('{"deep": true}')
            client = _fresh_client(Path(td))
            r = client.get("/file/a/b/c.json")
            self.assertEqual(r.status_code, 200)
            self.assertEqual(r.json(), {"deep": True})

    def test_get_directory_lists_entries(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td) / "json"
            d.mkdir()
            (d / "a.json").write_text("{}")
            (d / "b.json").write_text("{}")
            (d / "subdir").mkdir()
            client = _fresh_client(Path(td))
            r = client.get("/file/json")
            self.assertEqual(r.status_code, 200)
            j = r.json()
            self.assertTrue(j["ok"])
            self.assertEqual(j["dir"], "json")
            self.assertEqual(len(j["entries"]), 3)
            names = {e["name"] for e in j["entries"]}
            self.assertEqual(names, {"a.json", "b.json", "subdir"})
            # Entries are sorted
            self.assertEqual([e["name"] for e in j["entries"]],
                             ["a.json", "b.json", "subdir"])

    def test_get_missing_file_returns_404(self):
        with tempfile.TemporaryDirectory() as td:
            client = _fresh_client(Path(td))
            r = client.get("/file/does_not_exist.json")
            self.assertEqual(r.status_code, 404)

    def test_get_path_traversal_blocked(self):
        """An encoded `../` in the path must not escape PROJECT_ROOT."""
        with tempfile.TemporaryDirectory() as td:
            client = _fresh_client(Path(td))
            r = client.get("/file/..%2F..%2Fetc%2Fpasswd")
            self.assertEqual(r.status_code, 403)


class TestFileEndpointPOST(unittest.TestCase):
    def test_post_creates_file(self):
        with tempfile.TemporaryDirectory() as td:
            client = _fresh_client(Path(td))
            r = client.post("/file/written.txt", content=b"hello")
            self.assertEqual(r.status_code, 200)
            self.assertEqual(r.json(), {"ok": True, "path": "written.txt", "bytes": 5})
            self.assertEqual((Path(td) / "written.txt").read_text(), "hello")

    def test_post_creates_intermediate_dirs(self):
        with tempfile.TemporaryDirectory() as td:
            client = _fresh_client(Path(td))
            r = client.post("/file/a/b/c.json", content=b'{"deep": true}')
            self.assertEqual(r.status_code, 200)
            target = Path(td) / "a" / "b" / "c.json"
            self.assertTrue(target.exists())
            self.assertEqual(json.loads(target.read_text()), {"deep": True})

    def test_post_overwrites_existing(self):
        with tempfile.TemporaryDirectory() as td:
            (Path(td) / "f.txt").write_text("OLD")
            client = _fresh_client(Path(td))
            r = client.post("/file/f.txt", content=b"NEW")
            self.assertEqual(r.status_code, 200)
            self.assertEqual((Path(td) / "f.txt").read_text(), "NEW")

    def test_post_path_traversal_blocked(self):
        with tempfile.TemporaryDirectory() as td:
            client = _fresh_client(Path(td))
            r = client.post("/file/..%2Fevil.txt", content=b"x")
            self.assertEqual(r.status_code, 403)


class TestComputeEndpoint(unittest.TestCase):
    def test_echo_compute(self):
        with tempfile.TemporaryDirectory() as td:
            client = _fresh_client(Path(td))
            r = client.post("/compute/echo", json={"foo": "bar", "n": 42})
            self.assertEqual(r.status_code, 200)
            j = r.json()
            self.assertTrue(j["ok"])
            self.assertEqual(j["result"]["echo"], {"foo": "bar", "n": 42})
            self.assertIn("project_root", j["result"])
            self.assertIn("server_version", j["result"])

    def test_echo_no_body(self):
        """Empty POST should be treated as args={}."""
        with tempfile.TemporaryDirectory() as td:
            client = _fresh_client(Path(td))
            r = client.post("/compute/echo")
            self.assertEqual(r.status_code, 200)
            self.assertEqual(r.json()["result"]["echo"], {})

    def test_list_files_compute(self):
        with tempfile.TemporaryDirectory() as td:
            (Path(td) / "a.json").write_text("{}")
            (Path(td) / "b.json").write_text("{}")
            (Path(td) / "c.txt").write_text("hi")
            client = _fresh_client(Path(td))
            r = client.post("/compute/list_files", json={"pattern": "*.json"})
            self.assertEqual(r.status_code, 200)
            files = r.json()["result"]["files"]
            self.assertEqual(set(files), {"a.json", "b.json"})

    def test_unknown_compute_returns_404(self):
        with tempfile.TemporaryDirectory() as td:
            client = _fresh_client(Path(td))
            r = client.post("/compute/totally_made_up", json={})
            self.assertEqual(r.status_code, 404)
            # FastAPI wraps the dict in {"detail": …}
            detail = r.json()["detail"]
            self.assertIn("unknown compute", str(detail))

    def test_invalid_json_returns_400(self):
        with tempfile.TemporaryDirectory() as td:
            client = _fresh_client(Path(td))
            r = client.post("/compute/echo", content=b"this is not json{{{",
                            headers={"content-type": "application/json"})
            self.assertEqual(r.status_code, 400)


class TestSubsystemIsolation(unittest.TestCase):
    """The point of the merge: each subsystem is independently usable."""

    def test_file_works_when_popstats_disabled(self):
        with tempfile.TemporaryDirectory() as td:
            client = _fresh_client(Path(td))
            (Path(td) / "x.txt").write_text("hi")
            self.assertEqual(client.get("/file/x.txt").status_code, 200)
            self.assertEqual(client.post("/compute/echo", json={}).status_code, 200)
            # Popstats must 503
            self.assertEqual(client.get("/api/cache/keys").status_code, 503)

    def test_file_503_when_subsystem_disabled(self):
        """If neither subsystem is up, /file returns 503."""
        # Force-clear PROJECT_ROOT
        ps.PROJECT_ROOT = None
        client = TestClient(ps.app)
        r = client.get("/file/anything.json")
        self.assertEqual(r.status_code, 503)
        self.assertIn("subsystem not configured", r.json()["detail"])

    def test_compute_503_when_subsystem_disabled(self):
        ps.PROJECT_ROOT = None
        client = TestClient(ps.app)
        r = client.post("/compute/echo", json={})
        self.assertEqual(r.status_code, 503)


if __name__ == "__main__":
    unittest.main(verbosity=2)
