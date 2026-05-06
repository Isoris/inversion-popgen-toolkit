#!/usr/bin/env python3
# =============================================================================
# popstats_server.py
# =============================================================================
# Live-engine HTTP wrapper for inversion popstats. Sits on the cluster (LANTA)
# alongside BEAGLE / BAM / Engine B caches, exposes four endpoints to the atlas:
#
#   POST /api/popstats/groupwise        wraps region_popstats (Engine F)
#   POST /api/popstats/hobs_groupwise   wraps angsd_fixed_HWE -doHWE + hobs_windower
#   POST /api/ancestry/groupwise_q      reads instant_q local_Q_samples cache
#   POST /api/shelf_ld_test             optional server-side; atlas does this in JS
#   GET  /api/health                    engine versions + cache stats
#   GET  /api/cache/keys                debug: list cache hashes
#   DEL  /api/cache/keys/<hash>         debug: drop a single cache entry
#   GET  /api/jobs/<job_id>             progress polling for slow runs
#
# Architecture notes:
#
#   - The atlas is the source of truth for sample groupings. Every request
#     carries explicit member lists; the server NEVER runs k-means or
#     re-derives group assignments. This keeps the server stateless and
#     cleanly separates "what is a group" (atlas) from "compute stats given
#     a group" (server).
#
#   - Cache keys are content-addressable: sha256(chrom + region + sorted
#     groups + metric set + win/step + engine_binary_hash). When you
#     recompile region_popstats / hobs_windower / angsd, every dependent
#     cache entry invalidates automatically.
#
#   - Q07b first-run cost is real (5-15 s per group of ANGSD). Endpoint
#     returns 202 + job_id when a fresh ANGSD run is required; client polls
#     /api/jobs/<id> until done.
#
#   - This file is dependency-light by design: FastAPI + uvicorn + pyyaml
#     + numpy + pandas. It does NOT depend on any cluster R packages or the
#     atlas codebase. Run it as a standalone process.
# =============================================================================

from __future__ import annotations

import asyncio
import dataclasses
import gzip
import hashlib
import io
import json
import logging
import os
import re
import shutil
import string
import subprocess
import sys
import tempfile
import time
import uuid
from collections import OrderedDict
from glob import glob
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import yaml
from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, StreamingResponse, Response
from pydantic import BaseModel, Field, field_validator

# Turn 11a: LD split-heatmap endpoint (separate module to keep this file readable)
from ld_endpoint import LDSplitReq, handle_split_heatmap

# =============================================================================
# Logging
# =============================================================================

LOG_FMT = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
logging.basicConfig(level=logging.INFO, format=LOG_FMT, stream=sys.stderr)
log = logging.getLogger("popstats")

# =============================================================================
# Config loader
# =============================================================================
# YAML config with ${var} interpolation. ${var} can be either another config
# key (dotted path resolved against the same YAML), or an environment var.
# We resolve config-keys first, then env, leaving unresolved markers as-is
# (which then surface as a clear "no such file" error downstream).

_VAR_RE = re.compile(r"\$\{([a-zA-Z_][a-zA-Z0-9_]*)\}")


def _interpolate(value: Any, root: Dict[str, Any]) -> Any:
    """Recursively expand ${name} markers in strings inside dicts/lists."""
    if isinstance(value, str):
        # Two passes: substitute config-key references, then env vars.
        # Config-key references are flat (no dotted lookup) — sufficient
        # for our case where everything is hung off a single 'base'.
        prev = None
        cur = value
        # bound iterations to avoid pathological recursion in a malformed config
        for _ in range(10):
            if cur == prev:
                break
            prev = cur

            def _sub(m: re.Match[str]) -> str:
                key = m.group(1)
                if key in root and isinstance(root[key], (str, int, float)):
                    return str(root[key])
                env = os.environ.get(key)
                if env is not None:
                    return env
                return m.group(0)  # unresolved — leave for downstream error

            cur = _VAR_RE.sub(_sub, cur)
        return cur
    if isinstance(value, dict):
        return {k: _interpolate(v, root) for k, v in value.items()}
    if isinstance(value, list):
        return [_interpolate(v, root) for v in value]
    return value


def load_config(path: Path) -> Dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Config not found: {path}")
    with open(path) as fh:
        raw = yaml.safe_load(fh) or {}

    # Two-pass interpolation: first pass resolves nested keys against the
    # root, second pass picks up anything that became resolvable after pass 1.
    resolved = _interpolate(raw, raw)
    resolved = _interpolate(resolved, resolved)
    return resolved


# =============================================================================
# Engine registry (paths + content hashes, lazily computed)
# =============================================================================

class EngineRegistry:
    """Wraps engine-binary paths + their sha256s. Hashes are cached in-process
    and refreshed at most once per `ttl_sec`. We hash binary contents (not
    mtimes) so a host-specific recompile reliably invalidates downstream
    caches even when atime/ctime get clobbered by network FS.
    """

    def __init__(self, paths: Dict[str, str], ttl_sec: float = 60.0):
        self._paths = {k: Path(v) for k, v in paths.items()}
        self._hashes: Dict[str, str] = {}
        self._last: Dict[str, float] = {}
        self._ttl = ttl_sec

    def path(self, name: str) -> Path:
        if name not in self._paths:
            raise KeyError(f"engine '{name}' not in registry")
        return self._paths[name]

    def hash(self, name: str) -> str:
        p = self.path(name)
        last = self._last.get(name, 0.0)
        if (time.time() - last) < self._ttl and name in self._hashes:
            return self._hashes[name]
        if not p.exists():
            # Surface as a clear "not compiled" error rather than silently
            # returning an empty hash that would poison downstream cache keys.
            return f"MISSING:{p}"
        h = hashlib.sha256()
        with open(p, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        digest = h.hexdigest()[:16]  # truncate — collision risk negligible at this scale
        self._hashes[name] = digest
        self._last[name] = time.time()
        return digest

    def all_hashes(self) -> Dict[str, str]:
        return {name: self.hash(name) for name in self._paths}


# =============================================================================
# Cache layer (content-addressable, on-disk, LRU eviction)
# =============================================================================

class DiskCache:
    """Two-tier cache: in-memory OrderedDict for the hash → metadata index,
    on-disk JSON files for the actual payloads. Eviction is global LRU
    (touched on read) bounded by `max_bytes`.

    Layout:
        <root>/
            popstats/<hash>.json      ← popstats/groupwise results
            hobs/<hash>.json          ← hobs/hobs_groupwise results
            ancestry/<hash>.json      ← ancestry/groupwise_q results
            hwe/<hwe_hash>.hwe.gz     ← raw ANGSD HWE per-group output
                                        (separate dir; outlives popstats cache)
            index.jsonl               ← append-only log of {hash, kind, bytes, ts}
    """

    def __init__(self, root: Path, max_bytes: int):
        self.root = Path(root)
        self.max_bytes = int(max_bytes)
        self.root.mkdir(parents=True, exist_ok=True)
        for sub in ("popstats", "hobs", "ancestry", "hwe", "shelf_ld"):
            (self.root / sub).mkdir(exist_ok=True)
        # index: { hash -> {kind, bytes, ts_used} }
        self._index: "OrderedDict[str, Dict[str, Any]]" = OrderedDict()
        self._index_path = self.root / "index.jsonl"
        self._reload_index()

    def _reload_index(self) -> None:
        """Rebuild the in-memory index from existing on-disk files. The
        append-only log is a hint, not authoritative — we trust the
        filesystem. Order is "oldest first" by mtime so LRU eviction
        starts with the genuinely-old entries on first run."""
        for sub in ("popstats", "hobs", "ancestry", "shelf_ld"):
            d = self.root / sub
            if not d.exists():
                continue
            entries = []
            for f in d.glob("*.json"):
                try:
                    st = f.stat()
                    entries.append((f, st.st_mtime, st.st_size))
                except OSError:
                    continue
            entries.sort(key=lambda x: x[1])
            for f, mt, sz in entries:
                self._index[f.stem] = {"kind": sub, "bytes": sz, "ts_used": mt}
        log.info("cache: loaded %d entries (%.1f MB) from %s",
                 len(self._index), self.total_bytes() / 1e6, self.root)

    def total_bytes(self) -> int:
        return sum(int(v.get("bytes", 0)) for v in self._index.values())

    def _path_for(self, key: str, kind: str, ext: str = "json") -> Path:
        return self.root / kind / f"{key}.{ext}"

    def has(self, key: str) -> bool:
        return key in self._index

    def kind(self, key: str) -> Optional[str]:
        e = self._index.get(key)
        return e["kind"] if e else None

    def get_json(self, key: str, kind: str) -> Optional[Dict[str, Any]]:
        p = self._path_for(key, kind, "json")
        if not p.exists():
            self._index.pop(key, None)
            return None
        try:
            with open(p) as fh:
                payload = json.load(fh)
        except (OSError, json.JSONDecodeError) as e:
            log.warning("cache: corrupt %s, dropping: %s", p, e)
            try:
                p.unlink()
            except OSError:
                pass
            self._index.pop(key, None)
            return None
        # LRU touch
        self._index.move_to_end(key)
        self._index[key]["ts_used"] = time.time()
        return payload

    def put_json(self, key: str, kind: str, payload: Dict[str, Any]) -> None:
        p = self._path_for(key, kind, "json")
        # Atomic write: tmp + rename
        with tempfile.NamedTemporaryFile(
            mode="w", dir=str(p.parent), prefix=f".{key}.", suffix=".tmp",
            delete=False, encoding="utf-8") as tf:
            json.dump(payload, tf, separators=(",", ":"), allow_nan=False)
            tmp_path = Path(tf.name)
        try:
            tmp_path.replace(p)
        except OSError as e:
            tmp_path.unlink(missing_ok=True)
            raise e
        sz = p.stat().st_size
        self._index[key] = {"kind": kind, "bytes": sz, "ts_used": time.time()}
        self._index.move_to_end(key)
        self._evict_if_needed()

    def drop(self, key: str) -> bool:
        e = self._index.pop(key, None)
        if not e:
            return False
        p = self._path_for(key, e["kind"], "json")
        try:
            p.unlink()
        except OSError:
            pass
        return True

    def _evict_if_needed(self) -> None:
        while self.total_bytes() > self.max_bytes and self._index:
            oldest_key, _ = next(iter(self._index.items()))
            self.drop(oldest_key)

    def keys_with_prefix(self, prefix: str = "") -> List[Dict[str, Any]]:
        return [
            {"hash": k, **v}
            for k, v in self._index.items()
            if k.startswith(prefix)
        ]

    # -- HWE-specific (gzipped, not JSON) --
    # Q07b's per-group .hwe.gz is large and binary. Lives in the same cache
    # tree under hwe/ but tracked separately because its lifetime is
    # decoupled from the windowed JSON results.
    def hwe_path(self, hwe_key: str) -> Path:
        return self.root / "hwe" / f"{hwe_key}.hwe.gz"

    def hwe_has(self, hwe_key: str) -> bool:
        p = self.hwe_path(hwe_key)
        return p.exists() and p.stat().st_size > 0


# =============================================================================
# Cache key derivation
# =============================================================================
# Stable hashing requires canonical serialization. We sort group names + member
# lists so the atlas can send groups in any order with identical cache hits.

def _canonical_groups(groups: Dict[str, List[str]]) -> List[Tuple[str, List[str]]]:
    return sorted(((name, sorted(members)) for name, members in groups.items()),
                  key=lambda x: x[0])


def _hash_obj(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(payload.encode()).hexdigest()[:32]


def popstats_cache_key(
    chrom: str,
    region: Optional[Dict[str, int]],
    groups: Dict[str, List[str]],
    metrics: List[str],
    win_bp: int,
    step_bp: int,
    win_type: int,
    downsample: int,
    engines: Dict[str, str],
) -> str:
    return _hash_obj({
        "kind": "popstats_groupwise.v1",
        "chrom": chrom,
        "region": region,
        "groups": _canonical_groups(groups),
        "metrics": sorted(metrics),
        "win_bp": int(win_bp), "step_bp": int(step_bp),
        "type": int(win_type), "downsample": int(downsample),
        "engine_region_popstats": engines.get("region_popstats", "?"),
    })


def hobs_cache_key(
    chrom: str,
    region: Optional[Dict[str, int]],
    groups: Dict[str, List[str]],
    scales: List[str],
    engines: Dict[str, str],
) -> str:
    return _hash_obj({
        "kind": "hobs_groupwise.v1",
        "chrom": chrom,
        "region": region,
        "groups": _canonical_groups(groups),
        "scales": sorted(scales),
        "engine_hobs_windower": engines.get("hobs_windower", "?"),
        "engine_angsd": engines.get("angsd_patched", "?"),
    })


def ancestry_q_cache_key(
    chrom: str,
    region: Optional[Dict[str, int]],
    groups: Dict[str, List[str]],
    K: int,
    scale: str,
) -> str:
    return _hash_obj({
        "kind": "ancestry_groupwise_q.v1",
        "chrom": chrom,
        "region": region,
        "groups": _canonical_groups(groups),
        "K": int(K),
        "scale": scale,
    })


def hwe_cache_key(
    chrom: str,
    members: List[str],
    angsd_params: Dict[str, Any],
    engine_hash: str,
) -> str:
    """One HWE blob per (chrom, sorted member set, angsd config, binary hash)."""
    return _hash_obj({
        "kind": "hwe.v1",
        "chrom": chrom,
        "members": sorted(members),
        "params": {k: angsd_params[k] for k in sorted(angsd_params)},
        "engine_angsd": engine_hash,
    })


# =============================================================================
# Job manager (for slow Q07b ANGSD runs)
# =============================================================================
# Atlas POSTs hobs_groupwise → if HWE cache miss, server returns 202 + job_id.
# Client polls /api/jobs/<id> until status == "done" then re-POSTs (or waits
# until status carries the result inline). We keep jobs in-process; on server
# restart they're forgotten (atlas re-issues the request).

@dataclasses.dataclass
class JobState:
    job_id: str
    kind: str
    chrom: str
    status: str = "queued"      # queued | running | done | error
    progress: float = 0.0       # 0..1
    message: str = ""
    started_at: float = 0.0
    finished_at: float = 0.0
    result: Optional[Dict[str, Any]] = None
    error: Optional[str] = None
    _per_group: Dict[str, str] = dataclasses.field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "job_id": self.job_id,
            "kind": self.kind,
            "chrom": self.chrom,
            "status": self.status,
            "progress": round(self.progress, 3),
            "message": self.message,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
            "per_group": dict(self._per_group),
            "error": self.error,
            # We DON'T inline the result here — clients fetch it explicitly
            # via /api/jobs/<id>?include_result=1 once status == 'done'.
            "has_result": self.result is not None,
        }


class JobManager:
    def __init__(self, max_jobs: int = 200):
        self._jobs: "OrderedDict[str, JobState]" = OrderedDict()
        self._max = max_jobs

    def new(self, kind: str, chrom: str) -> JobState:
        jid = uuid.uuid4().hex[:12]
        st = JobState(job_id=jid, kind=kind, chrom=chrom, status="queued",
                      started_at=time.time())
        self._jobs[jid] = st
        # bound memory: drop oldest if too many
        while len(self._jobs) > self._max:
            self._jobs.popitem(last=False)
        return st

    def get(self, jid: str) -> Optional[JobState]:
        return self._jobs.get(jid)


# =============================================================================
# Sample-list resolver
# =============================================================================
# region_popstats requires a sample_list arg matching BAM-list order. The atlas
# sends sample IDs as group members; we have to validate they exist in the
# canonical sample_list and reject unknowns clearly.

class SampleListIndex:
    def __init__(self, path: Path):
        self.path = Path(path)
        if not self.path.exists():
            raise FileNotFoundError(f"sample_list not found: {self.path}")
        with open(self.path) as fh:
            self._ids = [line.strip() for line in fh if line.strip()]
        self._set = set(self._ids)
        log.info("sample_list: %d ids from %s", len(self._ids), self.path)

    def __len__(self) -> int:
        return len(self._ids)

    def has(self, sid: str) -> bool:
        return sid in self._set

    def filter_known(self, ids: List[str]) -> Tuple[List[str], List[str]]:
        kept: List[str] = []
        dropped: List[str] = []
        for s in ids:
            (kept if s in self._set else dropped).append(s)
        return kept, dropped


# =============================================================================
# BEAGLE resolver
# =============================================================================
# Mirrors run_chrom.sh + Q07's resolution: glob ${BEAGLE_DIR}/*<chrom>*beagle.gz,
# prefer the dense (non-thin) variant.

def resolve_beagle(beagle_dir: Path, chrom: str) -> Path:
    primary = beagle_dir / f"main_qcpass.{chrom}.beagle.gz"
    if primary.exists():
        return primary
    candidates = sorted(beagle_dir.glob(f"*{chrom}*.beagle.gz"))
    # Filter out thin if dense is also present
    dense = [c for c in candidates if "thin" not in c.name]
    pool = dense if dense else candidates
    if not pool:
        raise FileNotFoundError(
            f"No BEAGLE for chrom={chrom} in {beagle_dir}")
    return pool[0]


def resolve_chrom_size(ref_fai: Path, chrom: str) -> int:
    if not ref_fai.exists():
        raise FileNotFoundError(f"reference .fai not found: {ref_fai}")
    with open(ref_fai) as fh:
        for line in fh:
            parts = line.split("\t")
            if parts and parts[0] == chrom:
                return int(parts[1])
    raise KeyError(f"chrom {chrom} not in {ref_fai}")


# =============================================================================
# Pydantic request models
# =============================================================================

_CHROM_RE = re.compile(r"^[A-Za-z0-9_.\-]+$")
_GROUP_NAME_RE = re.compile(r"^[A-Za-z0-9_]+$")


def _validate_groups_shape(groups: Dict[str, List[str]], min_n: int) -> None:
    if not groups or len(groups) < 1:
        raise ValueError("groups must contain at least one named group")
    if len(groups) > 10:
        raise ValueError("at most 10 groups (engine F MAX_GROUPS=10)")
    for name, members in groups.items():
        if not _GROUP_NAME_RE.match(name):
            raise ValueError(f"invalid group name '{name}' "
                             "(use [A-Za-z0-9_])")
        if not isinstance(members, list) or not members:
            raise ValueError(f"group '{name}' must be non-empty list")
        if len(members) < min_n:
            raise ValueError(
                f"group '{name}' has n={len(members)} < min_group_n={min_n}")
        if len(set(members)) != len(members):
            raise ValueError(f"group '{name}' has duplicate member ids")


class Region(BaseModel):
    start_bp: int = Field(..., ge=0)
    end_bp: int = Field(..., ge=1)

    @field_validator("end_bp")
    @classmethod
    def _end_after_start(cls, v: int, info: Any) -> int:
        start = (info.data or {}).get("start_bp")
        if start is not None and v <= start:
            raise ValueError("end_bp must be > start_bp")
        return v


class PopstatsGroupwiseReq(BaseModel):
    chrom: str
    region: Optional[Region] = None
    groups: Dict[str, List[str]]
    metrics: List[str] = Field(default_factory=lambda: ["fst", "dxy", "theta_pi"])
    win_bp: int = 50000
    step_bp: int = 10000
    win_type: int = 2
    downsample: int = 1
    ncores: int = 4

    @field_validator("chrom")
    @classmethod
    def _chrom_ok(cls, v: str) -> str:
        if not _CHROM_RE.match(v):
            raise ValueError("chrom contains illegal characters")
        return v


class HobsGroupwiseReq(BaseModel):
    chrom: str
    region: Optional[Region] = None
    groups: Dict[str, List[str]]
    scales: List[str] = Field(default_factory=lambda: ["10kb"])

    @field_validator("chrom")
    @classmethod
    def _chrom_ok(cls, v: str) -> str:
        if not _CHROM_RE.match(v):
            raise ValueError("chrom contains illegal characters")
        return v


class AncestryGroupwiseQReq(BaseModel):
    chrom: str
    region: Optional[Region] = None
    groups: Dict[str, List[str]]
    K: int = 8
    scale: str = "dense"

    @field_validator("chrom")
    @classmethod
    def _chrom_ok(cls, v: str) -> str:
        if not _CHROM_RE.match(v):
            raise ValueError("chrom contains illegal characters")
        return v


class ShelfLDReq(BaseModel):
    chrom: str
    shelf: Region
    n_bins: int = 20
    invgt_assignments: Dict[str, str] = Field(default_factory=dict)


# =============================================================================
# Engine F (region_popstats) wrapper
# =============================================================================
# Workflow:
#   1. Validate group members exist in sample_list (reject unknowns)
#   2. Write per-group sample-id files under a per-request scratch dir
#   3. Build --groups "name1:path1,name2:path2,..." spec
#   4. Run binary with --range if region provided
#   5. Parse output TSV (skip the leading "# ..." comment line)
#   6. Convert to JSON-serializable shape


@dataclasses.dataclass
class EngineFOutcome:
    rows: List[Dict[str, Any]]
    columns: List[str]
    header_comment: str
    stderr_tail: str


def _write_group_files(scratch: Path, groups: Dict[str, List[str]]) -> Tuple[str, Dict[str, Path]]:
    """Write per-group .txt files, return (--groups spec string, paths dict)."""
    parts: List[str] = []
    paths: Dict[str, Path] = {}
    for name, members in sorted(groups.items()):
        gp = scratch / f"{name}.txt"
        with open(gp, "w") as fh:
            fh.write("\n".join(members) + "\n")
        parts.append(f"{name}:{gp}")
        paths[name] = gp
    return ",".join(parts), paths


def _parse_popstats_tsv(out_path: Path) -> EngineFOutcome:
    """Parse region_popstats output. First line is `# region_popstats ...`,
    second is the column header, then data rows."""
    rows: List[Dict[str, Any]] = []
    columns: List[str] = []
    header_comment = ""
    with open(out_path) as fh:
        first = fh.readline()
        if not first.startswith("#"):
            # Engine omitted the comment for some reason; treat first as header
            columns = first.strip().split("\t")
        else:
            header_comment = first[1:].strip()
            header_line = fh.readline()
            if not header_line:
                return EngineFOutcome([], [], header_comment, "")
            columns = header_line.strip().split("\t")
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            vals = line.split("\t")
            if len(vals) != len(columns):
                continue  # skip malformed
            row: Dict[str, Any] = {}
            for k, v in zip(columns, vals):
                # Numeric coercion: int for the count-y columns, float otherwise.
                if k in ("start", "end", "n_sites", "n_sites_used", "S"):
                    try:
                        row[k] = int(v)
                    except ValueError:
                        row[k] = None
                elif k in ("window_id", "chrom"):
                    row[k] = v
                else:
                    try:
                        f = float(v)
                        row[k] = f if np.isfinite(f) else None
                    except ValueError:
                        row[k] = None
            rows.append(row)
    return EngineFOutcome(rows, columns, header_comment, "")


async def run_region_popstats(
    binary: Path,
    beagle: Path,
    sample_list: Path,
    chrom: str,
    fixed_win: Tuple[int, int],
    win_type: int,
    downsample: int,
    ncores: int,
    groups: Dict[str, List[str]],
    region: Optional[Tuple[int, int]],
    scratch: Path,
) -> EngineFOutcome:
    """Run engine F asynchronously. Returns parsed outcome."""
    scratch.mkdir(parents=True, exist_ok=True)
    spec, _paths = _write_group_files(scratch, groups)
    out_path = scratch / "popstats.tsv"

    cmd: List[str] = [
        str(binary),
        "--beagle", str(beagle),
        "--sample_list", str(sample_list),
        "--groups", spec,
        "--fixed_win", f"{fixed_win[0]}:{fixed_win[1]}",
        "--type", str(int(win_type)),
        "--downsample", str(int(downsample)),
        "--chr", chrom,
        "--out", str(out_path),
        "--ncores", str(int(ncores)),
    ]
    if region is not None:
        cmd += ["--range", f"{region[0]}:{region[1]}"]

    log.info("engine F: %s", " ".join(cmd))
    proc = await asyncio.create_subprocess_exec(
        *cmd,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
    )
    _stdout, stderr = await proc.communicate()
    stderr_text = stderr.decode("utf-8", errors="replace")
    if proc.returncode != 0:
        tail = stderr_text[-4096:]
        raise RuntimeError(
            f"region_popstats exited {proc.returncode}.\n"
            f"--- stderr (last 4 KB) ---\n{tail}"
        )
    if not out_path.exists():
        tail = stderr_text[-4096:]
        raise RuntimeError(
            f"region_popstats produced no output file.\n"
            f"--- stderr ---\n{tail}"
        )
    outcome = _parse_popstats_tsv(out_path)
    outcome.stderr_tail = stderr_text[-2048:]
    return outcome


# =============================================================================
# Engine H (Q07b ANGSD + Q07c hobs_windower) wrapper
# =============================================================================

def _resolve_bam(bam_dir: Path, sample_id: str, suffix: str) -> Optional[Path]:
    """Match Q07b's resolution: try <id>.<suffix>, then <id>.bam, then <id>.markdup.bam."""
    for cand in (
        bam_dir / f"{sample_id}.{suffix}",
        bam_dir / f"{sample_id}.bam",
        bam_dir / f"{sample_id}.markdup.bam",
    ):
        if cand.exists():
            return cand
    return None


async def run_angsd_hwe(
    binary: Path,
    bam_list_path: Path,
    ref: Path,
    chrom: str,
    out_prefix: Path,
    threads: int,
    angsd_params: Dict[str, Any],
) -> Tuple[Path, str]:
    """Run patched ANGSD with -doHWE 1 on a single group's BAM list. Returns
    (path-to-.hwe.gz, stderr-text)."""
    cmd: List[str] = [
        str(binary),
        "-bam", str(bam_list_path),
        "-ref", str(ref),
        "-out", str(out_prefix),
        "-GL", "1",
        "-doMajorMinor", str(angsd_params.get("major_minor", 1)),
        "-doMaf", "1",
        "-SNP_pval", str(angsd_params.get("snp_pval", "1e-6")),
        "-doHWE", "1",
        "-maxHetFreq", str(angsd_params.get("max_het_freq", 1.0)),
        "-minMaf", str(angsd_params.get("min_maf", 0.05)),
        "-minMapQ", str(angsd_params.get("min_mapq", 20)),
        "-minQ", str(angsd_params.get("min_q", 20)),
        "-r", f"{chrom}:",
        "-nThreads", str(int(threads)),
    ]
    log.info("ANGSD: %s", " ".join(cmd))
    proc = await asyncio.create_subprocess_exec(
        *cmd, stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.PIPE)
    _stdout, stderr = await proc.communicate()
    stderr_text = stderr.decode("utf-8", errors="replace")
    if proc.returncode != 0:
        raise RuntimeError(
            f"angsd exited {proc.returncode}\n--- stderr (last 4 KB) ---\n"
            f"{stderr_text[-4096:]}")
    out_hwe = out_prefix.parent / f"{out_prefix.name}.hwe.gz"
    if not out_hwe.exists():
        raise RuntimeError(
            f"angsd produced no .hwe.gz at {out_hwe}\n"
            f"--- stderr ---\n{stderr_text[-4096:]}")
    return out_hwe, stderr_text[-2048:]


def _parse_hwe_sites(hwe_gz: Path) -> pd.DataFrame:
    """Parse the patched-ANGSD .hwe.gz output. Columns are
    chr pos major minor hweFreq freq F LRT pval [hetFreq] (10-col variant
    when -maxHetFreq < 1.0 mode, else 9). Mirrors hobs_windower's parser."""
    df = pd.read_csv(hwe_gz, sep="\t", compression="gzip", low_memory=False)
    # Normalize column casing — patched ANGSD writes "Chromo Position Major Minor hweFreq Freq F LRT p-value"
    rename = {c: c.strip() for c in df.columns}
    df = df.rename(columns=rename)
    # Lookup tolerant column names
    candidates_pos = ["Position", "position", "pos"]
    candidates_freq = ["Freq", "freq"]
    candidates_hwe = ["hweFreq", "Hwefreq", "hweFreqfix"]
    candidates_F = ["F", "f"]
    candidates_pval = ["p-value", "pval", "Pvalue", "pvalue"]
    def _pick(cands: List[str]) -> str:
        for c in cands:
            if c in df.columns:
                return c
        raise KeyError(f"none of {cands} in {list(df.columns)}")
    pos_col = _pick(candidates_pos)
    freq_col = _pick(candidates_freq)
    hwe_col = _pick(candidates_hwe)
    F_col = _pick(candidates_F)
    out = pd.DataFrame({
        "pos":     df[pos_col].astype(int),
        "freq":    df[freq_col].astype(float),
        "hweFreq": df[hwe_col].astype(float),
        "F":       df[F_col].astype(float),
    })
    # Mérot logic: Hexp = 2*p*(1-p), Hobs = Hexp*(1-F)
    p = out["hweFreq"].clip(0.0, 1.0)
    out["Hexp"] = 2.0 * p * (1.0 - p)
    out["Hobs"] = (out["Hexp"] * (1.0 - out["F"])).clip(0.0, 1.0)
    return out


_SCALE_RE = re.compile(r"^([A-Za-z0-9]+):(\d+):(\d+)$")
_NAMED_SCALES: Dict[str, str] = {
    "5kb":   "5kb:5000:1000",
    "10kb":  "10kb:10000:2000",
    "50kb":  "50kb:50000:10000",
    "100kb": "100kb:100000:20000",
    "250kb": "250kb:250000:50000",
    "500kb": "500kb:500000:100000",
    "1Mb":   "1Mb:1000000:200000",
}


def _normalize_scales(scales: List[str], default_specs: List[str]) -> List[Tuple[str, int, int]]:
    """Accept either named scales ('10kb') or full specs ('10kb:10000:2000')."""
    spec_map: Dict[str, Tuple[str, int, int]] = {}
    for spec in default_specs:
        m = _SCALE_RE.match(spec)
        if m:
            spec_map[m.group(1)] = (m.group(1), int(m.group(2)), int(m.group(3)))
    out: List[Tuple[str, int, int]] = []
    for s in scales:
        m = _SCALE_RE.match(s)
        if m:
            out.append((m.group(1), int(m.group(2)), int(m.group(3))))
            continue
        if s in spec_map:
            out.append(spec_map[s])
            continue
        if s in _NAMED_SCALES:
            mm = _SCALE_RE.match(_NAMED_SCALES[s])
            assert mm
            out.append((mm.group(1), int(mm.group(2)), int(mm.group(3))))
            continue
        raise ValueError(f"unknown scale '{s}' (try '10kb' or 'label:win:step')")
    return out


def _window_aggregate(sites: pd.DataFrame, win_bp: int, step_bp: int,
                      chrom_size: int) -> Dict[str, List[float]]:
    """Vectorized sliding-window mean of Hobs and Hexp.

    Mirrors hobs_windower's window grid: windows from 1 to chrom_size in
    `step_bp` steps; window endpoints clipped to chrom_size. Empty windows
    are NaN (the C binary skips them; we keep NaN-aligned vectors so all
    groups + scales share the same x-axis on the atlas side)."""
    n_sites = len(sites)
    if n_sites == 0:
        # empty result: still emit empty arrays
        return {"start_bp": [], "end_bp": [], "center_bp": [],
                "n_sites": [], "mean_Hobs": [], "mean_Hexp": [], "HoverE": []}
    pos = sites["pos"].values  # already sorted (ANGSD writes in order)
    Hobs = sites["Hobs"].values
    Hexp = sites["Hexp"].values

    starts: List[int] = []
    centers: List[int] = []
    ends: List[int] = []
    n_per: List[int] = []
    mh: List[float] = []
    me: List[float] = []
    he_ratio: List[float] = []

    # Two-pointer window walk
    ws = 1
    si = 0
    sj = 0
    while ws <= chrom_size:
        we = ws + win_bp - 1
        if we > chrom_size:
            we = chrom_size
        # advance si to first pos >= ws
        while si < n_sites and pos[si] < ws:
            si += 1
        # advance sj to first pos > we
        if sj < si:
            sj = si
        while sj < n_sites and pos[sj] <= we:
            sj += 1
        n = sj - si
        starts.append(int(ws))
        ends.append(int(we))
        centers.append(int((ws + we) // 2))
        n_per.append(int(n))
        if n > 0:
            mh.append(float(Hobs[si:sj].mean()))
            me.append(float(Hexp[si:sj].mean()))
            denom = me[-1]
            he_ratio.append(float(mh[-1] / denom) if denom > 1e-12 else float("nan"))
        else:
            mh.append(float("nan"))
            me.append(float("nan"))
            he_ratio.append(float("nan"))
        ws += step_bp
    # Replace NaN with None for clean JSON
    def _nn(xs: List[float]) -> List[Optional[float]]:
        return [None if (x != x) else x for x in xs]  # NaN check
    return {
        "start_bp": starts, "end_bp": ends, "center_bp": centers,
        "n_sites": n_per,
        "mean_Hobs": _nn(mh), "mean_Hexp": _nn(me), "HoverE": _nn(he_ratio),
    }


# =============================================================================
# Ancestry groupwise Q (no engine — pure pandas)
# =============================================================================

def _resolve_local_q_samples(local_q_dir: Path, chrom: str, scale: str, K: int) -> Path:
    """Mirrors Q06's resolution chain. Tries scaled+K layout, then canonical-K
    symlinks, then flat fallback."""
    kk = f"K{int(K):02d}"
    candidates = [
        local_q_dir / f"scale_{scale}" / kk / f"{chrom}.local_Q_samples.tsv.gz",
        local_q_dir / f"scale_{scale}" / f"{chrom}.local_Q_samples.tsv.gz",
        local_q_dir / scale / kk / f"{chrom}.local_Q_samples.tsv.gz",
        local_q_dir / kk / f"{chrom}.local_Q_samples.tsv.gz",
        local_q_dir / f"{chrom}.local_Q_samples.tsv.gz",
    ]
    for c in candidates:
        if c.exists():
            return c
    raise FileNotFoundError(
        f"local_Q_samples not found for chrom={chrom} scale={scale} K={K}\n"
        f"Tried:\n  " + "\n  ".join(str(c) for c in candidates))


def _per_group_mean_q(samples_path: Path, groups: Dict[str, List[str]],
                      region: Optional[Region]) -> Dict[str, Any]:
    """Read local_Q_samples.tsv.gz, group rows by sample-membership, compute
    per-window mean Q-vector per group.

    Expected schema (canonical Engine B): long-format rows with at least
        sample, chrom, window_mid_bp, Q1, Q2, ..., QK
    OR wide: window_mid_bp + one column per sample, values = K-vector.
    We sniff the columns and dispatch.
    """
    df = pd.read_csv(samples_path, sep="\t", compression="gzip", low_memory=False)
    cols = list(df.columns)
    q_cols = [c for c in cols if c.startswith("Q") and c[1:].isdigit()]
    if "sample" in cols and q_cols and "window_mid_bp" in cols:
        # Long format
        long_df = df
        if region is not None:
            mask = (long_df["window_mid_bp"] >= region.start_bp) & \
                   (long_df["window_mid_bp"] <= region.end_bp)
            long_df = long_df.loc[mask]
        out: Dict[str, Any] = {"K": len(q_cols), "groups": {}, "schema": "long"}
        # x-axis (windows) — derive from any single sample's row order; safer:
        # take sorted unique window_mid_bp.
        windows = np.sort(long_df["window_mid_bp"].unique())
        out["window_mid_bp"] = windows.astype(int).tolist()
        for gname, members in groups.items():
            sub = long_df[long_df["sample"].isin(members)]
            if sub.empty:
                out["groups"][gname] = {
                    "n": 0,
                    "Q_mean": [[None] * len(q_cols) for _ in windows],
                }
                continue
            grouped = sub.groupby("window_mid_bp")[q_cols].mean()
            grouped = grouped.reindex(windows)
            mat = grouped.values  # shape (n_win, K)
            # Replace NaN with None for JSON
            qmean = [
                [None if (v != v) else float(v) for v in row]
                for row in mat
            ]
            out["groups"][gname] = {"n": int(sub["sample"].nunique()), "Q_mean": qmean}
        return out
    # Wide format (window_mid_bp + per-sample columns of K-vector strings)
    if "window_mid_bp" in cols:
        # Heuristic: every column whose name is a sample id is a Q vector
        sample_cols = [c for c in cols if c not in ("chrom", "window_mid_bp",
                                                    "window_start_bp", "window_end_bp")]
        if region is not None:
            df = df[(df["window_mid_bp"] >= region.start_bp) &
                    (df["window_mid_bp"] <= region.end_bp)]
        windows = df["window_mid_bp"].astype(int).tolist()
        out = {"K": None, "groups": {}, "schema": "wide",
               "window_mid_bp": windows}
        # Each cell is "q1,q2,...,qK"; parse on demand per group
        for gname, members in groups.items():
            present = [s for s in members if s in df.columns]
            if not present:
                out["groups"][gname] = {"n": 0, "Q_mean": [[None]] * len(windows)}
                continue
            stacks: List[np.ndarray] = []
            K_seen: Optional[int] = None
            for sid in present:
                col_vals = df[sid].astype(str).values
                arr_list: List[List[float]] = []
                for v in col_vals:
                    parts = [p.strip() for p in v.split(",") if p.strip()]
                    try:
                        arr = [float(x) for x in parts]
                    except ValueError:
                        arr = []
                    arr_list.append(arr)
                if not arr_list:
                    continue
                K_here = max((len(a) for a in arr_list), default=0)
                K_seen = K_here if K_seen is None else max(K_seen, K_here)
                # Pad short rows
                arr_padded = [a + [float("nan")] * (K_here - len(a)) for a in arr_list]
                stacks.append(np.asarray(arr_padded))
            if not stacks:
                out["groups"][gname] = {"n": 0, "Q_mean": [[None]] * len(windows)}
                continue
            stack = np.stack(stacks, axis=0)  # shape (n_samp, n_win, K)
            mean = np.nanmean(stack, axis=0)  # shape (n_win, K)
            qmean = [[None if (v != v) else float(v) for v in row] for row in mean]
            out["groups"][gname] = {"n": len(present), "Q_mean": qmean}
            out["K"] = K_seen
        return out
    raise ValueError(
        f"local_Q_samples.tsv.gz schema not recognized (cols: {cols[:8]}...)")


# =============================================================================
# App + globals
# =============================================================================
# The server can be launched two ways:
#   1. `python3 popstats_server.py --config <path>` — main() calls _bootstrap
#   2. `uvicorn popstats_server:app` — startup hook reads POPSTATS_CONFIG env
# Both end up calling _bootstrap() exactly once before any request handler runs.

# Loaded at startup
CFG: Dict[str, Any] = {}
ENGINES: Optional[EngineRegistry] = None
CACHE: Optional[DiskCache] = None
SAMPLES: Optional[SampleListIndex] = None
JOBS: JobManager = JobManager()


from contextlib import asynccontextmanager


@asynccontextmanager
async def _lifespan(app: FastAPI):
    # On uvicorn-direct startup, _bootstrap may not have been called yet.
    # Honor POPSTATS_CONFIG=<path> if set so `uvicorn popstats_server:app`
    # works the same as `python popstats_server.py --config <path>`.
    if ENGINES is None:
        cfg_path = os.environ.get("POPSTATS_CONFIG")
        if cfg_path:
            try:
                _bootstrap(Path(cfg_path))
            except Exception as e:
                log.error("startup _bootstrap failed: %s", e)
    yield


app = FastAPI(title="popstats_server", version="1.0", lifespan=_lifespan)


# CORS must be added before the app starts handling requests. We default to
# permissive ('*') and let _bootstrap() rewrite the origin list onto the
# already-installed middleware via its `allow_origins` attribute. Atlas runs
# from file:// (origin "null") or 127.0.0.1, so a tight default is fine.
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],   # overwritten by _bootstrap() if config sets origins
    allow_credentials=False,
    allow_methods=["GET", "POST", "DELETE", "OPTIONS"],
    allow_headers=["*"],
)


def _ensure_ready() -> None:
    """Defensive guard for handlers — raise 503 if startup didn't complete."""
    if ENGINES is None or CACHE is None or SAMPLES is None:
        raise HTTPException(
            503,
            "server not initialized — set POPSTATS_CONFIG env var or run via "
            "`python popstats_server.py --config <path>`")


# =============================================================================
# Health
# =============================================================================

@app.get("/api/health")
async def health() -> Dict[str, Any]:
    _ensure_ready()
    return {
        "ok": True,
        "service": "popstats_server",
        "version": "1.0",
        "engines": ENGINES.all_hashes(),
        "engine_paths": {k: str(v) for k, v in ENGINES._paths.items()},
        "cache": {
            "entries": len(CACHE._index),
            "bytes_used": CACHE.total_bytes(),
            "max_bytes": CACHE.max_bytes,
            "root": str(CACHE.root),
        },
        "data": {
            "beagle_dir": CFG.get("beagle_dir"),
            "sample_list": str(SAMPLES.path),
            "n_samples": len(SAMPLES),
            "local_q_dir": CFG.get("local_q_dir"),
        },
        "limits": {
            "min_group_n": CFG.get("min_group_n", 10),
            "max_groups": 10,
        },
    }


@app.get("/api/cache/keys")
async def cache_keys(prefix: str = "", limit: int = 200) -> Dict[str, Any]:
    _ensure_ready()
    items = CACHE.keys_with_prefix(prefix)[-int(limit):]
    return {"n": len(items), "items": items}


@app.delete("/api/cache/keys/{key}")
async def cache_drop(key: str) -> Dict[str, Any]:
    _ensure_ready()
    ok = CACHE.drop(key)
    return {"ok": ok, "key": key}


@app.get("/api/jobs/{job_id}")
async def get_job(job_id: str, include_result: int = 0) -> JSONResponse:
    # NB: jobs are in-process and don't need _ensure_ready, but keeping the
    # check as a guard so unbootstrapped servers consistently 503.
    _ensure_ready()
    st = JOBS.get(job_id)
    if not st:
        raise HTTPException(404, f"unknown job_id: {job_id}")
    payload = st.to_dict()
    if int(include_result) and st.status == "done" and st.result is not None:
        payload["result"] = st.result
    return JSONResponse(payload)


# =============================================================================
# /api/popstats/groupwise
# =============================================================================

@app.post("/api/popstats/groupwise")
async def popstats_groupwise(req: PopstatsGroupwiseReq) -> Response:
    _ensure_ready()
    min_n = int(CFG.get("min_group_n", 10))
    try:
        _validate_groups_shape(req.groups, min_n)
    except ValueError as e:
        raise HTTPException(400, str(e))

    # Validate members exist in canonical sample list
    all_dropped: Dict[str, List[str]] = {}
    cleaned: Dict[str, List[str]] = {}
    for gname, members in req.groups.items():
        kept, dropped = SAMPLES.filter_known(members)
        if dropped:
            all_dropped[gname] = dropped
        if len(kept) < min_n:
            raise HTTPException(
                400,
                f"group '{gname}' has only {len(kept)} known samples after "
                f"filtering against canonical sample_list (min_group_n={min_n}). "
                f"Dropped (unknown): {dropped[:5]}"
            )
        cleaned[gname] = kept

    region_dict: Optional[Dict[str, int]] = None
    region_tuple: Optional[Tuple[int, int]] = None
    if req.region is not None:
        region_dict = {"start_bp": req.region.start_bp, "end_bp": req.region.end_bp}
        region_tuple = (req.region.start_bp, req.region.end_bp)

    engine_hashes = ENGINES.all_hashes()
    key = popstats_cache_key(
        chrom=req.chrom, region=region_dict, groups=cleaned,
        metrics=req.metrics, win_bp=req.win_bp, step_bp=req.step_bp,
        win_type=req.win_type, downsample=req.downsample,
        engines=engine_hashes,
    )

    cached = CACHE.get_json(key, "popstats")
    if cached is not None:
        cached["_cache"] = "hit"
        cached["_cache_key"] = key
        return JSONResponse(cached)

    # Cache miss → run engine F
    binary = ENGINES.path("region_popstats")
    if not binary.exists():
        raise HTTPException(500, f"engine not compiled: {binary}")
    sample_list = Path(SAMPLES.path)
    beagle_dir = Path(CFG["beagle_dir"])
    try:
        beagle = resolve_beagle(beagle_dir, req.chrom)
    except FileNotFoundError as e:
        raise HTTPException(404, str(e))

    with tempfile.TemporaryDirectory(prefix="popstats_", dir=str(CACHE.root)) as scratch:
        scratch_path = Path(scratch)
        try:
            outcome = await run_region_popstats(
                binary=binary, beagle=beagle, sample_list=sample_list,
                chrom=req.chrom,
                fixed_win=(req.win_bp, req.step_bp),
                win_type=req.win_type,
                downsample=req.downsample,
                ncores=req.ncores,
                groups=cleaned,
                region=region_tuple,
                scratch=scratch_path,
            )
        except RuntimeError as e:
            raise HTTPException(500, f"engine F failed: {e}")

    payload = {
        "kind": "popstats_groupwise.v1",
        "chrom": req.chrom,
        "region": region_dict,
        "groups": {n: len(m) for n, m in cleaned.items()},
        "header_comment": outcome.header_comment,
        "columns": outcome.columns,
        "windows": outcome.rows,
        "n_windows": len(outcome.rows),
        "metrics_requested": req.metrics,
        "engine_hash": engine_hashes.get("region_popstats"),
        "warnings": [
            f"dropped unknown samples in group '{g}': {ids[:10]}"
            for g, ids in all_dropped.items()
        ] if all_dropped else [],
    }
    CACHE.put_json(key, "popstats", payload)
    payload["_cache"] = "miss"
    payload["_cache_key"] = key
    return JSONResponse(payload)


# =============================================================================
# /api/popstats/hobs_groupwise
# =============================================================================
# Two-stage: ANGSD per group (slow if cache miss) → hobs_windower per group
# (fast). On HWE cache miss for any group we run ANGSD synchronously here
# because chat A is single-user. A future enhancement can return 202 + job_id
# and run ANGSD in a background task; the JobManager + endpoint are wired.

async def _ensure_hwe_for_group(
    chrom: str,
    members: List[str],
    angsd_params: Dict[str, Any],
    bam_dir: Path, bam_suffix: str,
    angsd_bin: Path, ref: Path,
    scratch: Path,
) -> Path:
    """Return path to .hwe.gz for this group, running ANGSD if cache misses."""
    angsd_hash = ENGINES.hash("angsd_patched")
    hwe_key = hwe_cache_key(chrom=chrom, members=members,
                            angsd_params=angsd_params, engine_hash=angsd_hash)
    cached_path = CACHE.hwe_path(hwe_key)
    if cached_path.exists() and cached_path.stat().st_size > 0:
        return cached_path

    # Resolve BAMs
    bam_paths: List[Path] = []
    missing: List[str] = []
    for sid in members:
        b = _resolve_bam(bam_dir, sid, bam_suffix)
        if b is None:
            missing.append(sid)
        else:
            bam_paths.append(b)
    if missing:
        raise RuntimeError(
            f"BAMs not found for {len(missing)} samples (e.g. {missing[:5]}) "
            f"under {bam_dir}")
    if not bam_paths:
        raise RuntimeError("no BAMs resolved for group")

    # Write bamlist
    scratch.mkdir(parents=True, exist_ok=True)
    bamlist_path = scratch / f"{hwe_key}.bamlist"
    with open(bamlist_path, "w") as fh:
        fh.write("\n".join(str(p) for p in bam_paths) + "\n")

    out_prefix = scratch / hwe_key
    threads = int(angsd_params.get("threads", 2))
    out_hwe, _stderr = await run_angsd_hwe(
        binary=angsd_bin, bam_list_path=bamlist_path,
        ref=ref, chrom=chrom, out_prefix=out_prefix,
        threads=threads, angsd_params=angsd_params)

    # Move into stable cache
    cached_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(out_hwe), str(cached_path))
    return cached_path


@app.post("/api/popstats/hobs_groupwise")
async def hobs_groupwise(req: HobsGroupwiseReq) -> Response:
    _ensure_ready()
    min_n = int(CFG.get("min_group_n", 10))
    try:
        _validate_groups_shape(req.groups, min_n)
    except ValueError as e:
        raise HTTPException(400, str(e))

    cleaned: Dict[str, List[str]] = {}
    all_dropped: Dict[str, List[str]] = {}
    for gname, members in req.groups.items():
        kept, dropped = SAMPLES.filter_known(members)
        if dropped:
            all_dropped[gname] = dropped
        if len(kept) < min_n:
            raise HTTPException(
                400,
                f"group '{gname}' n={len(kept)} < min_group_n={min_n} "
                f"after filtering")
        cleaned[gname] = kept

    region_dict: Optional[Dict[str, int]] = None
    if req.region is not None:
        region_dict = {"start_bp": req.region.start_bp, "end_bp": req.region.end_bp}

    engine_hashes = ENGINES.all_hashes()
    key = hobs_cache_key(
        chrom=req.chrom, region=region_dict, groups=cleaned,
        scales=req.scales, engines=engine_hashes,
    )
    cached = CACHE.get_json(key, "hobs")
    if cached is not None:
        cached["_cache"] = "hit"
        cached["_cache_key"] = key
        return JSONResponse(cached)

    # Resolve scales
    default_scale_specs: List[str] = CFG.get("hobs_defaults", {}).get(
        "scales", list(_NAMED_SCALES.values()))
    try:
        scales = _normalize_scales(req.scales, default_scale_specs)
    except ValueError as e:
        raise HTTPException(400, str(e))

    angsd_params = dict(CFG.get("angsd_defaults", {}))

    bam_dir = Path(CFG["bam_dir"])
    bam_suffix = CFG.get("bam_suffix", "markdup.bam")
    angsd_bin = ENGINES.path("angsd_patched")
    if not angsd_bin.exists():
        raise HTTPException(500, f"angsd_patched not found: {angsd_bin}")
    ref = Path(CFG["ref"])
    ref_fai = Path(CFG["ref_fai"])
    try:
        chrom_size = resolve_chrom_size(ref_fai, req.chrom)
    except (FileNotFoundError, KeyError) as e:
        raise HTTPException(404, str(e))

    # Run all groups' ANGSD in parallel (matching Q07b's 3-way parallelism).
    with tempfile.TemporaryDirectory(prefix="hobs_", dir=str(CACHE.root)) as scratch:
        scratch_path = Path(scratch)

        async def _one(gname: str, members: List[str]) -> Tuple[str, Path]:
            hwe = await _ensure_hwe_for_group(
                chrom=req.chrom, members=members, angsd_params=angsd_params,
                bam_dir=bam_dir, bam_suffix=bam_suffix,
                angsd_bin=angsd_bin, ref=ref, scratch=scratch_path)
            return gname, hwe

        tasks = [_one(g, m) for g, m in cleaned.items()]
        try:
            results = await asyncio.gather(*tasks)
        except RuntimeError as e:
            raise HTTPException(500, f"angsd failed: {e}")

        hwe_paths: Dict[str, Path] = {g: p for g, p in results}

        # Per-group, per-site -> windowed at each requested scale.
        # We do the windowing in-process (numpy) rather than calling the C
        # hobs_windower binary, because:
        #   (1) we need mean_Hexp alongside mean_Hobs to compute HoverE
        #       for the atlas (the C binary doesn't emit Hexp at window scale)
        #   (2) windowing 5e4 sites at 7 scales is microseconds in numpy
        result_groups: Dict[str, Dict[str, Any]] = {}
        for gname, hwe in hwe_paths.items():
            try:
                sites = _parse_hwe_sites(hwe)
            except (KeyError, ValueError, OSError) as e:
                raise HTTPException(500, f"failed to parse {hwe}: {e}")
            if req.region is not None:
                sites = sites[(sites["pos"] >= req.region.start_bp) &
                              (sites["pos"] <= req.region.end_bp)]
            scales_out: Dict[str, Dict[str, Any]] = {}
            for label, win, step in scales:
                scales_out[label] = _window_aggregate(
                    sites, win_bp=win, step_bp=step, chrom_size=chrom_size)
            result_groups[gname] = {
                "n": len(cleaned[gname]),
                "n_sites": int(len(sites)),
                "scales": scales_out,
            }

    payload = {
        "kind": "hobs_groupwise.v1",
        "chrom": req.chrom,
        "region": region_dict,
        "groups": result_groups,
        "scale_labels": [s[0] for s in scales],
        "engine_hash_hobs_windower": engine_hashes.get("hobs_windower"),
        "engine_hash_angsd": engine_hashes.get("angsd_patched"),
        "warnings": [
            f"dropped unknown samples in group '{g}': {ids[:10]}"
            for g, ids in all_dropped.items()
        ] if all_dropped else [],
    }
    CACHE.put_json(key, "hobs", payload)
    payload["_cache"] = "miss"
    payload["_cache_key"] = key
    return JSONResponse(payload)


# =============================================================================
# /api/ancestry/groupwise_q
# =============================================================================

@app.post("/api/ancestry/groupwise_q")
async def ancestry_groupwise_q(req: AncestryGroupwiseQReq) -> Response:
    _ensure_ready()
    min_n = int(CFG.get("min_group_n", 10))
    try:
        _validate_groups_shape(req.groups, min_n)
    except ValueError as e:
        raise HTTPException(400, str(e))

    cleaned: Dict[str, List[str]] = {}
    for gname, members in req.groups.items():
        kept, _ = SAMPLES.filter_known(members)
        if len(kept) < min_n:
            raise HTTPException(
                400,
                f"group '{gname}' n={len(kept)} < min_group_n={min_n}")
        cleaned[gname] = kept

    region_dict: Optional[Dict[str, int]] = None
    if req.region is not None:
        region_dict = {"start_bp": req.region.start_bp, "end_bp": req.region.end_bp}

    key = ancestry_q_cache_key(
        chrom=req.chrom, region=region_dict, groups=cleaned,
        K=req.K, scale=req.scale)
    cached = CACHE.get_json(key, "ancestry")
    if cached is not None:
        cached["_cache"] = "hit"
        cached["_cache_key"] = key
        return JSONResponse(cached)

    local_q_dir = Path(CFG["local_q_dir"])
    try:
        samples_path = _resolve_local_q_samples(local_q_dir, req.chrom, req.scale, req.K)
    except FileNotFoundError as e:
        raise HTTPException(404, str(e))

    try:
        agg = _per_group_mean_q(samples_path, cleaned, req.region)
    except (ValueError, KeyError, OSError) as e:
        raise HTTPException(500, f"ancestry aggregation failed: {e}")

    payload = {
        "kind": "ancestry_groupwise_q.v1",
        "chrom": req.chrom,
        "region": region_dict,
        "K": req.K,
        "scale": req.scale,
        "source": str(samples_path),
        **agg,
    }
    CACHE.put_json(key, "ancestry", payload)
    payload["_cache"] = "miss"
    payload["_cache_key"] = key
    return JSONResponse(payload)


# =============================================================================
# /api/shelf_ld_test  (optional; brief says JS port is preferred)
# =============================================================================
# Provided as a server fallback. The actual atlas implementation lives in JS.
# We keep this endpoint stub so existing R/Q09b workflows can still hit the
# server, but the rich diagnostic plot stays atlas-side.

@app.post("/api/shelf_ld_test")
async def shelf_ld_test(req: ShelfLDReq) -> Response:
    return JSONResponse({
        "kind": "shelf_ld_test.v1",
        "chrom": req.chrom,
        "shelf": {"start_bp": req.shelf.start_bp, "end_bp": req.shelf.end_bp},
        "n_bins": req.n_bins,
        "note": "server-side shelf LD test not implemented in chat A; "
                "atlas computes this in JS using dosage chunks. "
                "See atlas page-3 candidate-focus shelf-LD module.",
    }, status_code=501)


# =============================================================================
# /api/ld/split_heatmap  (turn 11a)
# =============================================================================
# Live ngsLD on a candidate region for two sample groups (lower triangle =
# all_samples, upper triangle = common_samples / e.g. HOM1). Returns binned
# r² matrices ready for the atlas-side split-triangle renderer.
#
# Cache: composite — lower and upper triangles are cached separately so
# swapping the upper-triangle group reuses the lower-triangle compute.
# Region capped at 8 Mb (page 3 candidate-zoom is ~3 Mb in practice).

@app.post("/api/ld/split_heatmap")
async def ld_split_heatmap(req: LDSplitReq) -> Response:
    _ensure_ready()

    # Resolve master BEAGLE for this chrom
    beagle_dir = Path(CFG["beagle_dir"])
    try:
        master_beagle = resolve_beagle(beagle_dir, req.chrom)
    except FileNotFoundError as e:
        raise HTTPException(404, str(e))

    # ngsLD binary path. Convention: configured under engines.ngsld.
    try:
        ngsld_bin = ENGINES.path("ngsld")
    except KeyError:
        raise HTTPException(
            500,
            "ngsLD binary not configured. Add 'ngsld: <path>' under engines "
            "in popstats_server_config.yaml. Typical path: "
            "<base>/ngsTools/ngsLD/ngsLD"
        )
    if not ngsld_bin.exists():
        raise HTTPException(500, f"ngsLD binary missing: {ngsld_bin}")
    ngsld_hash = ENGINES.hash("ngsld")

    # Adapter to the cache layer — match the popstats_groupwise idiom
    def cache_get(key: str) -> Optional[Dict[str, Any]]:
        return CACHE.get_json(key, "ld")
    def cache_put(key: str, payload: Dict[str, Any]) -> None:
        CACHE.put_json(key, "ld", payload)

    payload = await handle_split_heatmap(
        req,
        ngsld_bin=ngsld_bin,
        master_beagle=master_beagle,
        sample_filter=SAMPLES.filter_known,
        ngsld_hash=ngsld_hash,
        cache_get=cache_get,
        cache_put=cache_put,
        cache_root=CACHE.root,
    )
    return JSONResponse(payload)


# =============================================================================
# Startup
# =============================================================================

def _bootstrap(config_path: Path) -> None:
    global CFG, ENGINES, CACHE, SAMPLES
    CFG = load_config(config_path)
    log.info("config loaded from %s", config_path)
    # Update the already-installed CORS middleware's allowed origins.
    cors_origins = CFG.get("cors_origins", ["*"])
    for m in app.user_middleware:
        if m.cls is CORSMiddleware:
            # Starlette stores middleware kwargs on the marker; mutate them
            # so subsequent requests see the new origin list.
            m.kwargs["allow_origins"] = list(cors_origins)
            break
    # Engine registry
    engine_paths: Dict[str, str] = CFG.get("engines", {})
    if not engine_paths:
        raise RuntimeError("config missing 'engines' section")
    ENGINES = EngineRegistry(
        engine_paths,
        ttl_sec=float(CFG.get("cache_health_ttl_sec", 60)),
    )
    # Cache
    cache_dir = Path(CFG.get("cache_dir", "./popstats_cache"))
    CACHE = DiskCache(cache_dir, int(CFG.get("cache_max_bytes", 5 * 10**10)))
    # Sample list
    SAMPLES = SampleListIndex(Path(CFG["sample_list"]))
    log.info("startup complete: %d samples, %d engines, %d cached entries",
             len(SAMPLES), len(engine_paths), len(CACHE._index))


# =============================================================================
# Entry point
# =============================================================================

def main() -> None:
    import argparse
    ap = argparse.ArgumentParser(description="popstats_server")
    ap.add_argument("--config", required=True, type=Path,
                    help="path to popstats_server.config.yaml")
    ap.add_argument("--host", default=None,
                    help="override config bind.host")
    ap.add_argument("--port", default=None, type=int,
                    help="override config bind.port")
    ap.add_argument("--reload", action="store_true",
                    help="auto-reload on code change (dev only)")
    args = ap.parse_args()

    _bootstrap(args.config)
    bind = CFG.get("bind", {})
    host = args.host or bind.get("host", "127.0.0.1")
    port = int(args.port or bind.get("port", 8765))
    log.info("listening on http://%s:%d", host, port)

    import uvicorn
    uvicorn.run("popstats_server:app", host=host, port=port,
                reload=args.reload, log_level="info")


if __name__ == "__main__":
    main()
