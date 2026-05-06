#!/usr/bin/env python3
# =============================================================================
# common.py — shared utilities for phase_8_comparative_breakpoint_fragility
# =============================================================================
"""
Tiny, dependency-light helpers used across the pipeline:

  - logger setup (writes to output/logs/<step>__<timestamp>.log)
  - run-manifest sidecar writer (.run.json)
  - chromosome-name normalisation (Gar / Mac / EDTA / RM idiosyncrasies)
  - 1-based GFF -> 0-based BED conversion
  - safe TSV reader with header validation

Kept stdlib-only deliberately so this can run on any LANTA node without
loading a python module.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import logging
import os
import re
import sys
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Iterable, Iterator


# -----------------------------------------------------------------------------
# logging
# -----------------------------------------------------------------------------
def get_logger(step_name: str, log_dir: Path | str) -> logging.Logger:
    """Logger that writes to both stderr and output/logs/<step>__<ts>.log."""
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    ts = time.strftime("%Y%m%dT%H%M%S")
    log_path = log_dir / f"{step_name}__{ts}.log"

    logger = logging.getLogger(step_name)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    fmt = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        "%Y-%m-%d %H:%M:%S",
    )
    fh = logging.FileHandler(log_path)
    fh.setFormatter(fmt)
    sh = logging.StreamHandler(sys.stderr)
    sh.setFormatter(fmt)
    logger.addHandler(fh)
    logger.addHandler(sh)
    logger.info("log file: %s", log_path)
    return logger


# -----------------------------------------------------------------------------
# run manifest sidecar
# -----------------------------------------------------------------------------
def write_run_manifest(out_path: Path | str, **fields) -> Path:
    """Write a <out>.run.json beside any output, capturing inputs/params/counts."""
    out_path = Path(out_path)
    side = out_path.with_suffix(out_path.suffix + ".run.json")
    payload = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "host": os.uname().nodename,
        **fields,
    }
    side.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str))
    return side


def file_sha256(path: Path | str, max_bytes: int = 64 * 1024 * 1024) -> str:
    """Hash the first max_bytes for cheap provenance — full file would be slow."""
    h = hashlib.sha256()
    path = Path(path)
    with path.open("rb") as fh:
        h.update(fh.read(max_bytes))
    return h.hexdigest()[:16]


# -----------------------------------------------------------------------------
# chromosome-name normalisation
# -----------------------------------------------------------------------------
# Project canon (matches reference fClaHyb_Gar_LG.fa / Quentin's atlas):
#   focal Gar  -> "C_gar_LG##"
#   Mac        -> "C_mac_LG##"  (placeholder; adjust when Mac assembly is canonised)
# We accept many inputs and standardise. Unknown patterns are returned as-is
# with a log warning at the call site.
_CHR_PATTERNS = [
    # Already canonical
    (re.compile(r"^(C_gar_LG\d+)$"), r"\1"),
    (re.compile(r"^(C_mac_LG\d+)$"), r"\1"),
    # Common variants
    (re.compile(r"^Gar[_-]?(LG)?(\d+)$", re.I), r"C_gar_LG\2"),
    (re.compile(r"^Mac[_-]?(LG)?(\d+)$", re.I), r"C_mac_LG\2"),
    # NCBI / RefSeq style with a strip-prefix hint we can't resolve here
    # We keep it as-is and let the caller log.
]


def normalize_chrom(chrom: str, species_id: str | None = None) -> str:
    """Best-effort canonicalise a chromosome name. Unknown patterns pass through."""
    if chrom is None:
        return chrom
    c = chrom.strip()
    for pat, rep in _CHR_PATTERNS:
        m = pat.match(c)
        if m:
            return pat.sub(rep, c)
    return c  # unchanged; the caller can decide whether to log


# -----------------------------------------------------------------------------
# GFF/BED coordinate conversion
# -----------------------------------------------------------------------------
def gff_to_bed_coords(start_1based: int, end_1based_inclusive: int) -> tuple[int, int]:
    """GFF3 is 1-based inclusive; BED is 0-based half-open."""
    return start_1based - 1, end_1based_inclusive


def open_maybe_gz(path: Path | str, mode: str = "rt"):
    """Transparent handle for plain or .gz."""
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(p, mode)
    return open(p, mode)


# -----------------------------------------------------------------------------
# safe TSV reader
# -----------------------------------------------------------------------------
def read_tsv(path: Path | str, required_cols: Iterable[str] | None = None) -> list[dict]:
    """Read a TSV with header. Validates required columns."""
    rows: list[dict] = []
    with open_maybe_gz(path, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if required_cols:
            missing = [c for c in required_cols if c not in (reader.fieldnames or [])]
            if missing:
                raise ValueError(f"{path}: missing required columns: {missing}")
        for row in reader:
            rows.append(row)
    return rows


def write_tsv(path: Path | str, rows: list[dict], cols: list[str]) -> None:
    """Write a TSV with explicit column order. Missing keys -> empty string."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gz(path, "wt") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})


# -----------------------------------------------------------------------------
# canonical TE class buckets
# -----------------------------------------------------------------------------
TE_CLASS_BUCKETS = [
    "DNA", "LTR", "LINE", "SINE", "Helitron", "TIR", "MITE",
    "satellite", "low_complexity", "unknown",
]


def bucket_te_class(raw_class: str | None, raw_family: str | None = None) -> str:
    """Map a messy class/family string to one of TE_CLASS_BUCKETS."""
    s = (raw_class or "") + "/" + (raw_family or "")
    s_low = s.lower()

    # ordering matters: more specific first
    if any(k in s_low for k in ("ltr", "gypsy", "copia")):
        return "LTR"
    if "helitron" in s_low or "rc/" in s_low or "rc-" in s_low:
        return "Helitron"
    if "mite" in s_low:
        return "MITE"
    if "tir" in s_low or any(k in s_low for k in ("cacta", "mutator", "hat", "tc1", "mariner", "harbinger", "piggybac")):
        return "TIR"
    if "line" in s_low or "l1" in s_low or "rex" in s_low:
        return "LINE"
    if "sine" in s_low or "trna" in s_low or "5s" in s_low or "7sl" in s_low:
        return "SINE"
    if "dna" in s_low and "transposon" in s_low:
        return "DNA"
    if "satellite" in s_low or "satel" in s_low:
        return "satellite"
    if "low_complex" in s_low or "simple_repeat" in s_low or "tandem" in s_low:
        return "low_complexity"
    return "unknown"
