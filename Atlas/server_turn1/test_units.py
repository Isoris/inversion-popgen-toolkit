#!/usr/bin/env python3
# =============================================================================
# test_units.py — unit tests for popstats_server helpers
# =============================================================================
# Pure-Python tests of the things that aren't covered by the curl smoke tests:
#   - cache key canonicalization (group/member order independence)
#   - region_popstats TSV parser (header line, numeric coercion, NaN handling)
#   - HoverE windowing math (Mérot signature: HoverE=1 at HWE, 0 in homozygotes,
#     2 in heterozygotes)
#   - scale spec parsing
#
# Run:
#   python3 test_units.py
# =============================================================================
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

import popstats_server as ps


def assert_eq(a, b, label=""):
    if a != b:
        print(f"FAIL {label}: {a!r} != {b!r}")
        sys.exit(1)
    print(f"  ok  {label}")


def assert_close(a, b, label="", tol=1e-9):
    if not np.isclose(a, b, atol=tol):
        print(f"FAIL {label}: {a} ≉ {b}")
        sys.exit(1)
    print(f"  ok  {label}")


def test_canonical_groups():
    print("test_canonical_groups")
    a = {"B": ["s2", "s1"], "A": ["z", "a"]}
    b = {"A": ["a", "z"], "B": ["s1", "s2"]}
    assert_eq(ps._canonical_groups(a), ps._canonical_groups(b),
              "order-independent canonicalization")


def test_cache_key_stability():
    print("test_cache_key_stability")
    e = {"region_popstats": "abc"}
    k1 = ps.popstats_cache_key("LG28", None,
                               {"A": ["s1", "s2"], "B": ["s3", "s4"]},
                               ["fst"], 50000, 10000, 2, 1, e)
    k2 = ps.popstats_cache_key("LG28", None,
                               {"B": ["s4", "s3"], "A": ["s2", "s1"]},
                               ["fst"], 50000, 10000, 2, 1, e)
    assert_eq(k1, k2, "popstats key is order-invariant")

    # Different engine hash → different key
    e2 = {"region_popstats": "xyz"}
    k3 = ps.popstats_cache_key("LG28", None,
                               {"A": ["s1", "s2"], "B": ["s3", "s4"]},
                               ["fst"], 50000, 10000, 2, 1, e2)
    if k1 == k3:
        print(f"FAIL: engine hash change should change cache key but {k1} == {k3}")
        sys.exit(1)
    print("  ok  engine hash invalidates cache")


def test_popstats_parser():
    print("test_popstats_parser")
    fake = (
        "# region_popstats downsample=1 type=2\n"
        "window_id\tchrom\tstart\tend\tn_sites\tn_sites_used\tS\t"
        "theta_pi\ttheta_w\ttajD\thet\ttheta_pi_HOM1\tFst_HOM1_HET\n"
        "w1\tLG28\t10000\t60000\t450\t450\t320\t0.012\t0.011\t-1.5\t0.012\t0.010\t0.31\n"
        "w2\tLG28\t60000\t110000\t470\t470\t330\t0.013\t0.012\t-1.4\t0.013\t0.011\t0.30\n"
    )
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as f:
        f.write(fake); p = Path(f.name)
    try:
        out = ps._parse_popstats_tsv(p)
        assert_eq(len(out.rows), 2, "2 rows parsed")
        assert_eq(out.rows[0]["window_id"], "w1", "window_id is string")
        assert_eq(out.rows[0]["start"], 10000, "start is int")
        assert_close(out.rows[0]["theta_pi"], 0.012, "theta_pi is float")
        assert_close(out.rows[0]["Fst_HOM1_HET"], 0.31, "per-pair Fst parsed")
    finally:
        p.unlink()


def test_window_aggregate_hwe():
    print("test_window_aggregate Mérot signatures")
    # 100 sites at 1kb spacing, p=0.5, F=0 → Hexp=0.5, Hobs=0.5, HoverE=1
    sites_hwe = pd.DataFrame({
        "pos": [1000 * i for i in range(1, 101)],
        "Hexp": [0.5] * 100, "Hobs": [0.5] * 100, "F": [0.0] * 100,
        "freq": [0.5] * 100, "hweFreq": [0.5] * 100,
    })
    agg = ps._window_aggregate(sites_hwe, win_bp=10000, step_bp=10000, chrom_size=100000)
    he_clean = [v for v in agg["HoverE"] if v is not None]
    assert_close(np.mean(he_clean), 1.0, "HWE windows: HoverE → 1.0")

    # F=1 in HOMOZYGOTE: Hobs = Hexp*(1-1) = 0 → HoverE = 0
    sites_hom = sites_hwe.copy()
    sites_hom["F"] = 1.0
    sites_hom["Hobs"] = 0.0
    agg2 = ps._window_aggregate(sites_hom, win_bp=10000, step_bp=10000, chrom_size=100000)
    he_clean2 = [v for v in agg2["HoverE"] if v is not None]
    assert_close(np.mean(he_clean2), 0.0, "HOMOZYGOTE windows: HoverE → 0")

    # F = -1 (everyone het): Hobs = Hexp*(1+1) = 2*Hexp → HoverE = 2
    sites_het = sites_hwe.copy()
    sites_het["F"] = -1.0
    sites_het["Hobs"] = 1.0  # per Mérot: Hexp*(1-(-1)) = 2*Hexp = 1.0 (clamp)
    agg3 = ps._window_aggregate(sites_het, win_bp=10000, step_bp=10000, chrom_size=100000)
    he_clean3 = [v for v in agg3["HoverE"] if v is not None]
    assert_close(np.mean(he_clean3), 2.0, "HETEROKARYOTYPE windows: HoverE → 2.0")


def test_scale_normalize():
    print("test_scale_normalize")
    out = ps._normalize_scales(["10kb", "5kb:5000:1000"], ["10kb:10000:2000", "5kb:5000:1000"])
    assert_eq(out[0], ("10kb", 10000, 2000), "named scale → triple")
    assert_eq(out[1], ("5kb", 5000, 1000), "explicit spec passes through")
    try:
        ps._normalize_scales(["bogus"], ["10kb:10000:2000"])
        print("FAIL: should reject bogus scale")
        sys.exit(1)
    except ValueError:
        print("  ok  rejects unknown scale")


def test_validate_groups():
    print("test_validate_groups")
    # Min n
    try:
        ps._validate_groups_shape({"A": ["s1"]}, min_n=10)
        print("FAIL: should reject too-small group")
        sys.exit(1)
    except ValueError:
        print("  ok  rejects too-small group")
    # Bad name
    try:
        ps._validate_groups_shape({"bad name": ["s1"] * 20}, min_n=10)
        print("FAIL: should reject bad name")
        sys.exit(1)
    except ValueError:
        print("  ok  rejects spaces in group name")
    # Duplicate members
    try:
        ps._validate_groups_shape({"A": ["s1", "s1", "s1"] + ["s2"] * 12}, min_n=3)
        print("FAIL: should reject duplicate members")
        sys.exit(1)
    except ValueError:
        print("  ok  rejects duplicate members")
    # Happy path
    ps._validate_groups_shape({"A": [f"s{i}" for i in range(20)],
                                "B": [f"t{i}" for i in range(20)]}, min_n=10)
    print("  ok  happy-path validation")


def test_disk_cache_lru():
    print("test_disk_cache_lru")
    with tempfile.TemporaryDirectory() as td:
        # very small cap → forces eviction
        c = ps.DiskCache(Path(td), max_bytes=200)
        for i in range(10):
            c.put_json(f"k{i:02d}", "popstats", {"i": i, "data": list(range(10))})
        # Should have evicted older keys
        keys = [it["hash"] for it in c.keys_with_prefix()]
        assert_eq(len(keys) > 0 and len(keys) < 10, True,
                  "LRU evicted some entries")
        # Most recent (k09) must still be present
        assert_eq("k09" in keys, True, "newest entry survives eviction")


if __name__ == "__main__":
    test_canonical_groups()
    test_cache_key_stability()
    test_popstats_parser()
    test_window_aggregate_hwe()
    test_scale_normalize()
    test_validate_groups()
    test_disk_cache_lru()
    print("\nALL UNIT TESTS PASSED")
