#!/usr/bin/env python3
# Synthetic end-to-end test: build a tiny PAF + tiny TEfull JSONs, run
# STEP_CS01_extract_breakpoints.py, and verify the output JSON contains
# the breakpoints we engineered. Catches regressions without needing real
# wfmash data.

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path


def main():
    root = Path(__file__).parent
    script = root / "STEP_CS01_extract_breakpoints.py"
    if not script.exists():
        print(f"FATAL: script not found at {script}", file=sys.stderr)
        sys.exit(1)

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        # ---------------------------------------------------------------------
        # Synthesize a PAF for a single Gar chromosome (LG28) with three
        # syntenic blocks that look like:
        #   block A:  Gar LG28  0..3M       on Mac LG28  0..3M       (+)
        #   block B:  Gar LG28  3M..5M      on Mac LG28  5M..3M      (-)   <- INVERSION
        #   block C:  Gar LG28  5M..10M     on Mac LG28  5M..10M     (+)
        #
        # Plus a second Gar chrom (LG07) with a translocation:
        #   block A:  Gar LG07  0..5M       on Mac LG07  0..5M       (+)
        #   block B:  Gar LG07  5M..10M     on Mac LG10  20M..25M    (+)   <- TRANSLOCATION
        # ---------------------------------------------------------------------
        paf_path = tmp / "synthetic.paf"
        # PAF columns:
        # qname qlen qstart qend strand tname tlen tstart tend matches blocklen mapq
        rows = [
            # LG28 — three blocks with an inversion in the middle
            ("LG28", 10_000_000,        0, 3_000_000, "+", "LG28", 10_000_000,        0, 3_000_000, 2_900_000, 3_000_000, 60),
            ("LG28", 10_000_000, 3_000_000, 5_000_000, "-", "LG28", 10_000_000, 3_000_000, 5_000_000, 1_900_000, 2_000_000, 60),
            ("LG28", 10_000_000, 5_000_000,10_000_000, "+", "LG28", 10_000_000, 5_000_000,10_000_000, 4_900_000, 5_000_000, 60),
            # LG07 — translocation midway
            ("LG07", 10_000_000,        0, 5_000_000, "+", "LG07", 30_000_000,        0, 5_000_000, 4_900_000, 5_000_000, 60),
            ("LG07", 10_000_000, 5_000_000,10_000_000, "+", "LG10", 30_000_000,20_000_000,25_000_000, 4_900_000, 5_000_000, 60),
        ]
        with open(paf_path, "w") as f:
            for r in rows:
                f.write("\t".join(str(x) for x in r) + "\n")

        # ---------------------------------------------------------------------
        # Synthesize a Cgar TE-density JSON for LG28 with elevated all_TE
        # near the inversion breakpoint (Mb 3-5). This exercises the
        # flanking-repeat annotation path AND the Spalax-note generation.
        # ---------------------------------------------------------------------
        te_dir = tmp / "Cgar_TE"
        te_dir.mkdir()
        n_windows_28 = 20
        # Window centres: 0.25, 0.75, 1.25, ..., 9.75 Mb
        centers_28 = [(i + 0.5) * 0.5 for i in range(n_windows_28)]
        # all_TE: baseline 0.45, peak 0.92 around 3.0-5.0 Mb
        all_te = []
        for c in centers_28:
            if 2.5 <= c <= 5.5:
                all_te.append(0.92 - 0.05 * abs(c - 4.0))   # peak at 4 Mb
            else:
                all_te.append(0.45)
        young_te = [v * 0.6 for v in all_te]
        gar_lg28_json = {
            "version": 2,
            "species": "Clarias gariepinus",
            "binning_source": "scrubber_windows",
            "precomp_chrom": "LG28",
            "n_chromosomes": 1,
            "n_classes": 2,
            "classes": ["all_TE", "young_TE_all"],
            "default_class": "all_TE",
            "loess_span": 0.3,
            "generated_at": "2026-04-30T00:00:00Z",
            "chromosomes": [{
                "chrom": "LG28",
                "n_windows": n_windows_28,
                "window_centers_mb": centers_28,
                "window_start_bp":  [int(c * 1e6 - 250_000) for c in centers_28],
                "window_end_bp":    [int(c * 1e6 + 250_000) for c in centers_28],
                "by_class": {
                    "all_TE":       {"densities": all_te,   "loess": all_te,   "max_density": max(all_te),   "loess_span": 0.3},
                    "young_TE_all": {"densities": young_te, "loess": young_te, "max_density": max(young_te), "loess_span": 0.3},
                },
            }],
        }
        with open(te_dir / "LG28_repeat_density_TEfull.json", "w") as f:
            json.dump(gar_lg28_json, f)

        # LG07 — flat baseline (no enrichment expected at the translocation)
        n_windows_07 = 20
        centers_07 = [(i + 0.5) * 0.5 for i in range(n_windows_07)]
        gar_lg07_json = {
            **gar_lg28_json,
            "precomp_chrom": "LG07",
            "chromosomes": [{
                "chrom": "LG07",
                "n_windows": n_windows_07,
                "window_centers_mb": centers_07,
                "window_start_bp":  [int(c * 1e6 - 250_000) for c in centers_07],
                "window_end_bp":    [int(c * 1e6 + 250_000) for c in centers_07],
                "by_class": {
                    "all_TE":       {"densities": [0.45]*n_windows_07, "loess": [0.45]*n_windows_07, "max_density": 0.45, "loess_span": 0.3},
                    "young_TE_all": {"densities": [0.27]*n_windows_07, "loess": [0.27]*n_windows_07, "max_density": 0.27, "loess_span": 0.3},
                },
            }],
        }
        with open(te_dir / "LG07_repeat_density_TEfull.json", "w") as f:
            json.dump(gar_lg07_json, f)

        # Mac TE dir (sparse — just enough to exercise the path)
        mac_te_dir = tmp / "Cmac_TE"
        mac_te_dir.mkdir()
        for chrom_name in ("LG28", "LG07", "LG10"):
            n = 20
            cc = [(i + 0.5) * 1.5 for i in range(n)]
            block = {
                "version": 2,
                "species": "Clarias macrocephalus",
                "binning_source": "scrubber_windows",
                "precomp_chrom": chrom_name,
                "classes": ["all_TE"],
                "default_class": "all_TE",
                "chromosomes": [{
                    "chrom": chrom_name,
                    "n_windows": n,
                    "window_centers_mb": cc,
                    "window_start_bp":  [int(c * 1e6 - 750_000) for c in cc],
                    "window_end_bp":    [int(c * 1e6 + 750_000) for c in cc],
                    "by_class": {
                        "all_TE": {"densities": [0.40]*n, "loess": [0.40]*n, "max_density": 0.40, "loess_span": 0.3},
                    },
                }],
            }
            with open(mac_te_dir / f"{chrom_name}_repeat_density_TEfull.json", "w") as f:
                json.dump(block, f)

        # Optional candidates TSV — overlap one of the engineered breakpoints
        cand_path = tmp / "candidates.tsv"
        with open(cand_path, "w") as f:
            f.write("candidate_id\tchrom\tstart_bp\tend_bp\n")
            f.write("CAND001\tLG28\t2500000\t5500000\n")   # contains both LG28 breakpoints

        out_path = tmp / "cs_breakpoints_v1.json"
        cmd = [
            sys.executable, str(script),
            "--paf", str(paf_path),
            "--out", str(out_path),
            "--te-dir", str(te_dir),
            "--mac-te-dir", str(mac_te_dir),
            "--candidates", str(cand_path),
            "--min-mapq", "20",
            "--min-block-bp", "1000",
            "--cluster-radius-bp", "100000",
            "--flank-bp", "1000000",
            "--classes", "all_TE,young_TE_all",
        ]
        print("[test] running:", " ".join(cmd))
        proc = subprocess.run(cmd, capture_output=True, text=True)
        print("[test] stderr:")
        print(proc.stderr)
        if proc.returncode != 0:
            print("[test] FAIL: script returned non-zero")
            print(proc.stdout)
            sys.exit(1)

        with open(out_path) as f:
            result = json.load(f)

        # ------------------ Assertions -------------------------------------
        assert result["tool"] == "cross_species_breakpoints_v1", result["tool"]
        assert result["schema_version"] == 2, f"expected schema v2, got {result['schema_version']}"
        assert result["n_breakpoints"] >= 2, f"expected at least 2 breakpoints, got {result['n_breakpoints']}"

        # v2 additions: synteny blocks + chrom lengths
        assert "synteny_blocks" in result, "v2 must emit synteny_blocks"
        assert "n_synteny_blocks" in result, "v2 must emit n_synteny_blocks"
        assert result["n_synteny_blocks"] == len(result["synteny_blocks"])
        # We engineered 5 PAF rows (3 on LG28, 2 on LG07), all post-filter
        assert result["n_synteny_blocks"] == 5, \
            f"expected 5 synteny blocks, got {result['n_synteny_blocks']}"
        for blk in result["synteny_blocks"]:
            for k in ("gar_chr", "gar_start", "gar_end", "mac_chr", "mac_start",
                      "mac_end", "strand", "block_size_bp", "mapping_quality"):
                assert k in blk, f"synteny block missing {k}: {blk}"
        # Sorted by gar_chr then gar_start
        sb = result["synteny_blocks"]
        for i in range(1, len(sb)):
            if sb[i-1]["gar_chr"] == sb[i]["gar_chr"]:
                assert sb[i-1]["gar_start"] <= sb[i]["gar_start"], \
                    f"synteny_blocks not sorted by gar_start within {sb[i]['gar_chr']}"
        assert "chrom_lengths_query" in result and "chrom_lengths_target" in result
        # LG28 query length should be 10M (per synthetic), LG07 same
        assert result["chrom_lengths_query"].get("LG28") == 10_000_000, \
            f"unexpected LG28 query length: {result['chrom_lengths_query'].get('LG28')}"
        # LG10 should appear in target lengths (translocation row)
        assert "LG10" in result["chrom_lengths_target"], \
            "LG10 should appear in target chrom_lengths"
        print(f"[test] v2 schema check OK ({result['n_synteny_blocks']} synteny blocks, "
              f"{len(result['chrom_lengths_query'])} gar chroms, "
              f"{len(result['chrom_lengths_target'])} mac chroms)")

        # We engineered 2 inversion-edge breakpoints on LG28 (entering and leaving
        # the inverted block) — they may cluster into one with the default
        # cluster radius. And 1 translocation on LG07. So expect 2 or 3 total.
        bps = result["breakpoints"]
        bp_by_chrom = {}
        for bp in bps:
            bp_by_chrom.setdefault(bp["gar_chr"], []).append(bp)

        assert "LG28" in bp_by_chrom, f"LG28 breakpoints missing; got chroms {list(bp_by_chrom)}"
        assert "LG07" in bp_by_chrom, f"LG07 breakpoints missing; got chroms {list(bp_by_chrom)}"

        # LG28 should have at least one inversion-type breakpoint
        lg28_inv = [b for b in bp_by_chrom["LG28"]
                    if b["event_type"] == "inversion" or "inversion" in str(b.get("event_type"))]
        assert len(lg28_inv) >= 1, f"expected an inversion bp on LG28; got {bp_by_chrom['LG28']}"

        # LG07 should have a translocation (not fission, since one Mac chrom
        # contributes 50% — that's actually our edge case, see reclassify rule)
        lg07_trans = [b for b in bp_by_chrom["LG07"]
                      if b["event_type"] in ("translocation_or_fission", "mixed")]
        assert len(lg07_trans) >= 1, f"expected translocation bp on LG07; got {bp_by_chrom['LG07']}"
        # Check refined classification ran
        assert "event_type_refined" in lg07_trans[0], \
            f"missing event_type_refined on translocation; got {lg07_trans[0]}"

        # Spalax-note check: the LG28 breakpoint at ~3 Mb should be inside the
        # all_TE peak (0.92), so flanking_repeat_density_gar.all_TE.mean should
        # be elevated. The manuscript_note is gated by a 1.5x chrom-wide ratio;
        # with a tall narrow peak in a small synthetic chrom, the chrom-wide
        # mean is itself elevated so the ratio may not fire — we don't assert
        # the note's presence, just that the local flank is enriched vs
        # the LG07 (flat) baseline.
        for bp in lg28_inv:
            gar_rd = bp.get("flanking_repeat_density_gar", {})
            all_te = gar_rd.get("all_TE", {})
            assert all_te, f"missing all_TE flanking density on LG28 inversion bp: {bp}"
            print(f"[test] LG28 inv bp {bp['id']}: gar all_TE flank mean = {all_te['mean']}, max = {all_te['max']}")
            if 2_500_000 <= (bp["gar_pos_start"] + bp["gar_pos_end"]) // 2 <= 5_500_000:
                # Flanking mean should be visibly elevated above the LG07 baseline (0.45)
                assert all_te["mean"] >= 0.55, \
                    f"expected elevated all_TE flanking mean (>0.55); got {all_te['mean']}"
                # Max should reach near the engineered peak (0.92)
                assert all_te["max"] >= 0.85, \
                    f"expected all_TE flanking max near peak (>=0.85); got {all_te['max']}"

        # LG07 control: flanking should be flat at the engineered baseline
        for bp in bp_by_chrom["LG07"]:
            all_te = bp.get("flanking_repeat_density_gar", {}).get("all_TE", {})
            if all_te:
                assert 0.40 <= all_te["mean"] <= 0.50, \
                    f"expected flat baseline on LG07; got mean={all_te['mean']}"

        # Candidate overlap check: at least one LG28 breakpoint should overlap CAND001
        lg28_overlap = [b for b in bp_by_chrom["LG28"] if "candidate_overlap" in b and "CAND001" in b["candidate_overlap"]]
        assert len(lg28_overlap) >= 1, "expected LG28 breakpoint to overlap CAND001"

        # Schema check: required top-level fields
        for k in ("tool", "schema_version", "generated_at", "species_query", "species_target",
                  "input_paf", "params", "n_breakpoints", "n_by_event_type", "breakpoints"):
            assert k in result, f"missing top-level field: {k}"

        # Per-bp schema check
        for bp in bps:
            for k in ("id", "event_type", "gar_chr", "gar_pos_start", "gar_pos_end",
                      "gar_pos_mb", "prev_block", "next_block",
                      "flanking_repeat_density_gar", "flanking_repeat_density_mac"):
                assert k in bp, f"bp {bp.get('id')} missing field: {k}"

        print(f"[test] PASS — n_breakpoints={result['n_breakpoints']}, "
              f"event_types={result['n_by_event_type']}")
        for bp in bps:
            note = bp.get("manuscript_note", "")
            print(f"  {bp['id']:<25} {bp['gar_chr']:<6} @ {bp['gar_pos_mb']:>7.2f} Mb  "
                  f"{bp.get('event_type_refined', bp['event_type']):<28}  {note}")


if __name__ == "__main__":
    main()
