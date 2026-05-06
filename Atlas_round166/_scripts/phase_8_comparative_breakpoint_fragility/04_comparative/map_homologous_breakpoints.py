#!/usr/bin/env python3
# =============================================================================
# 04_comparative/map_homologous_breakpoints.py
# =============================================================================
"""
Reads `config/synteny_mapping.tsv`, validates it against the species and
candidate manifests, and emits a per-candidate JSON describing the homologous
region(s) in each target species:

    output/comparative/map/<candidate_id>.homologs.json

This script does not compute density. It only resolves and validates the
mapping, so step 03 (compute_te_density.py) — which already runs against
arbitrary breakpoint windows — can be reused on the comparative window file
produced by step 02 with --synteny.

Validation:
  - candidate_id must exist in candidate_breakpoints.tsv
  - target_species must exist in species_manifest.tsv
  - target_chrom must exist in that species' .fai
  - target_start < target_end after swapping if needed
  - orientation in {+, -}
  - boundary_type in {internal, terminal, block_edge, chr_context_change} or empty

If a candidate has no synteny rows, the output file still gets written with an
empty `targets` list — useful so downstream code never has to special-case the
"no homologs at all" case.

Usage:
    python map_homologous_breakpoints.py \
        --candidates config/candidate_breakpoints.tsv \
        --species_manifest config/species_manifest.tsv \
        --synteny config/synteny_mapping.tsv \
        --out_dir output/comparative/map
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

_THIS = Path(__file__).resolve()
sys.path.insert(0, str(_THIS.parent.parent))
from common import (  # noqa: E402
    get_logger,
    normalize_chrom,
    read_tsv,
    write_run_manifest,
)


VALID_BOUNDARY = {"", "internal", "terminal", "block_edge", "chr_context_change"}


def load_chr_lens(fai_path: Path) -> dict[str, int]:
    out: dict[str, int] = {}
    with fai_path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            try:
                out[normalize_chrom(parts[0])] = int(parts[1])
            except ValueError:
                continue
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidates", required=True, type=Path)
    ap.add_argument("--species_manifest", required=True, type=Path)
    ap.add_argument("--synteny", required=True, type=Path)
    ap.add_argument("--out_dir", required=True, type=Path)
    ap.add_argument("--log_dir", type=Path, default=Path("output/logs"))
    args = ap.parse_args()

    log = get_logger("04_map_homologous_breakpoints", args.log_dir)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    cands = read_tsv(args.candidates,
                     required_cols=["candidate_id", "chrom",
                                    "left_breakpoint", "right_breakpoint"])
    cand_ids = {r["candidate_id"] for r in cands}
    cand_meta = {r["candidate_id"]: r for r in cands}

    species_rows = read_tsv(args.species_manifest,
                            required_cols=["species_id", "fai_path", "is_focal"])
    species_idx = {r["species_id"]: r for r in species_rows}
    focal_id = next(r["species_id"] for r in species_rows
                    if str(r["is_focal"]).lower() == "true")
    log.info("focal=%s; %d species in manifest", focal_id, len(species_idx))

    # cache fai per species lazily
    fai_cache: dict[str, dict[str, int]] = {}

    def chr_lens_for(sp: str) -> dict[str, int]:
        if sp not in fai_cache:
            fai_cache[sp] = load_chr_lens(Path(species_idx[sp]["fai_path"]))
        return fai_cache[sp]

    synt = read_tsv(args.synteny, required_cols=[
        "candidate_id", "target_species", "target_chrom",
        "target_start", "target_end", "orientation",
    ])
    log.info("loaded %d synteny rows", len(synt))

    by_cand: dict[str, list[dict]] = defaultdict(list)
    n_invalid = 0
    for r in synt:
        cid = r["candidate_id"]
        tsp = r["target_species"]
        tch = normalize_chrom(r["target_chrom"])

        if cid not in cand_ids:
            log.warning("synteny row references unknown candidate_id=%s; dropped", cid)
            n_invalid += 1
            continue
        if tsp not in species_idx:
            log.warning("synteny row references unknown target_species=%s; dropped", tsp)
            n_invalid += 1
            continue
        if tsp == focal_id:
            log.warning("synteny row maps candidate %s to focal species %s; dropped",
                        cid, tsp)
            n_invalid += 1
            continue

        try:
            ts = int(r["target_start"])
            te = int(r["target_end"])
        except ValueError:
            log.warning("non-integer target_start/end for %s -> %s; dropped",
                        cid, tsp)
            n_invalid += 1
            continue
        if ts > te:
            ts, te = te, ts

        chr_lens = chr_lens_for(tsp)
        if tch not in chr_lens:
            log.warning("target_chrom %s not in %s .fai; dropped",
                        tch, tsp)
            n_invalid += 1
            continue
        if te > chr_lens[tch]:
            log.warning("target_end %d > chrom_len %d for %s/%s; clipping",
                        te, chr_lens[tch], tsp, tch)
            te = chr_lens[tch]

        orient = (r.get("orientation") or "").strip()
        if orient not in {"+", "-"}:
            log.warning("invalid orientation=%r for %s -> %s; coercing to '+'",
                        orient, cid, tsp)
            orient = "+"

        btype = (r.get("boundary_type") or "").strip()
        if btype not in VALID_BOUNDARY:
            log.warning("unknown boundary_type=%r; coercing to 'internal'", btype)
            btype = "internal"

        by_cand[cid].append({
            "species_id": tsp,
            "assembly_name": species_idx[tsp].get("assembly_name") or "",
            "homologous_region": {
                "chrom": tch,
                "start": ts,
                "end": te,
                "orientation_relative_to_focal": orient,
                "synteny_block_id": (r.get("synteny_block_id") or "") or None,
                "boundary_type": btype or "internal",
                "mapping_method": (r.get("mapping_method") or "") or None,
            },
        })

    # write one file per candidate (always, even empty)
    n_written = 0
    for cid, fr in cand_meta.items():
        targets = by_cand.get(cid, [])
        out = {
            "schema_version": "comparative_breakpoint_fragility_v0.1",
            "candidate_id": cid,
            "focal_species": focal_id,
            "focal_chrom": normalize_chrom(fr["chrom"]),
            "focal_left_breakpoint":  int(fr["left_breakpoint"]),
            "focal_right_breakpoint": int(fr["right_breakpoint"]),
            "targets": targets,
        }
        out_path = args.out_dir / f"{cid}.homologs.json"
        out_path.write_text(json.dumps(out, indent=2, sort_keys=False))
        n_written += 1

    write_run_manifest(args.out_dir / "_map.run.json",
                       n_synteny_rows=len(synt),
                       n_dropped=n_invalid,
                       n_candidates=n_written,
                       n_with_targets=sum(1 for v in by_cand.values() if v))
    log.info("wrote %d homolog files (%d with targets, %d dropped rows)",
             n_written, sum(1 for v in by_cand.values() if v), n_invalid)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
