#!/usr/bin/env python3
"""
STEP_SV_EVID_COMB_emit_combinations.py
--------------------------------------
Build sv_evidence_combinations_v1.json (the UpSet-panel layer) for one
candidate. The atlas page renders one bar per evidence combination, with
height ∝ intersection_size and click → selectedSamples handshake.

Inputs:
  --evidence    path to per-sample × per-evidence-type 0/1 matrix (TSV)
                Header row: sample_id <tab> ev_id1 <tab> ev_id2 ...
                Body row:   sample_id <tab> 0|1 <tab> 0|1 ...
                Cells are 0/1; missing data should be encoded as 0.
                The producer (typically MODULE_5A2 / S7) writes this TSV.

  --evidence-types  path to evidence-types JSON: a list of dicts with
                    {id, label, side, kind, tier}. Order is preserved as
                    the canonical row order in the UpSet plot.

  --candidate-id    candidate identifier (will be embedded in output JSON
                    AND used for the per-candidate folder name).

  --chrom           pseudochromosome label (used in output folder path).

  --out-root        output root directory.

  --top-n           keep only the top-N combinations by intersection_size
                    (default 20). Smaller = cleaner UpSet panel.

  --indent          JSON indent (default: compact).

Output:
  <out>/<chrom>/candidates/<candidate_id>/sv_evidence_combinations.json

Worked example (matches the test fixture in tests/sv_evidence/):
  evidence_types = [left_SA, right_SA, left_PE, right_PE,
                    Manta_INV_GT, DELLY_INV_GT, MAPQ0_left, MAPQ0_right]
  226 fish; 4 dominant combinations:
    42  fish: left_SA + right_SA + left_PE + right_PE  (strongest support)
    21  fish: MAPQ0 only                                (likely artefact)
    17  fish: Manta_INV_GT + DELLY_INV_GT               (caller-only)
     8  fish: left breakpoint only                      (single-sided)

Author: Quentin Andres
Project: MS_Inversions_North_african_catfish
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

# Local sibling import
sys.path.insert(0, str(Path(__file__).parent))
from write_candidate_folder import write_candidate_layer, update_chrom_manifest, log_write


def _read_evidence_matrix(path: str | Path) -> tuple[list[str], dict[str, list[int]]]:
    """
    Read the per-sample × per-evidence-type 0/1 matrix from a TSV file.
    Returns (evidence_type_ids, dict(sample_id → [0/1, 0/1, ...]))
    where the inner list has one entry per evidence_type_id in the header
    column order.

    Validation:
      - First column must be sample_id
      - Every cell in body rows must be 0 or 1 (else cast and warn)
    """
    ev_ids: list[str] = []
    rows: dict[str, list[int]] = {}
    with open(path) as f:
        header = None
        for ln, line in enumerate(f, 1):
            if line.startswith("#") or not line.strip():
                continue
            cells = line.rstrip("\n").split("\t")
            if header is None:
                header = cells
                if header[0] != "sample_id":
                    print(f"[STEP_SV_EVID_COMB] expected first column 'sample_id', "
                          f"got {header[0]!r}", file=sys.stderr)
                ev_ids = header[1:]
                continue
            sid = cells[0]
            vec = []
            for c in cells[1:]:
                try:
                    v = int(c)
                except ValueError:
                    v = 0
                vec.append(1 if v else 0)
            # Pad / truncate to evidence-type count if mismatched
            if len(vec) < len(ev_ids):
                vec.extend([0] * (len(ev_ids) - len(vec)))
            elif len(vec) > len(ev_ids):
                vec = vec[: len(ev_ids)]
            rows[sid] = vec
    return ev_ids, rows


def build_combinations(ev_ids: list[str], matrix: dict[str, list[int]],
                       evidence_types: list[dict], top_n: int = 20) -> dict:
    """
    Walk samples; for each sample compute the set of evidence_ids it has
    (where the bit is 1), then group samples by their evidence-set, then
    sort by group size desc and keep the top_n.

    Per-evidence totals are computed independently (not just within
    top_n bars — every evidence type's total carrier count).

    Returns the sv_evidence_combinations_v1 payload.
    """
    # Reindex evidence_types to make sure order matches the matrix header
    # (the producer trusts the JSON's order, so we reorder if needed).
    ev_id_to_meta = {e["id"]: e for e in evidence_types}
    # If the matrix header has a different order, we'll project to the
    # canonical evidence_types order
    ev_ids_canonical = [e["id"] for e in evidence_types]
    # Map header index → canonical index (or -1 for unknown header columns)
    proj = []
    for i, hid in enumerate(ev_ids):
        if hid in ev_id_to_meta:
            proj.append(ev_ids_canonical.index(hid))
        else:
            proj.append(-1)
            print(f"[STEP_SV_EVID_COMB] header column {hid!r} not in evidence_types; "
                  f"will be ignored", file=sys.stderr)

    n_canonical = len(ev_ids_canonical)
    # Group samples by their bit-pattern (canonical order)
    by_pattern: dict[tuple[int, ...], list[str]] = defaultdict(list)
    per_evidence_carriers: dict[str, int] = defaultdict(int)

    for sid, vec in matrix.items():
        canon = [0] * n_canonical
        for i, v in enumerate(vec):
            ci = proj[i]
            if ci >= 0:
                canon[ci] = canon[ci] | v
        # Skip samples with no evidence at all (they're not informative for
        # the UpSet panel; they'd land in the "empty intersection" bucket
        # which is conventionally hidden).
        if not any(canon):
            continue
        by_pattern[tuple(canon)].append(sid)
        for i, v in enumerate(canon):
            if v:
                per_evidence_carriers[ev_ids_canonical[i]] += 1

    # Convert patterns to combinations
    combos = []
    for pattern, sids in by_pattern.items():
        members = [ev_ids_canonical[i] for i, v in enumerate(pattern) if v]
        combos.append({
            "members":           members,
            "intersection_size": len(sids),
            "samples":           sorted(sids),
        })
    # Sort by size desc, then by tier (lower = stronger evidence first
    # in stable order), then alphabetic
    combos.sort(key=lambda c: (-c["intersection_size"],
                               len(c["members"]),
                               c["members"]))
    if top_n and len(combos) > top_n:
        combos = combos[:top_n]

    payload = {
        "format_version":     "sv_evidence_combinations_v1",
        "n_samples_total":    len(matrix),
        "evidence_types":     evidence_types,
        "combinations":       combos,
        "per_evidence_totals": {
            eid: {"n_samples": per_evidence_carriers.get(eid, 0)}
            for eid in ev_ids_canonical
        },
    }
    return payload


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--evidence", required=True,
                    help="per-sample × per-evidence-type 0/1 TSV")
    ap.add_argument("--evidence-types", required=True,
                    help="evidence-types JSON (list of {id,label,side,kind,tier})")
    ap.add_argument("--candidate-id", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--out-root", required=True)
    ap.add_argument("--top-n", type=int, default=20)
    ap.add_argument("--indent", type=int, default=None)
    ap.add_argument("--no-manifest", action="store_true")
    args = ap.parse_args()

    evidence_types = json.loads(Path(args.evidence_types).read_text())
    if not isinstance(evidence_types, list):
        ap.error("--evidence-types must be a JSON list")

    ev_ids_header, matrix = _read_evidence_matrix(args.evidence)
    print(f"[STEP_SV_EVID_COMB] {len(matrix)} samples × {len(ev_ids_header)} evidence types",
          file=sys.stderr)

    payload = build_combinations(ev_ids_header, matrix, evidence_types, top_n=args.top_n)
    payload["candidate_id"] = args.candidate_id
    print(f"[STEP_SV_EVID_COMB] {len(payload['combinations'])} combinations after top-N",
          file=sys.stderr)

    out = write_candidate_layer(payload, args.out_root, args.chrom,
                                args.candidate_id, indent=args.indent)
    log_write(out)

    if not args.no_manifest:
        m = update_chrom_manifest(args.out_root, args.chrom)
        log_write(m)


if __name__ == "__main__":
    main()
