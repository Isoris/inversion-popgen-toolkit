#!/usr/bin/env python3
"""
STEP_SV_SUPPORT_emit_support_by_sample.py
-----------------------------------------
Build sv_support_by_sample_v1.json (the third per-candidate layer) for
one candidate. This is the per-sample × per-SV dosage matrix the atlas's
heatmap (step 6) and the deferred-feature unlocks read from:

  - heatmap rendering (compact + large overlay)
  - selection-scoped genotype-count recount in the SV table
  - dim non-carriers in the locus-glyph track when an UpSet bar is active
  - per-cell hover tooltip on the heatmap

Inputs:
  --tsv         Same calls TSV format STEP_SV_GT_AGG accepts:
                  sv_id, chrom, position_bp, sv_type, sample_id, GT
                Optional: end_bp, quality, callers
  --candidate   Path to candidate JSON (must include candidate_id, chrom).
                The candidate window (boundary_left/right ± flank) is used
                only as a positional filter to keep the matrix small.
  --karyotype   Karyotype labels TSV (sample_id, label). Same format as the
                gt aggregator: HOMO_1 / HET / HOMO_2.
  --out-root    Output root. JSON written to:
                  <out>/<chrom>/candidates/<candidate_id>/sv_support_by_sample.json
  --sv-list     Optional: path to a one-sv_id-per-line file. When present,
                only those SVs are included (in file order). When absent,
                all SVs in the calls TSV that fall in the candidate window
                are included, sorted by position_bp.
  --indent      JSON indent (default: compact one-line).

Output schema (compact-string row form, decided in STEP6_DESIGN_NOTE.md):
  {
    "format_version": "sv_support_by_sample_v1",
    "candidate_id":   "...",
    "encoding":       "0=AA, 1=AB, 2=BB, .=miss",
    "samples":        [...],            # ordered: H1/H1 → H1/H2 → H2/H2
    "sv_ids":         [...],
    "row_groups":     {"H1/H1": [start_idx, end_idx_inclusive], ...},
    "dosage_compact": ["00...0...", ...] # one string per sample; len == n_svs
  }

Cohort note: this producer is for the **226-sample pure C. gariepinus
hatchery cohort** only. K clusters reflect hatchery broodlines, NOT
species admixture. Do NOT use against the F1 hybrid cohort or the
C. macrocephalus wild cohort.

Author:  Quentin Andres
Project: MS_Inversions_North_african_catfish (atlas SV evidence page)
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# Local sibling import — same layout as the other producers in this folder.
sys.path.insert(0, str(Path(__file__).parent))
from write_candidate_folder import write_candidate_layer, update_chrom_manifest, log_write


# Karyotype label remap — must match STEP_SV_GT_AGG and the atlas.
KARYOTYPE_REMAP = {
    "HOMO_1": "H1/H1",
    "HET":    "H1/H2",
    "HOMO_2": "H2/H2",
    # Accept already-remapped labels too (idempotent).
    "H1/H1":  "H1/H1",
    "H1/H2":  "H1/H2",
    "H2/H2":  "H2/H2",
}

# Genotype char encoding: 0=AA, 1=AB, 2=BB, '.'=miss.
_GT_TO_CH = {
    "AA": "0",
    "AB": "1",
    "BB": "2",
    "miss": ".",
}


def _gt_to_code(gt: str | None) -> str:
    """
    Map a VCF GT string to one of: 'AA' (0/0), 'AB' (0/1 or 1/0),
    'BB' (1/1), 'miss' (./. or ambiguous).

    Same logic as STEP_SV_GT_AGG so the dosage encoding agrees with
    the per-group genotype counts.
    """
    if not gt or gt in (".", "./.", ".|."):
        return "miss"
    g = gt.replace("|", "/")
    parts = g.split("/")
    if len(parts) != 2 or any(p == "." for p in parts):
        return "miss"
    try:
        a, b = int(parts[0]), int(parts[1])
    except ValueError:
        return "miss"
    if a == 0 and b == 0:
        return "AA"
    if a == 1 and b == 1:
        return "BB"
    return "AB"


def _read_karyotype(path: str | Path) -> dict[str, str]:
    """
    Read karyotype TSV (sample_id, label). Returns dict of
    sample_id → group ('H1/H1', 'H1/H2', or 'H2/H2'). Skips header.
    Skips samples with unrecognized labels (logs a warning to stderr).
    """
    out = {}
    seen_header = False
    with open(path) as f:
        for ln, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            # Detect a "sample_id\tlabel" header (case-insensitive).
            if not seen_header and parts[0].lower() in ("sample_id", "sample", "id"):
                seen_header = True
                continue
            seen_header = True
            sid, label = parts[0], parts[1]
            grp = KARYOTYPE_REMAP.get(label)
            if grp is None:
                print(f"[STEP_SV_SUPPORT] line {ln}: skipping unknown label "
                      f"{label!r} for {sid}", file=sys.stderr)
                continue
            out[sid] = grp
    return out


def _read_calls_tsv(path: str | Path):
    """
    Read a tabular TSV of SV calls. Required columns:
      sv_id, chrom, position_bp, sv_type, sample_id, GT

    Yields one record per SV (grouped by sv_id):
      {sv_id, chrom, position_bp, end_bp, sv_type,
       per_sample: dict(sample_id → 'AA'|'AB'|'BB'|'miss')}

    Same parser semantics as STEP_SV_GT_AGG._read_calls_tsv. Kept inline
    here (not factored out into a shared module) to keep producers
    self-contained on the LANTA HPC clone.
    """
    by_sv: dict[str, dict] = {}
    header = None
    with open(path) as f:
        for ln, line in enumerate(f, 1):
            if line.startswith("#"):
                continue
            line = line.rstrip("\n")
            if header is None:
                header = line.split("\t")
                continue
            cells = line.split("\t")
            row = dict(zip(header, cells))
            sv_id = row["sv_id"]
            rec = by_sv.get(sv_id)
            if rec is None:
                rec = {
                    "sv_id":       sv_id,
                    "chrom":       row["chrom"],
                    "position_bp": int(row["position_bp"]),
                    "end_bp":      int(row["end_bp"]) if row.get("end_bp") not in (None, "", ".") else None,
                    "sv_type":     row["sv_type"],
                    "per_sample":  {},
                }
                by_sv[sv_id] = rec
            rec["per_sample"][row["sample_id"]] = _gt_to_code(row.get("GT"))
    return list(by_sv.values())


def _candidate_window(cand: dict, flank_bp: int = 500_000) -> tuple[int, int]:
    """
    Compute the [start, end] bp window for filtering SV calls. Prefers
    the candidate's zone_definitions_bp.left_flank[0] / right_flank[1]
    if present; otherwise falls back to boundary ± flank_bp.
    """
    z = cand.get("zone_definitions_bp")
    if z and z.get("left_flank") and z.get("right_flank"):
        return int(z["left_flank"][0]), int(z["right_flank"][1])
    bL = int(cand["boundary_left_bp"])
    bR = int(cand["boundary_right_bp"])
    return bL - flank_bp, bR + flank_bp


def _order_samples_by_group(karyotype: dict[str, str]) -> tuple[list[str], dict[str, list[int]]]:
    """
    Re-order karyotype labels into the canonical H1/H1 → H1/H2 → H2/H2
    sample order expected by the atlas. Returns:
      (sample_order, row_groups)
    where row_groups is {group: [start_idx, end_idx_inclusive]}.

    Empty groups still produce a [start, start-1] empty range so the
    consumer can render group separators consistently.
    """
    grouped = {"H1/H1": [], "H1/H2": [], "H2/H2": []}
    # Stable lexicographic sort within group for reproducibility.
    for sid, g in karyotype.items():
        if g in grouped:
            grouped[g].append(sid)
    for g in grouped:
        grouped[g].sort()
    sample_order = grouped["H1/H1"] + grouped["H1/H2"] + grouped["H2/H2"]
    row_groups: dict[str, list[int]] = {}
    cumul = 0
    for g in ("H1/H1", "H1/H2", "H2/H2"):
        n = len(grouped[g])
        if n > 0:
            row_groups[g] = [cumul, cumul + n - 1]
        else:
            row_groups[g] = [cumul, cumul - 1]
        cumul += n
    return sample_order, row_groups


def _encode_row(per_sample_calls: list[dict], sv_ids: list[str], sample_id: str) -> str:
    """
    Build the compact-string row for one sample across the canonical
    sv_ids list. Each char: '0'=AA, '1'=AB, '2'=BB, '.'=miss. If the
    sample has no entry for an SV, encode as missing.
    """
    chars: list[str] = []
    for sv in per_sample_calls:
        code = sv["per_sample"].get(sample_id, "miss")
        chars.append(_GT_TO_CH.get(code, "."))
    return "".join(chars)


def build_support_layer(calls: list[dict], karyotype: dict[str, str],
                        cand: dict,
                        sv_filter: list[str] | None = None,
                        flank_bp: int = 500_000) -> dict:
    """
    Build the sv_support_by_sample_v1 payload from raw SV calls + karyotype
    + candidate JSON.

    `calls` is the output of _read_calls_tsv (list of per-SV records).
    `karyotype` is the output of _read_karyotype.
    `sv_filter`, if provided, restricts the matrix to the named SVs in
    the order given (else: all SVs in the candidate window, position-sorted).
    """
    win_start, win_end = _candidate_window(cand, flank_bp=flank_bp)
    # Filter calls to the candidate window and (optionally) the SV whitelist.
    in_window = [c for c in calls if win_start <= c["position_bp"] <= win_end]
    if sv_filter is not None:
        wanted = set(sv_filter)
        in_window = [c for c in in_window if c["sv_id"] in wanted]
        # Preserve the user's sv_filter order
        order = {sid: i for i, sid in enumerate(sv_filter)}
        in_window.sort(key=lambda c: order.get(c["sv_id"], 1 << 30))
    else:
        in_window.sort(key=lambda c: c["position_bp"])

    sv_ids = [c["sv_id"] for c in in_window]
    sample_order, row_groups = _order_samples_by_group(karyotype)

    dosage_compact = []
    for sid in sample_order:
        dosage_compact.append(_encode_row(in_window, sv_ids, sid))

    return {
        "format_version": "sv_support_by_sample_v1",
        "candidate_id":   cand["candidate_id"],
        "encoding":       "0=AA, 1=AB, 2=BB, .=miss",
        "samples":        sample_order,
        "sv_ids":         sv_ids,
        "row_groups":     row_groups,
        "dosage_compact": dosage_compact,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tsv",        required=True,
                    help="SV calls TSV (sv_id, chrom, position_bp, sv_type, sample_id, GT)")
    ap.add_argument("--candidate",  required=True, help="candidate JSON")
    ap.add_argument("--karyotype",  required=True, help="karyotype TSV (sample_id, label)")
    ap.add_argument("--out-root",   required=True, help="output root directory")
    ap.add_argument("--sv-list",    default=None,
                    help="optional one-sv_id-per-line whitelist (preserves order)")
    ap.add_argument("--flank-bp",   type=int, default=500_000,
                    help="bp flank around boundaries when zone_definitions_bp absent")
    ap.add_argument("--indent",     type=int, default=None,
                    help="JSON indent (default: compact one-line)")
    args = ap.parse_args()

    cand = json.loads(Path(args.candidate).read_text())
    if "candidate_id" not in cand or "chrom" not in cand:
        raise SystemExit("candidate JSON must include candidate_id and chrom")

    karyotype = _read_karyotype(args.karyotype)
    if not karyotype:
        raise SystemExit(f"no usable karyotype labels in {args.karyotype}")

    calls = _read_calls_tsv(args.tsv)

    sv_filter = None
    if args.sv_list:
        sv_filter = [ln.strip() for ln in Path(args.sv_list).read_text().splitlines()
                     if ln.strip() and not ln.startswith("#")]

    payload = build_support_layer(calls, karyotype, cand,
                                  sv_filter=sv_filter, flank_bp=args.flank_bp)

    written = write_candidate_layer(
        payload,
        out_root=args.out_root,
        chrom=cand["chrom"],
        candidate_id=cand["candidate_id"],
        indent=args.indent,
    )
    log_write(written)

    # Update chrom-level manifest so the atlas knows this candidate now
    # has the support layer too.
    update_chrom_manifest(args.out_root, cand["chrom"])

    n_svs = len(payload["sv_ids"])
    n_samples = len(payload["samples"])
    print(f"[STEP_SV_SUPPORT] {payload['candidate_id']}: "
          f"{n_samples} samples × {n_svs} SVs → {written}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
