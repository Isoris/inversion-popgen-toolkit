"""
producers/STEP_SV_SUPPORT_emit_support_by_sample.py
====================================================
Emit `sv_support_by_sample_v1.json` for one candidate.

This is the per-sample × per-SV dosage matrix that powers the atlas's
heatmap (step 6) and unlocks "dim non-selected SV glyphs when an UpSet
bar is active" in step 5.

Reads (per candidate):
  - The same merged DELLY+Manta VCF that STEP_SV_GT_AGG reads
  - The candidate's surviving SV list (post-min-carriers filter) —
    typically the SVs already emitted in sv_genotype_counts.json
  - locked_labels for sample order + group split

Writes:
  - <out_root>/<chrom>/candidates/<cand_id>/sv_support_by_sample.json

Format (compact-string row encoding):
    {
      "format_version": "sv_support_by_sample_v1",
      "candidate_id":   "...",
      "encoding":       "0=AA, 1=AB, 2=BB, .=miss",
      "samples":        ["s_..", "s_..", ...],   // ordered: H1/H1 → H1/H2 → H2/H2
      "sv_ids":         ["SV001", ...],
      "row_groups":     {"H1/H1": [0, 60], "H1/H2": [61, 163], "H2/H2": [164, 225]},
      "dosage_compact": ["00..0...", ...]         // one string per sample
    }

Char encoding chosen so:
  - producer side is trivial: ''.join(map(_d_to_ch, doses))
  - JSON stays small (1 byte/cell vs 3 for [0,0,1,...])
  - parser is O(n) charAt in JS

THIS IS A STUB — wire the per-sample dosage extraction against your
actual VCF library. The JSON shape is locked.

Usage:
    python STEP_SV_SUPPORT_emit_support_by_sample.py \\
        --vcf            merged_delly_manta.vcf.gz \\
        --sv-list        cand_LG28_15Mb__sv_ids.txt \\
        --locks          locks/cand_LG28_15Mb.json \\
        --candidate-id   cand_LG28_15Mb \\
        --chrom          C_gar_LG28 \\
        --out-root       data/
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


# Dosage → char
_DOSAGE_CH = {0: '0', 1: '1', 2: '2', -1: '.', None: '.'}


def encode_row(dosages: list[int]) -> str:
    """Compact-string encoding for one sample's row across all SVs."""
    return ''.join(_DOSAGE_CH.get(d, '.') for d in dosages)


def build_support_layer(
    sample_dosages: dict[str, list[int]],
    sample_order:   list[str],
    sv_ids:         list[str],
    group_ranges:   dict[str, tuple[int, int]],
    candidate_id:   str,
) -> dict:
    """
    `sample_dosages`: {sample_id: [dosage for each sv_id, in `sv_ids` order]}
    `sample_order`:   list of sample_ids ordered as H1/H1 → H1/H2 → H2/H2.
                      This MUST match the order of locked_labels grouped by label.
    `sv_ids`:         column order for each row's dosage vector.
    `group_ranges`:   {'H1/H1': (start_idx, end_idx), 'H1/H2': ..., 'H2/H2': ...}
                      inclusive end. Indices into `sample_order`.
    """
    rows = []
    for sid in sample_order:
        doses = sample_dosages.get(sid)
        if doses is None:
            doses = [-1] * len(sv_ids)
        if len(doses) != len(sv_ids):
            raise ValueError(
                f"row length mismatch for sample {sid!r}: "
                f"got {len(doses)} dosages, expected {len(sv_ids)}"
            )
        rows.append(encode_row(doses))

    return {
        'format_version': 'sv_support_by_sample_v1',
        'candidate_id':   candidate_id,
        'encoding':       '0=AA, 1=AB, 2=BB, .=miss',
        'samples':        sample_order,
        'sv_ids':         sv_ids,
        'row_groups': {
            g: list(rng) for g, rng in group_ranges.items()
        },
        'dosage_compact': rows,
    }


def order_samples_by_group(locks: list[dict]) -> tuple[list[str], dict[str, tuple[int, int]]]:
    """
    Re-order locked_labels into the canonical H1/H1 → H1/H2 → H2/H2 order
    expected by the atlas. Returns (sample_order, group_ranges).
    """
    label_to_group = {'HOMO_1': 'H1/H1', 'HET': 'H1/H2', 'HOMO_2': 'H2/H2'}
    grouped = {'H1/H1': [], 'H1/H2': [], 'H2/H2': []}
    for entry in locks:
        g = label_to_group.get(entry['label'])
        if g in grouped:
            grouped[g].append(entry['sample_id'])
    sample_order = grouped['H1/H1'] + grouped['H1/H2'] + grouped['H2/H2']
    ranges = {}
    cumul = 0
    for g in ('H1/H1', 'H1/H2', 'H2/H2'):
        n = len(grouped[g])
        if n > 0:
            ranges[g] = (cumul, cumul + n - 1)
        else:
            # Empty group — still record an empty range so the atlas can
            # render group separators consistently.
            ranges[g] = (cumul, cumul - 1)
        cumul += n
    return sample_order, ranges


# ---------------------------------------------------------------------------
# VCF dosage extraction — STUB
# ---------------------------------------------------------------------------

def extract_dosages_from_vcf(vcf_path: str, sv_ids: list[str],
                              sample_order: list[str]) -> dict[str, list[int]]:
    """
    Replace with cyvcf2/pysam:

        from cyvcf2 import VCF
        v = VCF(vcf_path)
        sample_idx = {s: i for i, s in enumerate(v.samples)}
        out = {sid: [] for sid in sample_order}
        wanted = set(sv_ids)
        # Iterate VCF; for each record matching wanted: read GT codes
        # and write them into out[sample] in the order of sv_ids.
        ...

    Returns {sample_id: [dosage_for_sv_in_order_of_sv_ids]}.
    """
    raise SystemExit("STEP_SV_SUPPORT: VCF parsing not implemented. "
                     "Wire to cyvcf2/pysam; see top-of-file docstring.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--vcf',          required=True)
    ap.add_argument('--sv-list',      required=True,
                    help='one sv_id per line; the SVs that survived '
                         'min_carriers filtering in STEP_SV_GT_AGG')
    ap.add_argument('--locks',        required=True)
    ap.add_argument('--candidate-id', required=True)
    ap.add_argument('--chrom',        required=True)
    ap.add_argument('--out-root',     default='data/')
    args = ap.parse_args()

    with open(args.locks) as f:
        locks = json.load(f)
    if isinstance(locks, dict) and 'locked_labels' in locks:
        locks = locks['locked_labels']
    sample_order, group_ranges = order_samples_by_group(locks)

    with open(args.sv_list) as f:
        sv_ids = [ln.strip() for ln in f if ln.strip()]

    sample_dosages = extract_dosages_from_vcf(args.vcf, sv_ids, sample_order)

    payload = build_support_layer(
        sample_dosages, sample_order, sv_ids, group_ranges,
        candidate_id=args.candidate_id,
    )

    from sv_evidence_io import CandidateDir, write_layer
    cand_dir = CandidateDir(Path(args.out_root), args.chrom, args.candidate_id)
    written = write_layer(cand_dir, 'sv_support_by_sample', payload)
    print(f"wrote {written} ({len(payload['samples'])} samples × {len(payload['sv_ids'])} SVs)")


if __name__ == '__main__':
    main()
