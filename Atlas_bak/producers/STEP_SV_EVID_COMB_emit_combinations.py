"""
producers/STEP_SV_EVID_COMB_emit_combinations.py
================================================
Emit `sv_evidence_combinations_v1.json` for one candidate.

This is the layer behind the atlas's UpSet panel (step 5 of the SV
evidence page). Rows = per-SV evidence types; top bars = fish counts
per evidence combination; click a bar → those samples become the
`selectedSamples` set in the atlas.

Reads (per candidate):
  - DELLY2 + Manta GT-call status per sample
    → Manta_INV_GT, DELLY_INV_GT booleans per sample
  - BAM-evidence track JSON from MODULE_4 / S7 / phase_8_comparative_breakpoint_fragility:
      - split-read clip-counts at left/right boundary → left_SA, right_SA
      - paired-end discordant orientation reads at left/right → left_PE, right_PE
      - low-MAPQ regions overlapping each boundary → MAPQ0_left, MAPQ0_right
  - locked_labels JSON (for sample ordering and consistency check)
  - candidate metadata

Writes:
  - <out_root>/<chrom>/candidates/<cand_id>/sv_evidence_combinations.json

Pipeline:
  1. For each sample, build a vector of (evidence_type → 0/1).
  2. Group samples by their identical evidence-vector → "combinations".
  3. Sort combinations by intersection_size desc.
  4. Emit top N (default 30; spec'd cap is 30).
  5. Per-evidence-type totals for the right-side mini-bars.

THIS IS A STUB — the actual evidence-extractor functions need to be
wired against your S7 BAM-evidence pipeline. The JSON shape is locked.

Usage:
    python STEP_SV_EVID_COMB_emit_combinations.py \\
        --bam-evidence  s7_phase8/cand_LG28_15Mb.json \\
        --vcf-status    merged_delly_manta_status.json \\
        --locks         locks/cand_LG28_15Mb.json \\
        --candidate-id  cand_LG28_15Mb \\
        --chrom         C_gar_LG28 \\
        --out-root      data/
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path
from typing import Any


# Canonical evidence type definitions — the atlas displays them in this
# order (tier-by-tier: split-read first, PE, callers, MAPQ0). Producer
# is free to add more types but the atlas only knows these.
CANONICAL_EVIDENCE_TYPES = [
    {'id': 'left_SA',       'label': 'Left split-read',  'side': 'left',  'kind': 'SA',     'tier': 1},
    {'id': 'right_SA',      'label': 'Right split-read', 'side': 'right', 'kind': 'SA',     'tier': 1},
    {'id': 'left_PE',       'label': 'Left PE',          'side': 'left',  'kind': 'PE',     'tier': 2},
    {'id': 'right_PE',      'label': 'Right PE',         'side': 'right', 'kind': 'PE',     'tier': 2},
    {'id': 'Manta_INV_GT',  'label': 'Manta INV',        'side': None,    'kind': 'caller', 'tier': 3},
    {'id': 'DELLY_INV_GT',  'label': 'DELLY INV',        'side': None,    'kind': 'caller', 'tier': 3},
    {'id': 'MAPQ0_left',    'label': 'Left MAPQ0',       'side': 'left',  'kind': 'mapq0',  'tier': 4},
    {'id': 'MAPQ0_right',   'label': 'Right MAPQ0',      'side': 'right', 'kind': 'mapq0',  'tier': 4},
]


def build_evidence_combinations(
    per_sample_evidence: dict[str, dict[str, bool]],
    candidate_id: str,
    n_samples_total: int | None = None,
    top_n: int = 30,
) -> dict:
    """
    Roll up per-sample evidence → combinations.

    `per_sample_evidence`: {sample_id: {evidence_type_id: bool}}
        Every sample appears, even if all-False (UpSet's "no evidence"
        bar is meaningful — those are the samples that don't carry
        the candidate by any layer).

    Returns the sv_evidence_combinations_v1 payload.
    """
    if n_samples_total is None:
        n_samples_total = len(per_sample_evidence)

    # 1. Build canonical vector per sample (fixed order = CANONICAL_EVIDENCE_TYPES order)
    type_ids = [t['id'] for t in CANONICAL_EVIDENCE_TYPES]

    # 2. Group by vector
    combos: dict[tuple[bool, ...], list[str]] = defaultdict(list)
    for sid, ev in per_sample_evidence.items():
        vec = tuple(bool(ev.get(t, False)) for t in type_ids)
        combos[vec].append(sid)

    # 3. Convert to atlas-shape combinations, drop the all-False group
    #    (it adds noise to UpSet — display "no evidence" via a separate
    #    counter if needed). Sort by size desc.
    out_combos = []
    for vec, samples in combos.items():
        members = [tid for tid, b in zip(type_ids, vec) if b]
        if not members:
            continue
        out_combos.append({
            'members':           members,
            'intersection_size': len(samples),
            'samples':           sorted(samples),
        })
    out_combos.sort(key=lambda c: -c['intersection_size'])
    out_combos = out_combos[:top_n]

    # 4. Per-evidence-type totals
    totals = {tid: 0 for tid in type_ids}
    for sid, ev in per_sample_evidence.items():
        for tid in type_ids:
            if ev.get(tid):
                totals[tid] += 1

    return {
        'format_version':  'sv_evidence_combinations_v1',
        'candidate_id':    candidate_id,
        'n_samples_total': n_samples_total,
        'evidence_types':  CANONICAL_EVIDENCE_TYPES,
        'combinations':    out_combos,
        'per_evidence_totals': {tid: {'n_samples': totals[tid]} for tid in type_ids},
    }


# ---------------------------------------------------------------------------
# Evidence extractors — STUBS to be wired against your S7 / MODULE_4 outputs
# ---------------------------------------------------------------------------

def extract_evidence_from_bam_layer(bam_evidence: dict, sample_id: str,
                                     boundary_left_bp: int, boundary_right_bp: int,
                                     window_bp: int = 500) -> dict[str, bool]:
    """
    Look up split-read / paired-end / MAPQ0 evidence for one sample at the
    candidate's boundaries.

    Replace this stub with calls into your phase_8_comparative_breakpoint_fragility
    layer's actual API. Expected per-sample fields:
        bam_evidence[sample_id]['left_SA_count']  → int
        bam_evidence[sample_id]['right_SA_count'] → int
        bam_evidence[sample_id]['left_PE_count']  → int
        bam_evidence[sample_id]['right_PE_count'] → int
        bam_evidence[sample_id]['mapq0_left']     → bool
        bam_evidence[sample_id]['mapq0_right']    → bool

    Threshold defaults (overridable via config in your pipeline):
        SA threshold: ≥ 3 split reads at boundary ± `window_bp`
        PE threshold: ≥ 3 discordant PE reads at boundary ± `window_bp`
    """
    s = bam_evidence.get(sample_id, {})
    return {
        'left_SA':      s.get('left_SA_count',  0) >= 3,
        'right_SA':     s.get('right_SA_count', 0) >= 3,
        'left_PE':      s.get('left_PE_count',  0) >= 3,
        'right_PE':     s.get('right_PE_count', 0) >= 3,
        'MAPQ0_left':   bool(s.get('mapq0_left',  False)),
        'MAPQ0_right':  bool(s.get('mapq0_right', False)),
    }


def extract_caller_evidence(vcf_status: dict, sample_id: str) -> dict[str, bool]:
    """
    From a per-sample status JSON (built upstream from DELLY+Manta VCFs),
    determine whether each caller's INV genotype passes for this sample.

    `vcf_status[sample_id]['delly_inv_gt']` → 'AB' | 'BB' | 'AA' | None
    `vcf_status[sample_id]['manta_inv_gt']` → same
    """
    s = vcf_status.get(sample_id, {})
    delly = s.get('delly_inv_gt')
    manta = s.get('manta_inv_gt')
    return {
        'Manta_INV_GT': manta in ('AB', 'BB'),
        'DELLY_INV_GT': delly in ('AB', 'BB'),
    }


def merge_evidence(*partials: dict[str, bool]) -> dict[str, bool]:
    """Union-merge multiple per-sample evidence dicts, e.g. BAM ⊎ caller."""
    out: dict[str, bool] = {}
    for p in partials:
        out.update(p)
    return out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bam-evidence', required=True)
    ap.add_argument('--vcf-status',   required=True)
    ap.add_argument('--locks',        required=True)
    ap.add_argument('--candidate-id', required=True)
    ap.add_argument('--chrom',        required=True)
    ap.add_argument('--left',  type=int, required=True)
    ap.add_argument('--right', type=int, required=True)
    ap.add_argument('--out-root', default='data/')
    ap.add_argument('--top-n',    type=int, default=30)
    args = ap.parse_args()

    # Load
    with open(args.bam_evidence) as f:
        bam = json.load(f)
    with open(args.vcf_status) as f:
        vcf = json.load(f)
    with open(args.locks) as f:
        locks = json.load(f)
    if isinstance(locks, dict) and 'locked_labels' in locks:
        locks = locks['locked_labels']

    # Build per-sample evidence vector
    per_sample = {}
    for entry in locks:
        sid = entry['sample_id']
        per_sample[sid] = merge_evidence(
            extract_evidence_from_bam_layer(bam, sid, args.left, args.right),
            extract_caller_evidence(vcf, sid),
        )

    # Build payload
    payload = build_evidence_combinations(
        per_sample, args.candidate_id, n_samples_total=len(locks),
        top_n=args.top_n,
    )

    # Write to per-candidate folder
    from sv_evidence_io import CandidateDir, write_layer
    cand_dir = CandidateDir(Path(args.out_root), args.chrom, args.candidate_id)
    written = write_layer(cand_dir, 'sv_evidence_combinations', payload)
    print(f"wrote {written} ({len(payload['combinations'])} combinations)")


if __name__ == '__main__':
    main()
