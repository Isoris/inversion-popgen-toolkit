"""
producers/STEP_SV_GT_AGG_aggregate_genotype_counts.py
=====================================================
Emit `sv_genotype_counts_v1.json` for one candidate.

Reads:
  - DELLY2 + Manta VCFs (or merged VCF) for the candidate's region
  - locked_labels JSON (karyotype K=3 assignment per sample)
  - candidate metadata: candidate_id, chrom, boundary_left_bp,
                        boundary_right_bp, zone definitions

Writes:
  - <out_root>/<chrom>/candidates/<cand_id>/sv_genotype_counts.json

Per-SV pipeline:
    1. Parse VCF → per-sample 0/1/2/miss dosage
    2. Filter by `min_carriers_per_band: 5` (Quentin's hard filter at the
       producer level — pre-emit so the JSON stays small AND the atlas
       only sees variants that actually segregate). An SV passes if at
       least one karyotype band has ≥ 5 carriers (HET or HOM-ALT).
    3. Per group (H1/H1, H1/H2, H2/H2): compute AA/AB/BB/miss counts
    4. Fisher exact on the most-discriminating pair (H1/H1 vs H2/H2)
    5. Benjamini-Hochberg FDR across all SVs in the candidate
    6. Pattern label: deterministic decision rule on (OR, FDR, zone, group counts)
    7. Optional: a `samples_with_call` field per SV (list of carrier
       sample_ids) so the atlas can dim non-selected glyphs without
       loading the full support layer.

Usage:
    python STEP_SV_GT_AGG_aggregate_genotype_counts.py \\
        --vcf            merged_delly_manta.vcf.gz \\
        --locks          locks/cand_LG28_15Mb.json \\
        --candidate-id   cand_LG28_15Mb \\
        --chrom          C_gar_LG28 \\
        --left           15142000 \\
        --right          18124000 \\
        --out-root       data/

THIS IS A STUB. The classification, Fisher, and FDR logic must be
written against your actual VCF library (cyvcf2 or pysam). The JSON
shape it emits is what the atlas validates against and is locked.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

# Producer-side filter — Quentin's spec'd hard threshold
MIN_CARRIERS_PER_BAND = 5

# Pattern-label rules (deterministic; producer side; atlas trusts the label)
def classify_pattern(
    fdr: float,
    odds_ratio: float,
    zone: str,
    group_counts: dict[str, dict[str, int]],
) -> str:
    """
    Return the pattern_label key for this SV.

    Rules (Quentin's classification framework):
      - canonical_breakpoint_marker: at boundary zone, FDR < 0.05, OR > 20,
                                     mostly AA in REF group + BB in ALT group
      - dominant_presence_marker:    body or boundary, FDR < 0.05, OR > 5
      - het_specific_marker:         AB enriched in HET group, FDR < 0.05
      - sub_haplotype_marker:        body, FDR < 0.05, partial penetrance
      - internal_linked_marker:      body, FDR < 0.10, lower OR but still signal
      - uninformative:               FDR ≥ 0.10 or ambiguous distribution
    """
    if fdr is None or fdr >= 0.10:
        return 'uninformative'
    is_boundary = zone in ('left_boundary', 'right_boundary')
    is_body     = zone == 'inversion_body'
    if is_boundary and fdr < 0.05 and odds_ratio > 20:
        return 'canonical_breakpoint_marker'
    if fdr < 0.05 and odds_ratio > 5:
        return 'dominant_presence_marker'
    # HET specificity: |AB|/|H1H2 total| significantly higher than other groups
    h1h2 = group_counts.get('H1/H2', {})
    h1h1 = group_counts.get('H1/H1', {})
    h2h2 = group_counts.get('H2/H2', {})
    h1h2_n = sum(h1h2.values()) or 1
    if h1h2.get('AB', 0) / h1h2_n > 0.6 and \
       h1h1.get('AB', 0) / max(sum(h1h1.values()), 1) < 0.2 and \
       h2h2.get('AB', 0) / max(sum(h2h2.values()), 1) < 0.2 and \
       fdr < 0.05:
        return 'het_specific_marker'
    if is_body and fdr < 0.05:
        return 'sub_haplotype_marker'
    if is_body and fdr < 0.10:
        return 'internal_linked_marker'
    return 'uninformative'


def passes_min_carriers_filter(
    group_dosages: dict[str, list[int]],
    min_carriers: int = MIN_CARRIERS_PER_BAND,
) -> bool:
    """
    Quentin's pre-emit filter: keep an SV only if at least one karyotype
    band has ≥ `min_carriers` carriers (HET or HOM-ALT). For the 226-sample
    cohort with bands of ~60-100, n=5 is the floor for "this SV actually
    segregates with the band". Below that, it's noise.

    `group_dosages` is { 'H1/H1': [0, 0, 1, ...], 'H1/H2': [...], ... }
    where each value is a per-sample dosage (0/1/2/-1).
    """
    for grp, doses in group_dosages.items():
        n_carriers = sum(1 for d in doses if d == 1 or d == 2)
        if n_carriers >= min_carriers:
            return True
    return False


def aggregate_genotype_counts(
    sv_records: list[dict],
    locked_labels: list[dict],
    candidate: dict,
) -> dict:
    """
    Build the sv_genotype_counts_v1 payload from raw per-SV records.

    `sv_records`: list of {sv_id, sv_type, position_bp, end_bp, zone,
                            quality, callers, dosages: {sample_id: 0|1|2|-1}}
    `locked_labels`: list of {sample_id, label} (HOMO_1 / HET / HOMO_2)
    `candidate`: dict with candidate_id, chrom, boundary_left_bp,
                 boundary_right_bp, zone_definitions_bp.

    Pipeline:
      1. Build sample → group map.
      2. For each SV: split dosages by group; apply min_carriers filter.
      3. Compute genotype_counts (AA/AB/BB/miss) per group.
      4. Fisher exact + Benjamini-Hochberg FDR over the surviving SVs.
      5. Classify pattern_label.
      6. Compose boundary_summary aggregate (count by zone+sv_type).
    """
    # 1. sample → group
    label_to_group = {'HOMO_1': 'H1/H1', 'HET': 'H1/H2', 'HOMO_2': 'H2/H2'}
    sample_group = {l['sample_id']: label_to_group.get(l['label']) for l in locked_labels}
    groups = ['H1/H1', 'H1/H2', 'H2/H2']
    group_n = {g: sum(1 for v in sample_group.values() if v == g) for g in groups}

    # 2-5. iterate SVs
    surviving = []  # list of dict matching atlas sv_calls schema
    pvalues   = []
    for r in sv_records:
        # Split dosages by group
        gd = {g: [] for g in groups}
        for sid, dose in r.get('dosages', {}).items():
            g = sample_group.get(sid)
            if g is not None:
                gd[g].append(dose)
        if not passes_min_carriers_filter(gd):
            continue
        # AA/AB/BB/miss counts per group
        gc = {}
        for g in groups:
            doses = gd[g]
            gc[g] = {
                'AA':   sum(1 for d in doses if d == 0),
                'AB':   sum(1 for d in doses if d == 1),
                'BB':   sum(1 for d in doses if d == 2),
                'miss': sum(1 for d in doses if d == -1),
            }
        # Fisher H1/H1 vs H2/H2 on (AA+AB) vs (BB) — classic INV vs REF
        a, b = gc['H1/H1']['AA'] + gc['H1/H1']['AB'], gc['H1/H1']['BB']
        c, d = gc['H2/H2']['AA'] + gc['H2/H2']['AB'], gc['H2/H2']['BB']
        odds_ratio = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
        # Real implementation uses scipy.stats.fisher_exact here:
        # from scipy.stats import fisher_exact
        # _, p = fisher_exact([[a, b], [c, d]])
        # Stub uses a placeholder:
        p_value = float(_placeholder_fisher_p(a, b, c, d))
        pvalues.append(p_value)
        n_carriers_total = sum(
            1 for d in r.get('dosages', {}).values() if d in (1, 2)
        )
        carriers_list = sorted(
            sid for sid, dose in r.get('dosages', {}).items() if dose in (1, 2)
        )
        surviving.append({
            'sv_id':                  r['sv_id'],
            'sv_type':                r['sv_type'],
            'chrom':                  candidate['chrom'],
            'position_bp':            r['position_bp'],
            'end_bp':                 r.get('end_bp'),
            'zone':                   r.get('zone'),
            'distance_to_edge_bp':    r.get('distance_to_edge_bp'),
            'n_samples_with_call':    n_carriers_total,
            'samples_with_call':      carriers_list,
            'quality':                r.get('quality', 'PASS'),
            'callers':                r.get('callers', []),
            'genotype_counts':        gc,
            'fisher': {
                'comparison':         'H1/H1_vs_H2/H2',
                'odds_ratio':         odds_ratio,
                'p_value':            p_value,
                'fdr_bh':             None,  # filled in step 5 below
            },
            'pattern_label':          None,  # filled below
            'notes':                  r.get('notes', ''),
        })

    # 5. Benjamini-Hochberg FDR across surviving p-values
    fdrs = _benjamini_hochberg(pvalues)
    for sv, fdr in zip(surviving, fdrs):
        sv['fisher']['fdr_bh'] = fdr
        sv['pattern_label']    = classify_pattern(
            fdr, sv['fisher']['odds_ratio'], sv.get('zone', ''),
            sv['genotype_counts']
        )

    # 6. boundary_summary aggregate
    boundary_summary = _build_boundary_summary(surviving, candidate)

    return {
        'format_version': 'sv_genotype_counts_v1',
        'candidate_id':   candidate['candidate_id'],
        'chrom':          candidate['chrom'],
        'boundary_left_bp':  candidate['boundary_left_bp'],
        'boundary_right_bp': candidate['boundary_right_bp'],
        'zone_definitions_bp': candidate['zone_definitions_bp'],
        'groups_used':    {
            g: {'n': group_n[g], 'members': []} for g in groups
        },
        'min_carriers_filter': {
            'min_carriers_per_band': MIN_CARRIERS_PER_BAND,
            'n_sv_pre_filter':  len(sv_records),
            'n_sv_post_filter': len(surviving),
        },
        'sv_calls':         surviving,
        'boundary_summary': boundary_summary,
        'upset_top_combinations': [],  # produced separately by STEP_SV_EVID_COMB
    }


# ---------------------------------------------------------------------------
# Stat helpers (placeholders; swap with scipy in production)
# ---------------------------------------------------------------------------

def _placeholder_fisher_p(a: int, b: int, c: int, d: int) -> float:
    """
    Z-approximation to Fisher's exact, just so the stub produces sensible
    numbers. Replace with `scipy.stats.fisher_exact` in production.
    """
    n = a + b + c + d
    if n == 0:
        return 1.0
    expected_a = (a + b) * (a + c) / n
    var = max(((a + b) * (c + d) * (a + c) * (b + d)) / (n * n * (n - 1)), 1e-12) if n > 1 else 1.0
    z = (a - expected_a) / math.sqrt(var)
    # two-sided p-value via normal approximation
    p = 2 * (1 - _phi(abs(z)))
    return max(min(p, 1.0), 1e-300)


def _phi(z: float) -> float:
    """Standard normal CDF."""
    return 0.5 * (1 + math.erf(z / math.sqrt(2)))


def _benjamini_hochberg(pvalues: list[float]) -> list[float]:
    """
    Benjamini-Hochberg FDR control. Returns adjusted p-values in the same
    order as the input. Equivalent to R's p.adjust(p, method='BH').
    """
    n = len(pvalues)
    if n == 0:
        return []
    indexed = sorted(enumerate(pvalues), key=lambda t: t[1])
    fdrs = [0.0] * n
    prev = 1.0
    for rank in range(n - 1, -1, -1):
        original_idx, p = indexed[rank]
        adj = p * n / (rank + 1)
        if adj > prev:
            adj = prev
        prev = adj
        fdrs[original_idx] = min(adj, 1.0)
    return fdrs


def _build_boundary_summary(svs: list[dict], cand: dict) -> dict:
    """Build boundary_summary.{left,right}.by_sv_type.{n_total, n_associated_fdr_lt_0_05}."""
    out = {
        'left':  {'interval_bp': None, 'by_sv_type': {}},
        'right': {'interval_bp': None, 'by_sv_type': {}},
    }
    types = ['BND', 'INV', 'DEL', 'DUP', 'Other']
    z = cand.get('zone_definitions_bp', {})
    out['left']['interval_bp']  = z.get('left_boundary')
    out['right']['interval_bp'] = z.get('right_boundary')
    for side, zone_key in [('left', 'left_boundary'), ('right', 'right_boundary')]:
        for t in types:
            n_total = sum(1 for s in svs if s['zone'] == zone_key and s['sv_type'] == t)
            n_assoc = sum(1 for s in svs if s['zone'] == zone_key and s['sv_type'] == t
                          and (s['fisher']['fdr_bh'] is not None
                               and s['fisher']['fdr_bh'] < 0.05))
            out[side]['by_sv_type'][t] = {
                'n_total': n_total,
                'n_associated_fdr_lt_0_05': n_assoc,
            }
    return out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--vcf',          required=True, help='merged DELLY+Manta VCF')
    ap.add_argument('--locks',        required=True, help='locked_labels JSON')
    ap.add_argument('--candidate-id', required=True)
    ap.add_argument('--chrom',        required=True)
    ap.add_argument('--left',  type=int, required=True, help='boundary_left_bp')
    ap.add_argument('--right', type=int, required=True, help='boundary_right_bp')
    ap.add_argument('--out-root', default='data/')
    args = ap.parse_args()

    # ---- VCF parsing block — replace with cyvcf2/pysam ---------------------
    # sv_records = parse_vcf_for_region(args.vcf, args.chrom, args.left, args.right)
    # For the stub, we expect upstream to provide the records directly.
    # ------------------------------------------------------------------------
    raise SystemExit("STEP_SV_GT_AGG: VCF parsing block not implemented. "
                     "See top-of-file docstring; swap in cyvcf2 or pysam.")


if __name__ == '__main__':
    main()
