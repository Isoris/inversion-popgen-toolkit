#!/usr/bin/env python3
"""
STEP_C60_cargo_inventory.py — 6A inventory & diagnostic gate

For each inversion candidate, produces:
  - <CARGO_INVENTORY_DIR>/<cid>/genes.tsv
      gene_id, chrom, start, end, strand, length_bp, cds_length_bp,
      n_transcripts, distance_to_nearest_breakpoint_bp,
      preferred_name, family, go_terms, kegg_pathway, repeat_overlap_frac
  - <CARGO_INVENTORY_DIR>/<cid>/paralog_summary.tsv
      family, n_genes_in_inversion, n_genes_genome_wide, enrichment_ratio
  - <CARGO_INVENTORY_DIR>/diagnostic_table.tsv  (one row per candidate)
      candidate_id, chrom, start_bp, end_bp, length_bp,
      n_HOM_REF, n_HET, n_HOM_INV,
      n_genes_inside, gene_density_per_mb,
      n_missense_sites_inside, median_NS_per_gene, max_NS_per_gene,
      n_lof_sites_inside,
      level_1_ok, level_2_ok, notes

This is the gating script. Downstream cargo modules (STEP_C61, STEP_C62) read
diagnostic_table.tsv to decide which candidates to process at which level.

Usage:
    python3 STEP_C60_cargo_inventory.py [--candidate-id CID]

Reads paths from environment via 00_cargo_config.sh.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ── Resolve paths from environment ──
BASE = os.environ['BASE']
CARGO_DIR = Path(os.environ['CARGO_DIR'])
CARGO_INVENTORY_DIR = Path(os.environ['CARGO_INVENTORY_DIR'])
GENE_BED = Path(os.environ['GENE_BED'])
GENE_LENGTH_TSV = Path(os.environ['GENE_LENGTH_TSV'])
GENE_FUNCTION_TSV = Path(os.environ.get('GENE_FUNCTION_TSV', ''))
REPEAT_BED = os.environ.get('REPEAT_BED', '') or None
SNAKE_CAND_FILE = Path(os.environ.get('SNAKE_CAND_FILE', ''))
INVDIR = Path(os.environ['INVDIR'])
SAMPLE_REGISTRY = Path(os.environ['SAMPLE_REGISTRY'])
VARIANT_MASTER = Path(os.environ['VARIANT_MASTER'])

LEVEL1_MIN_HOMO = int(os.environ.get('CARGO_LEVEL1_MIN_HOMO', '5'))
LEVEL2_MIN_HOMO = int(os.environ.get('CARGO_LEVEL2_MIN_HOMO', '10'))
LEVEL2_MIN_NS = int(os.environ.get('CARGO_LEVEL2_MIN_NS_PER_GENE', '5'))


# =============================================================================
# I/O helpers
# =============================================================================

def open_maybe_gz(path: str | Path, mode: str = 'rt'):
    p = str(path)
    if p.endswith('.gz'):
        return gzip.open(p, mode)
    return open(p, mode)


def read_tsv_dicts(path: str | Path) -> List[Dict[str, str]]:
    out = []
    with open_maybe_gz(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < len(header):
                parts += [''] * (len(header) - len(parts))
            out.append(dict(zip(header, parts)))
    return out


# =============================================================================
# Load reference annotation
# =============================================================================

def load_genes_bed() -> List[Tuple[str, int, int, str, str]]:
    """Returns list of (chrom, start, end, gene_id, strand)."""
    rows = []
    with open_maybe_gz(GENE_BED) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 6:
                continue
            rows.append((parts[0], int(parts[1]), int(parts[2]),
                        parts[3], parts[5]))
    return rows


def load_gene_lengths() -> Dict[str, Dict[str, str]]:
    return {r['gene_id']: r for r in read_tsv_dicts(GENE_LENGTH_TSV)}


def load_gene_function() -> Dict[str, Dict[str, str]]:
    if not GENE_FUNCTION_TSV or not Path(GENE_FUNCTION_TSV).exists():
        return {}
    return {r['gene_id']: r for r in read_tsv_dicts(GENE_FUNCTION_TSV)}


# =============================================================================
# Repeat overlap (optional)
# =============================================================================

def load_repeat_intervals() -> Dict[str, List[Tuple[int, int]]]:
    """Returns chrom → sorted list of (start, end). Empty if no REPEAT_BED."""
    out = defaultdict(list)
    if not REPEAT_BED or not Path(REPEAT_BED).exists():
        return out
    with open_maybe_gz(REPEAT_BED) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            try:
                out[parts[0]].append((int(parts[1]), int(parts[2])))
            except ValueError:
                continue
    for chrom in out:
        out[chrom].sort()
    return out


def repeat_overlap_frac(chrom: str, start: int, end: int,
                        repeats: Dict[str, List[Tuple[int, int]]]) -> float:
    glen = end - start
    if glen <= 0:
        return 0.0
    rs = repeats.get(chrom, [])
    if not rs:
        return float('nan')
    # Binary search for first repeat with end >= start
    import bisect
    starts = [r[0] for r in rs]
    i = bisect.bisect_left(starts, start) - 1
    if i < 0:
        i = 0
    overlap = 0
    for j in range(i, len(rs)):
        rstart, rend = rs[j]
        if rstart >= end:
            break
        if rend <= start:
            continue
        overlap += min(rend, end) - max(rstart, start)
    return min(1.0, overlap / glen)


# =============================================================================
# Karyotype counts from registered sample groups
# =============================================================================

def load_karyotype_counts(cid: str) -> Tuple[int, int, int]:
    """Read sample_registry to count members of inv_<cid>_HOM_REF / HET / HOM_INV.

    Falls back to the per-candidate sample_karyotypes.tsv if groups not registered.
    """
    n = {'HOM_REF': 0, 'HET': 0, 'HOM_INV': 0}
    groups_dir = SAMPLE_REGISTRY / 'groups'
    found_any = False
    for kary in n:
        gfile = groups_dir / f'inv_{cid}_{kary}.txt'
        if gfile.exists():
            with open(gfile) as f:
                n[kary] = sum(1 for line in f if line.strip())
            found_any = True

    if not found_any:
        # Fallback: parse the per-candidate sample_karyotypes.tsv
        for cand_dir in INVDIR.glob(f'**/candidate_{cid}/data/sample_karyotypes.tsv'):
            try:
                with open(cand_dir) as f:
                    header = f.readline().rstrip('\n').split('\t')
                    kidx = header.index('karyotype') if 'karyotype' in header else None
                    if kidx is None:
                        break
                    for line in f:
                        parts = line.rstrip('\n').split('\t')
                        if len(parts) <= kidx:
                            continue
                        k = parts[kidx]
                        # Translate REF/HET/INV → HOM_REF/HET/HOM_INV
                        if k == 'REF':
                            n['HOM_REF'] += 1
                        elif k == 'INV':
                            n['HOM_INV'] += 1
                        elif k == 'HET':
                            n['HET'] += 1
                break
            except Exception as e:
                print(f"[WARN] fallback karyotype read failed for {cid}: {e}", file=sys.stderr)
    return n['HOM_REF'], n['HET'], n['HOM_INV']


# =============================================================================
# Missense / LoF density inside an interval
# =============================================================================

# SnpEff annotations we count as missense / LoF
MISSENSE_KEYWORDS = ('missense_variant',)
LOF_KEYWORDS = ('stop_gained', 'frameshift_variant', 'splice_acceptor_variant',
                'splice_donor_variant', 'start_lost', 'stop_lost')


def annotation_class(ann: str) -> str:
    a = (ann or '').lower()
    for kw in LOF_KEYWORDS:
        if kw in a:
            return 'LoF'
    if 'missense' in a:
        return 'missense'
    return 'other'


def count_variants_in_interval(variants: List[Dict[str, str]],
                                chrom: str, start: int, end: int,
                                gene_index: Dict[str, str]
                                ) -> Tuple[int, int, Dict[str, Dict[str, int]]]:
    """Returns (n_missense, n_lof, per_gene_counts).

    per_gene_counts: gene_id → {'NS': int, 'LoF': int}
    """
    n_miss = 0
    n_lof = 0
    per_gene = defaultdict(lambda: {'NS': 0, 'LoF': 0})
    for v in variants:
        c = v.get('chr', '')
        if c != chrom:
            continue
        try:
            p = int(v.get('pos', '0'))
        except ValueError:
            continue
        if not (start <= p < end):
            continue
        cls = annotation_class(v.get('snpeff_annotation', '') or v.get('bcsq_consequence', ''))
        if cls == 'missense':
            n_miss += 1
        elif cls == 'LoF':
            n_lof += 1
        else:
            continue
        gid = v.get('gene_id', '')
        if gid:
            if cls == 'missense':
                per_gene[gid]['NS'] += 1
            elif cls == 'LoF':
                per_gene[gid]['LoF'] += 1
    return n_miss, n_lof, dict(per_gene)


# =============================================================================
# Per-candidate inventory writer
# =============================================================================

def inventory_one_candidate(cand: Dict[str, str],
                             genes_by_chrom: Dict[str, List[Tuple[int, int, str, str]]],
                             gene_lengths: Dict[str, Dict[str, str]],
                             gene_function: Dict[str, Dict[str, str]],
                             repeats: Dict[str, List[Tuple[int, int]]],
                             variants: List[Dict[str, str]]) -> Dict[str, str]:
    cid = cand['candidate_id']
    chrom = cand['chrom']
    start_bp = int(float(cand['start_bp']))
    end_bp = int(float(cand['end_bp']))
    length_bp = end_bp - start_bp

    out_dir = CARGO_INVENTORY_DIR / cid
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── Genes overlapping the interval ──
    chrom_genes = genes_by_chrom.get(chrom, [])
    inside = []
    for gstart, gend, gid, strand in chrom_genes:
        # Overlap test: gene midpoint inside, or any overlap
        if gend <= start_bp or gstart >= end_bp:
            continue
        inside.append((gstart, gend, gid, strand))

    # ── Per-gene rows ──
    gene_rows = []
    for gstart, gend, gid, strand in inside:
        glen_info = gene_lengths.get(gid, {})
        gfun = gene_function.get(gid, {})
        d_left = abs(gstart - start_bp)
        d_right = abs(gend - end_bp)
        d_breakpoint = min(d_left, d_right)
        rep_frac = repeat_overlap_frac(chrom, gstart, gend, repeats)
        gene_rows.append({
            'gene_id': gid,
            'chrom': chrom,
            'start': gstart,
            'end': gend,
            'strand': strand,
            'length_bp': glen_info.get('length_bp', str(gend - gstart)),
            'cds_length_bp': glen_info.get('cds_length_bp', '0'),
            'n_transcripts': glen_info.get('n_transcripts', '0'),
            'distance_to_nearest_breakpoint_bp': d_breakpoint,
            'preferred_name': gfun.get('preferred_name', ''),
            'family': gfun.get('family', ''),
            'go_terms': gfun.get('go_terms', ''),
            'kegg_pathway': gfun.get('kegg_pathway', ''),
            'description': gfun.get('description', ''),
            'repeat_overlap_frac': '' if rep_frac != rep_frac else f"{rep_frac:.4f}",  # NaN check
        })

    genes_path = out_dir / 'genes.tsv'
    if gene_rows:
        cols = list(gene_rows[0].keys())
    else:
        cols = ['gene_id', 'chrom', 'start', 'end', 'strand', 'length_bp',
                'cds_length_bp', 'n_transcripts',
                'distance_to_nearest_breakpoint_bp',
                'preferred_name', 'family', 'go_terms', 'kegg_pathway',
                'description', 'repeat_overlap_frac']
    with open(genes_path, 'w') as out:
        out.write('\t'.join(cols) + '\n')
        for r in gene_rows:
            out.write('\t'.join(str(r[c]) for c in cols) + '\n')

    # ── Family / paralog summary inside this inversion ──
    fam_counts = defaultdict(int)
    for r in gene_rows:
        fam = r['family']
        if fam:
            fam_counts[fam] += 1
    fam_path = out_dir / 'paralog_summary.tsv'
    with open(fam_path, 'w') as out:
        out.write("family\tn_genes_in_inversion\n")
        for fam in sorted(fam_counts.keys(), key=lambda f: -fam_counts[f]):
            out.write(f"{fam}\t{fam_counts[fam]}\n")

    # ── Variants inside the interval ──
    n_miss, n_lof, per_gene_counts = count_variants_in_interval(
        variants, chrom, start_bp, end_bp,
        {r['gene_id']: r['gene_id'] for r in gene_rows})

    ns_per_gene = [per_gene_counts.get(r['gene_id'], {'NS': 0})['NS']
                   for r in gene_rows]
    median_ns = (sorted(ns_per_gene)[len(ns_per_gene) // 2]
                 if ns_per_gene else 0)
    max_ns = max(ns_per_gene) if ns_per_gene else 0

    # ── Karyotype counts ──
    n_ref, n_het, n_inv = load_karyotype_counts(cid)

    # ── Diagnostic gate ──
    level_1_ok = (n_ref >= LEVEL1_MIN_HOMO and n_inv >= LEVEL1_MIN_HOMO)
    level_2_ok = (n_ref >= LEVEL2_MIN_HOMO and n_inv >= LEVEL2_MIN_HOMO
                  and max_ns >= LEVEL2_MIN_NS)

    notes = []
    if n_ref < LEVEL1_MIN_HOMO:
        notes.append(f'n_HOM_REF<{LEVEL1_MIN_HOMO}')
    if n_inv < LEVEL1_MIN_HOMO:
        notes.append(f'n_HOM_INV<{LEVEL1_MIN_HOMO}')
    if level_1_ok and not level_2_ok:
        if n_ref < LEVEL2_MIN_HOMO or n_inv < LEVEL2_MIN_HOMO:
            notes.append(f'sample_size_below_L2_threshold')
        if max_ns < LEVEL2_MIN_NS:
            notes.append(f'max_NS_per_gene<{LEVEL2_MIN_NS}')

    return {
        'candidate_id': cid,
        'chrom': chrom,
        'start_bp': start_bp,
        'end_bp': end_bp,
        'length_bp': length_bp,
        'n_HOM_REF': n_ref,
        'n_HET': n_het,
        'n_HOM_INV': n_inv,
        'n_genes_inside': len(gene_rows),
        'gene_density_per_mb': round(len(gene_rows) / max(length_bp / 1e6, 1e-6), 3),
        'n_missense_sites_inside': n_miss,
        'median_NS_per_gene': median_ns,
        'max_NS_per_gene': max_ns,
        'n_lof_sites_inside': n_lof,
        'level_1_ok': int(level_1_ok),
        'level_2_ok': int(level_2_ok),
        'notes': ';'.join(notes) if notes else '',
    }


# =============================================================================
# Main
# =============================================================================

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--candidate-id', default=None,
                    help='Process only this candidate (default: all)')
    args = ap.parse_args()

    print(f"[STEP_C60] CARGO_INVENTORY_DIR = {CARGO_INVENTORY_DIR}")
    print(f"[STEP_C60] Loading gene BED: {GENE_BED}")
    genes = load_genes_bed()
    genes_by_chrom = defaultdict(list)
    for chrom, start, end, gid, strand in genes:
        genes_by_chrom[chrom].append((start, end, gid, strand))
    for chrom in genes_by_chrom:
        genes_by_chrom[chrom].sort()
    print(f"[STEP_C60] {sum(len(v) for v in genes_by_chrom.values())} protein-coding genes across {len(genes_by_chrom)} chromosomes")

    print(f"[STEP_C60] Loading gene lengths: {GENE_LENGTH_TSV}")
    gene_lengths = load_gene_lengths()
    print(f"[STEP_C60] Loading gene function: {GENE_FUNCTION_TSV}")
    gene_function = load_gene_function()

    print(f"[STEP_C60] Loading repeat intervals: {REPEAT_BED or '(none)'}")
    repeats = load_repeat_intervals()

    print(f"[STEP_C60] Loading variants: {VARIANT_MASTER}")
    if not VARIANT_MASTER.exists():
        print(f"[ERROR] Missing variant_master_scored.tsv at {VARIANT_MASTER}")
        sys.exit(1)
    variants = read_tsv_dicts(VARIANT_MASTER)
    print(f"[STEP_C60] {len(variants)} annotated variants loaded")

    print(f"[STEP_C60] Loading candidates: {SNAKE_CAND_FILE}")
    cands = read_tsv_dicts(SNAKE_CAND_FILE)
    if args.candidate_id:
        cands = [c for c in cands if c['candidate_id'] == args.candidate_id]
        if not cands:
            print(f"[ERROR] Candidate {args.candidate_id} not found")
            sys.exit(1)
    print(f"[STEP_C60] Processing {len(cands)} candidate(s)")

    diag_rows = []
    for i, cand in enumerate(cands, 1):
        cid = cand['candidate_id']
        print(f"  [{i}/{len(cands)}] {cid} ...")
        try:
            row = inventory_one_candidate(
                cand, genes_by_chrom, gene_lengths, gene_function,
                repeats, variants)
            diag_rows.append(row)
        except Exception as e:
            print(f"    ERROR processing {cid}: {e}")
            import traceback; traceback.print_exc()

    # ── Write diagnostic table ──
    diag_path = CARGO_INVENTORY_DIR / 'diagnostic_table.tsv'
    if diag_rows:
        cols = list(diag_rows[0].keys())
        with open(diag_path, 'w') as out:
            out.write('\t'.join(cols) + '\n')
            for r in diag_rows:
                out.write('\t'.join(str(r[c]) for c in cols) + '\n')
        print(f"\n[STEP_C60] Diagnostic table → {diag_path}")
        n_l1 = sum(1 for r in diag_rows if r['level_1_ok'])
        n_l2 = sum(1 for r in diag_rows if r['level_2_ok'])
        print(f"  Candidates analyzable at Level 1 (per-arrangement burden): {n_l1}/{len(diag_rows)}")
        print(f"  Candidates analyzable at Level 2 (configuration spectrum): {n_l2}/{len(diag_rows)}")


if __name__ == '__main__':
    main()
