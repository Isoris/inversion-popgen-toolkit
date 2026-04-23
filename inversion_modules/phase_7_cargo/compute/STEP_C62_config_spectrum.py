#!/usr/bin/env python3
"""
STEP_C62_config_spectrum.py — 6C Level 2 (configuration spectrum)

For each gene with sufficient missense density on each arrangement, compute:
  - pairwise mutual information between missense sites (within arrangement)
  - mean MI across pairs (the "coupling" statistic)
  - Z-score and p-value vs a frequency-preserving permutation null
  - number of distinct observed configurations (chromosomes)
  - per-arrangement entropy of configuration distribution

Two nulls:
  null_independent: permute sites independently within arrangement (preserves
                    per-site allele frequency; destroys all coupling). Tests
                    whether sites are independent given their frequencies.
  null_inversion_bg: gene-set background — compare each gene's MI to the
                     distribution of MI across all genes in the same inversion.
                     This implicitly controls for the suppressed-recombination
                     baseline that all genes inside the inversion experience.

Output stats (per gene per arrangement):
  cargo_config_mi: gene_id, n_sites_used, mean_pairwise_mi, mi_z_independent,
                   mi_p_independent, mi_rank_in_inversion, n_distinct_configs,
                   config_entropy_bits

Usage:
    python3 STEP_C62_config_spectrum.py --candidate-id <CID>
    python3 STEP_C62_config_spectrum.py --all
"""
from __future__ import annotations

import argparse
import gzip
import math
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

CARGO_DIR = Path(os.environ['CARGO_DIR'])
CARGO_INVENTORY_DIR = Path(os.environ['CARGO_INVENTORY_DIR'])
CARGO_CONFIG_DIR = Path(os.environ['CARGO_CONFIG_DIR'])
SAMPLE_REGISTRY = Path(os.environ['SAMPLE_REGISTRY'])
VARIANT_MASTER = Path(os.environ['VARIANT_MASTER'])
NORMALIZED_VCF_DIR = Path(os.environ['NORMALIZED_VCF_DIR'])

CARGO_MIN_DP = int(os.environ.get('CARGO_MIN_DP', '3'))
CARGO_PERM_N = int(os.environ.get('CARGO_PERM_N', '1000'))
CARGO_MAF_FLOOR = float(os.environ.get('CARGO_MAF_FLOOR', '0.0'))
LEVEL2_MIN_HOMO = int(os.environ.get('CARGO_LEVEL2_MIN_HOMO', '10'))
LEVEL2_MIN_NS = int(os.environ.get('CARGO_LEVEL2_MIN_NS_PER_GENE', '5'))

sys.path.insert(0, str(Path(os.environ['REGISTRY_LOADER_PY']).parent))
from registry_loader import load_registry  # noqa: E402

SCRIPT_NAME = 'STEP_C62_config_spectrum.py'

LOF_KEYWORDS = ('stop_gained', 'frameshift_variant', 'splice_acceptor_variant',
                'splice_donor_variant', 'start_lost', 'stop_lost')


# =============================================================================
# I/O helpers (duplicated from C61 to keep modules independent)
# =============================================================================

def open_maybe_gz(path: str | Path, mode: str = 'rt'):
    p = str(path)
    return gzip.open(p, mode) if p.endswith('.gz') else open(p, mode)


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


def load_sample_group(group_id: str) -> List[str]:
    p = SAMPLE_REGISTRY / 'groups' / f'{group_id}.txt'
    if not p.exists():
        raise FileNotFoundError(f'Group not registered: {group_id}')
    return [line.strip() for line in p.read_text().splitlines() if line.strip()]


def annotation_class(ann: str) -> str:
    a = (ann or '').lower()
    for kw in LOF_KEYWORDS:
        if kw in a:
            return 'LoF'
    if 'missense' in a:
        return 'missense'
    return 'other'


# =============================================================================
# Build per-gene per-arrangement genotype matrices
# =============================================================================

def build_gene_matrices(chrom: str, start_bp: int, end_bp: int,
                         samples: List[str],
                         missense_var_keys_by_gene: Dict[str, List[str]]
                         ) -> Dict[str, np.ndarray]:
    """For each gene, return a (n_samples × n_sites) dosage matrix.

    Missing genotypes are encoded as -1 and later ignored by pairwise MI.
    """
    sample_set = set(samples)
    s_to_idx = {s: i for i, s in enumerate(samples)}
    vcf = NORMALIZED_VCF_DIR / f'{chrom}.clair3.norm.vcf.gz'
    if not vcf.exists():
        return {}

    # All target var_keys we want
    target_vk_to_gene = {}
    for gid, vks in missense_var_keys_by_gene.items():
        for vk in vks:
            target_vk_to_gene[vk] = gid

    # Column indices per gene
    gene_columns: Dict[str, List[str]] = defaultdict(list)
    for gid, vks in missense_var_keys_by_gene.items():
        gene_columns[gid] = list(vks)
    gene_col_idx: Dict[str, Dict[str, int]] = {}
    for gid, vks in gene_columns.items():
        gene_col_idx[gid] = {vk: i for i, vk in enumerate(vks)}

    # Allocate matrices
    matrices: Dict[str, np.ndarray] = {}
    for gid, vks in gene_columns.items():
        if not vks:
            continue
        m = np.full((len(samples), len(vks)), -1, dtype=np.int8)
        matrices[gid] = m

    with gzip.open(vcf, 'rt') as f:
        vcf_samples: List[str] = []
        kept_idx: List[int] = []
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                parts = line.rstrip('\n').split('\t')
                vcf_samples = parts[9:]
                kept_idx = [(i, s_to_idx[s])
                            for i, s in enumerate(vcf_samples)
                            if s in sample_set]
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 10:
                continue
            c, pos = parts[0], parts[1]
            if c != chrom:
                continue
            try:
                pos_i = int(pos)
            except ValueError:
                continue
            if pos_i < start_bp or pos_i >= end_bp:
                continue
            ref, alts_str = parts[3], parts[4]
            alts = alts_str.split(',')
            fmt = parts[8].split(':')
            try:
                gt_idx = fmt.index('GT')
            except ValueError:
                continue
            dp_idx = fmt.index('DP') if 'DP' in fmt else -1
            for alt_n, this_alt in enumerate(alts, start=1):
                vk = f"{c}:{pos_i}:{ref}:{this_alt}"
                gid = target_vk_to_gene.get(vk)
                if gid is None:
                    continue
                m = matrices.get(gid)
                if m is None:
                    continue
                col = gene_col_idx[gid][vk]
                for vcf_i, sample_row in kept_idx:
                    field = parts[9 + vcf_i]
                    fparts = field.split(':')
                    gt = fparts[gt_idx] if gt_idx < len(fparts) else './.'
                    dp = 0
                    if dp_idx >= 0 and dp_idx < len(fparts):
                        try:
                            dp = int(fparts[dp_idx])
                        except ValueError:
                            dp = 0
                    if dp < CARGO_MIN_DP or gt in ('./.', '.'):
                        continue
                    a1, a2 = '.', '.'
                    for sep in '/|':
                        if sep in gt:
                            try:
                                a1, a2 = gt.split(sep)[:2]
                            except Exception:
                                pass
                            break
                    if a1 == '.' or a2 == '.':
                        continue
                    try:
                        a1i, a2i = int(a1), int(a2)
                    except ValueError:
                        continue
                    dosage = (1 if a1i == alt_n else 0) + (1 if a2i == alt_n else 0)
                    m[sample_row, col] = dosage
    return matrices


# =============================================================================
# Mutual information statistics
# =============================================================================

def site_freq(col: np.ndarray) -> float:
    """Allele frequency of column, ignoring missing (-1)."""
    valid = col[col >= 0]
    if valid.size == 0:
        return 0.0
    return float(valid.sum()) / (2 * valid.size)


def maf_filter(matrix: np.ndarray, maf_floor: float) -> np.ndarray:
    """Drop columns with MAF below floor. Returns subset matrix."""
    if maf_floor <= 0:
        return matrix
    keep = []
    for j in range(matrix.shape[1]):
        f = site_freq(matrix[:, j])
        f = min(f, 1 - f)
        if f >= maf_floor:
            keep.append(j)
    return matrix[:, keep] if keep else matrix[:, :0]


def pairwise_mutual_information(matrix: np.ndarray) -> float:
    """Mean pairwise MI across columns (sites), excluding self-pairs.

    Treats dosage 0/1/2 as 3-state categorical. Missing (-1) entries are
    excluded pair-wise.
    Returns NaN if fewer than 2 valid columns.
    """
    n_samples, n_sites = matrix.shape
    if n_sites < 2:
        return float('nan')
    mi_vals = []
    for i in range(n_sites):
        for j in range(i + 1, n_sites):
            ci = matrix[:, i]
            cj = matrix[:, j]
            mask = (ci >= 0) & (cj >= 0)
            if mask.sum() < 4:
                continue
            xi = ci[mask]; xj = cj[mask]
            mi_vals.append(_mi_pair(xi, xj))
    if not mi_vals:
        return float('nan')
    return float(np.mean(mi_vals))


def _mi_pair(x: np.ndarray, y: np.ndarray) -> float:
    """Discrete MI for two categorical vectors. Uses natural log → bits/log2."""
    n = x.size
    # Joint and marginal counts
    joint = defaultdict(int)
    mx = defaultdict(int); my = defaultdict(int)
    for a, b in zip(x.tolist(), y.tolist()):
        joint[(a, b)] += 1
        mx[a] += 1; my[b] += 1
    h = 0.0
    for (a, b), c in joint.items():
        pxy = c / n
        px = mx[a] / n
        py = my[b] / n
        if pxy > 0 and px > 0 and py > 0:
            h += pxy * math.log2(pxy / (px * py))
    return h


def n_distinct_configurations(matrix: np.ndarray) -> int:
    """Count distinct rows (treating missing as a separate category).

    For configuration spectrum, this gives the number of distinct 'protein
    haplotypes' observed at the sample-pair (chromosome) level.
    """
    if matrix.size == 0:
        return 0
    rows = {tuple(r.tolist()) for r in matrix}
    return len(rows)


def configuration_entropy(matrix: np.ndarray) -> float:
    """Shannon entropy (bits) of the configuration distribution across rows."""
    if matrix.size == 0:
        return float('nan')
    counts: Dict[tuple, int] = defaultdict(int)
    for r in matrix:
        counts[tuple(r.tolist())] += 1
    n = sum(counts.values())
    h = 0.0
    for c in counts.values():
        p = c / n
        if p > 0:
            h -= p * math.log2(p)
    return h


# =============================================================================
# Permutation null
# =============================================================================

def permutation_mi_null(matrix: np.ndarray, n_perm: int,
                          rng_seed: int = 1) -> Tuple[float, float, float]:
    """Returns (z, p_one_sided, mean_null) for mean pairwise MI under
    site-frequency-preserving permutation.

    Each column is independently shuffled, preserving the per-site allele
    distribution but destroying all coupling. p is the fraction of permutations
    with null MI >= observed.
    """
    obs = pairwise_mutual_information(matrix)
    if obs != obs:  # NaN
        return (float('nan'),) * 3
    rng = np.random.default_rng(rng_seed)
    null_vals = []
    for _ in range(n_perm):
        perm = matrix.copy()
        for j in range(perm.shape[1]):
            rng.shuffle(perm[:, j])
        null_vals.append(pairwise_mutual_information(perm))
    null_arr = np.array([v for v in null_vals if v == v])
    if null_arr.size < 2:
        return (float('nan'),) * 3
    mean = float(np.mean(null_arr))
    sd = float(np.std(null_arr, ddof=1))
    z = (obs - mean) / sd if sd > 0 else float('nan')
    p = (np.sum(null_arr >= obs) + 1) / (null_arr.size + 1)
    return z, float(p), mean


# =============================================================================
# Per-candidate driver
# =============================================================================

def process_candidate(cid: str, chrom: str, start_bp: int, end_bp: int,
                       reg) -> None:
    g_REF = f'inv_{cid}_HOM_REF'
    g_INV = f'inv_{cid}_HOM_INV'
    samples_REF = load_sample_group(g_REF)
    samples_INV = load_sample_group(g_INV)
    n_REF, n_INV = len(samples_REF), len(samples_INV)
    if n_REF < LEVEL2_MIN_HOMO or n_INV < LEVEL2_MIN_HOMO:
        print(f"  [{cid}] SKIP — below Level 2 sample threshold ({n_REF},{n_INV})")
        return

    # ── Identify missense variants per gene inside the candidate ──
    print(f"  [{cid}] Loading variant master to find missense per gene...")
    variants = read_tsv_dicts(VARIANT_MASTER)
    miss_by_gene: Dict[str, List[str]] = defaultdict(list)
    for v in variants:
        if v.get('chr') != chrom:
            continue
        try:
            p = int(v.get('pos', '0'))
        except ValueError:
            continue
        if not (start_bp <= p < end_bp):
            continue
        ann = v.get('snpeff_annotation', '') or v.get('bcsq_consequence', '')
        if annotation_class(ann) != 'missense':
            continue
        gid = v.get('gene_id', '') or v.get('bcsq_gene', '')
        if not gid:
            continue
        miss_by_gene[gid].append(v['var_key'])

    # Filter genes with enough sites
    miss_by_gene = {gid: vks for gid, vks in miss_by_gene.items()
                    if len(vks) >= LEVEL2_MIN_NS}
    if not miss_by_gene:
        print(f"  [{cid}] No genes with >= {LEVEL2_MIN_NS} missense sites")
        return
    print(f"  [{cid}] {len(miss_by_gene)} genes meet missense density threshold")

    # ── Build per-gene matrices on each arrangement ──
    print(f"  [{cid}] Building genotype matrices (HOM_REF)...")
    mats_A = build_gene_matrices(chrom, start_bp, end_bp, samples_REF, miss_by_gene)
    print(f"  [{cid}] Building genotype matrices (HOM_INV)...")
    mats_B = build_gene_matrices(chrom, start_bp, end_bp, samples_INV, miss_by_gene)

    # ── Compute MI + null per gene per arrangement ──
    def gene_stats(matrix: np.ndarray, gid: str, group: str) -> Dict[str, str]:
        mat = maf_filter(matrix, CARGO_MAF_FLOOR)
        n_sites = mat.shape[1]
        if n_sites < 2:
            return {'gene_id': gid, 'n_sites_used': n_sites,
                    'mean_pairwise_mi': '', 'mi_z_independent': '',
                    'mi_p_independent': '', 'n_distinct_configs': '',
                    'config_entropy_bits': ''}
        z, p, mean_null = permutation_mi_null(mat, CARGO_PERM_N, rng_seed=hash((gid, group)) % (2**31))
        obs = pairwise_mutual_information(mat)
        n_distinct = n_distinct_configurations(mat)
        ent = configuration_entropy(mat)
        return {
            'gene_id': gid,
            'n_sites_used': n_sites,
            'mean_pairwise_mi': '' if obs != obs else f"{obs:.6f}",
            'mi_mean_null': '' if mean_null != mean_null else f"{mean_null:.6f}",
            'mi_z_independent': '' if z != z else f"{z:.4f}",
            'mi_p_independent': '' if p != p else f"{p:.6f}",
            'n_distinct_configs': n_distinct,
            'config_entropy_bits': '' if ent != ent else f"{ent:.4f}",
        }

    rows_A = [gene_stats(mats_A[gid], gid, g_REF)
              for gid in sorted(miss_by_gene) if gid in mats_A]
    rows_B = [gene_stats(mats_B[gid], gid, g_INV)
              for gid in sorted(miss_by_gene) if gid in mats_B]

    # ── Add inversion-background rank: where does each gene's MI sit
    # relative to all other genes in the SAME inversion on the same arrangement?
    def add_bg_rank(rows):
        mis = [(r['gene_id'], float(r['mean_pairwise_mi']))
               for r in rows
               if r['mean_pairwise_mi'] not in ('', None)]
        mis.sort(key=lambda x: -x[1])
        rank = {gid: i + 1 for i, (gid, _) in enumerate(mis)}
        n = len(mis)
        for r in rows:
            gid = r['gene_id']
            r['mi_rank_in_inversion'] = str(rank.get(gid, ''))
            r['mi_total_genes_in_inversion_with_mi'] = str(n)

    add_bg_rank(rows_A)
    add_bg_rank(rows_B)

    # ── Write per-arrangement to registry ──
    out_dir_cand = CARGO_CONFIG_DIR / cid
    out_dir_cand.mkdir(parents=True, exist_ok=True)

    def write_table(group_id: str, rows: List[Dict[str, str]]):
        if not rows:
            return
        cols = list(rows[0].keys())
        try:
            reg.results.put_interval_summary(
                chrom=chrom, start_bp=start_bp, end_bp=end_bp,
                group_1=group_id, stat='cargo_config_mi',
                rows=rows, fieldnames=cols,
                source_script=SCRIPT_NAME,
                upstream_files=[str(VARIANT_MASTER)],
            )
        except (ValueError, Exception) as e:
            print(f"    [WARN] registry rejected cargo_config_mi: {e}", file=sys.stderr)
            fp = out_dir_cand / f'config_mi__{group_id}.tsv'
            with open(fp, 'w') as out:
                out.write('\t'.join(cols) + '\n')
                for r in rows:
                    out.write('\t'.join(str(r[c]) for c in cols) + '\n')

    write_table(g_REF, rows_A)
    write_table(g_INV, rows_B)
    print(f"  [{cid}] Done — {len(rows_A)} genes scored on HOM_REF, {len(rows_B)} on HOM_INV")


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--candidate-id')
    g.add_argument('--all', action='store_true')
    args = ap.parse_args()

    reg = load_registry()
    diag = read_tsv_dicts(CARGO_INVENTORY_DIR / 'diagnostic_table.tsv')

    if args.all:
        targets = [(r['candidate_id'], r['chrom'],
                    int(r['start_bp']), int(r['end_bp']))
                   for r in diag if r.get('level_2_ok') == '1']
    else:
        rows = [r for r in diag if r['candidate_id'] == args.candidate_id]
        if not rows:
            print(f"[ERROR] candidate {args.candidate_id} not in diagnostic_table.tsv")
            sys.exit(1)
        r = rows[0]
        if r.get('level_2_ok') != '1':
            print(f"[SKIP] {args.candidate_id} fails level_2 gate: {r.get('notes')}")
            sys.exit(0)
        targets = [(r['candidate_id'], r['chrom'],
                    int(r['start_bp']), int(r['end_bp']))]

    print(f"[STEP_C62] Processing {len(targets)} candidate(s)")
    for cid, chrom, s, e in targets:
        try:
            process_candidate(cid, chrom, s, e, reg)
        except Exception as exc:
            print(f"  [ERROR] {cid}: {exc}", file=sys.stderr)
            import traceback; traceback.print_exc()


if __name__ == '__main__':
    main()
