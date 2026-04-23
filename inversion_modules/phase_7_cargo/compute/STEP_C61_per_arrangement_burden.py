#!/usr/bin/env python3
"""
STEP_C61_per_arrangement_burden.py — 6C Level 1

Per-candidate, per-arrangement, per-gene aggregation of evolutionary signatures.
Reuses MODULE_CONSERVATION's variant annotations (no re-annotation), subsets
genotypes by HOM_REF / HOM_INV sample groups, and writes results to the
results_registry.

Per gene per arrangement (one HOM_REF row, one HOM_INV row):
  - cargo_burden       — VESM LLR sum and mean across alt alleles carried
  - cargo_lof_count    — count of HIGH-impact variants (LoF)
  - cargo_nonsyn_synon — segregating NS and SS counts on this arrangement
                          (the beetle-paper Table S4 split, one row per arrangement)
  - cargo_sfs          — joint missense SFS (folded by arrangement freq)

Per gene as a pairwise A-vs-B comparison:
  - cargo_fst_per_site  — Hudson FST per missense site between HOM_REF and HOM_INV
                            (the beetle-paper Table S5 as continuous distribution)
  - cargo_burden_diff   — Δburden_per_arrangement with permutation p-value

Usage:
    python3 STEP_C61_per_arrangement_burden.py --candidate-id <CID>
    python3 STEP_C61_per_arrangement_burden.py --all   # iterate all level_1_ok candidates
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ── Resolve paths from environment ──
CARGO_DIR = Path(os.environ['CARGO_DIR'])
CARGO_INVENTORY_DIR = Path(os.environ['CARGO_INVENTORY_DIR'])
CARGO_BURDEN_DIR = Path(os.environ['CARGO_BURDEN_DIR'])
GENE_BED = Path(os.environ['GENE_BED'])
SNAKE_CAND_FILE = Path(os.environ.get('SNAKE_CAND_FILE', ''))
SAMPLE_REGISTRY = Path(os.environ['SAMPLE_REGISTRY'])
VARIANT_MASTER = Path(os.environ['VARIANT_MASTER'])
VESM_SCORES = Path(os.environ['VESM_SCORES'])
NORMALIZED_VCF_DIR = Path(os.environ['NORMALIZED_VCF_DIR'])

CARGO_MIN_DP = int(os.environ.get('CARGO_MIN_DP', '3'))
CARGO_PERM_N = int(os.environ.get('CARGO_PERM_N', '1000'))
CARGO_THREADS = int(os.environ.get('CARGO_THREADS', '8'))

# Registry hooks
sys.path.insert(0, str(Path(os.environ['REGISTRY_LOADER_PY']).parent))
from registry_loader import load_registry  # noqa: E402

SCRIPT_NAME = 'STEP_C61_per_arrangement_burden.py'

# =============================================================================
# Constants
# =============================================================================

LOF_KEYWORDS = ('stop_gained', 'frameshift_variant', 'splice_acceptor_variant',
                'splice_donor_variant', 'start_lost', 'stop_lost')


def annotation_class(ann: str) -> str:
    a = (ann or '').lower()
    for kw in LOF_KEYWORDS:
        if kw in a:
            return 'LoF'
    if 'missense' in a:
        return 'missense'
    if 'synonymous' in a and 'missense' not in a:
        return 'synonymous'
    return 'other'


# =============================================================================
# I/O
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


def load_sample_group(group_id: str) -> List[str]:
    p = SAMPLE_REGISTRY / 'groups' / f'{group_id}.txt'
    if not p.exists():
        raise FileNotFoundError(f'Sample group not registered: {group_id} (expected {p})')
    return [line.strip() for line in p.read_text().splitlines() if line.strip()]


# =============================================================================
# Genotype extraction from per-chromosome cohort VCFs
# =============================================================================

def extract_genotypes(chrom: str, start_bp: int, end_bp: int,
                       sample_subset: List[str]
                       ) -> Tuple[List[str], List[Dict[str, int]]]:
    """Returns (var_keys, genotypes) where genotypes is parallel to var_keys.

    Each entry in genotypes is {sample_id: dosage} with dosage in {0,1,2}.
    Samples failing depth filter or missing → not in dict.
    Only emits sites with at least one alt allele in the subset.
    """
    sample_set = set(sample_subset)
    vcf = NORMALIZED_VCF_DIR / f'{chrom}.clair3.norm.vcf.gz'
    if not vcf.exists():
        print(f"  [WARN] no VCF for {chrom}: {vcf}", file=sys.stderr)
        return [], []

    var_keys = []
    gts = []

    with gzip.open(vcf, 'rt') as f:
        vcf_samples: List[str] = []
        kept_idx: List[int] = []   # indices in vcf_samples that we care about
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                parts = line.rstrip('\n').split('\t')
                vcf_samples = parts[9:]
                kept_idx = [i for i, s in enumerate(vcf_samples) if s in sample_set]
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 10:
                continue
            c = parts[0]
            if c != chrom:
                continue
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            if pos < start_bp or pos >= end_bp:
                continue
            ref, alt = parts[3], parts[4]

            # Multiallelic handling: emit one row per alt allele
            alts = alt.split(',')
            fmt = parts[8].split(':')
            try:
                gt_idx = fmt.index('GT')
            except ValueError:
                continue
            dp_idx = fmt.index('DP') if 'DP' in fmt else -1

            for alt_n, this_alt in enumerate(alts, start=1):
                vk = f"{c}:{pos}:{ref}:{this_alt}"
                gt_dict: Dict[str, int] = {}
                for i in kept_idx:
                    sample = vcf_samples[i]
                    field = parts[9 + i]
                    fparts = field.split(':')
                    gt = fparts[gt_idx] if gt_idx < len(fparts) else './.'
                    # Depth filter
                    dp = 0
                    if dp_idx >= 0 and dp_idx < len(fparts):
                        try:
                            dp = int(fparts[dp_idx])
                        except ValueError:
                            dp = 0
                    if dp < CARGO_MIN_DP:
                        continue
                    if gt in ('./.', '.'):
                        continue
                    # Count copies of THIS specific alt allele (alt_n)
                    seps = '/|'
                    a1, a2 = '.', '.'
                    for sep in seps:
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
                    gt_dict[sample] = dosage
                if any(d > 0 for d in gt_dict.values()):
                    var_keys.append(vk)
                    gts.append(gt_dict)
    return var_keys, gts


# =============================================================================
# Annotation join
# =============================================================================

def index_annotations(variants: List[Dict[str, str]],
                      vesm: Dict[str, float]
                      ) -> Dict[str, Dict[str, str]]:
    """var_key → annotation dict including resolved 'class' and 'vesm_llr'."""
    out: Dict[str, Dict[str, str]] = {}
    for v in variants:
        vk = v.get('var_key', '')
        if not vk:
            continue
        ann_se = v.get('snpeff_annotation', '')
        ann_cq = v.get('bcsq_consequence', '')
        # Use the more severe class (max of two)
        cls = annotation_class(ann_se)
        cls_cq = annotation_class(ann_cq)
        if cls == 'other' and cls_cq != 'other':
            cls = cls_cq
        elif cls_cq == 'LoF':
            cls = 'LoF'
        out[vk] = {
            'gene_id': v.get('gene_id', '') or v.get('bcsq_gene', ''),
            'class': cls,
            'snpeff_annotation': ann_se,
            'snpeff_impact': v.get('snpeff_impact', ''),
            'vesm_llr': str(vesm.get(vk, '')),
        }
    return out


def load_vesm() -> Dict[str, float]:
    out: Dict[str, float] = {}
    if not VESM_SCORES.exists():
        print(f"  [WARN] No VESM scores file: {VESM_SCORES}", file=sys.stderr)
        return out
    for r in read_tsv_dicts(VESM_SCORES):
        try:
            out[r['var_key']] = float(r['vesm_llr'])
        except (ValueError, KeyError):
            pass
    return out


# =============================================================================
# Per-arrangement aggregation
# =============================================================================

def per_arrangement_stats(var_keys: List[str], gts: List[Dict[str, int]],
                          ann: Dict[str, Dict[str, str]],
                          n_samples: int
                          ) -> Tuple[Dict[str, Dict[str, float]], Dict[str, List[float]]]:
    """Returns (per_gene_summary, sfs_bins).

    per_gene_summary: gene_id → {
        'n_NS': int, 'n_SS': int, 'n_LoF': int,
        'burden_sum': float, 'burden_mean': float, 'n_carriers_total': int
    }
    sfs_bins: 'NS_freqs', 'SS_freqs', 'LoF_freqs' → list of frequencies (per site)
    """
    by_gene: Dict[str, Dict[str, float]] = defaultdict(
        lambda: {'n_NS': 0, 'n_SS': 0, 'n_LoF': 0,
                 'burden_sum': 0.0, 'burden_n_scored': 0, 'n_carriers_total': 0}
    )
    sfs = {'NS_freqs': [], 'SS_freqs': [], 'LoF_freqs': []}
    n_chr = max(2 * n_samples, 1)
    for vk, gt in zip(var_keys, gts):
        a = ann.get(vk)
        if a is None:
            continue
        cls = a['class']
        if cls not in ('missense', 'synonymous', 'LoF'):
            continue
        gid = a['gene_id']
        if not gid:
            continue
        ac = sum(gt.values())
        if ac == 0:
            continue
        af = ac / n_chr
        if cls == 'missense':
            by_gene[gid]['n_NS'] += 1
            sfs['NS_freqs'].append(af)
            try:
                v = float(a['vesm_llr'])
                by_gene[gid]['burden_sum'] += v * ac
                by_gene[gid]['burden_n_scored'] += 1
            except (ValueError, TypeError):
                pass
        elif cls == 'synonymous':
            by_gene[gid]['n_SS'] += 1
            sfs['SS_freqs'].append(af)
        elif cls == 'LoF':
            by_gene[gid]['n_LoF'] += 1
            sfs['LoF_freqs'].append(af)
        by_gene[gid]['n_carriers_total'] += sum(1 for d in gt.values() if d > 0)

    # Compute mean
    for gid, d in by_gene.items():
        n = d.pop('burden_n_scored')
        d['burden_mean'] = d['burden_sum'] / n if n > 0 else 0.0

    return dict(by_gene), sfs


# =============================================================================
# Per-site FST
# =============================================================================

def hudson_fst(p1: float, n1: int, p2: float, n2: int) -> float:
    """Hudson's FST estimator from allele frequencies and chromosome counts.

    Returns NaN if either side has insufficient data.
    """
    if n1 < 2 or n2 < 2:
        return float('nan')
    num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
    den = p1 * (1 - p2) + p2 * (1 - p1)
    if den <= 0:
        return float('nan')
    return num / den


def per_site_fst_table(vk_list: List[str],
                        gts_A: Dict[str, Dict[str, int]],
                        gts_B: Dict[str, Dict[str, int]],
                        ann: Dict[str, Dict[str, str]],
                        n_A: int, n_B: int) -> List[Dict[str, str]]:
    """For each missense site present in either arrangement, compute FST."""
    out = []
    nA_chr = 2 * n_A
    nB_chr = 2 * n_B
    union = set(vk_list)
    for vk in sorted(union):
        a = ann.get(vk)
        if not a or a['class'] != 'missense':
            continue
        gA = gts_A.get(vk, {})
        gB = gts_B.get(vk, {})
        acA = sum(gA.values()); acB = sum(gB.values())
        if (acA + acB) == 0:
            continue
        pA = acA / nA_chr if nA_chr > 0 else 0.0
        pB = acB / nB_chr if nB_chr > 0 else 0.0
        fst = hudson_fst(pA, nA_chr, pB, nB_chr)
        out.append({
            'var_key': vk,
            'gene_id': a['gene_id'],
            'freq_HOM_REF': f"{pA:.6f}",
            'freq_HOM_INV': f"{pB:.6f}",
            'n_chr_HOM_REF': nA_chr,
            'n_chr_HOM_INV': nB_chr,
            'fst_hudson': '' if fst != fst else f"{fst:.6f}",  # NaN check
            'vesm_llr': a['vesm_llr'],
            'snpeff_annotation': a['snpeff_annotation'],
        })
    return out


# =============================================================================
# Burden permutation test (per gene, paired)
# =============================================================================

def burden_permutation_test(burden_A: float, burden_B: float,
                              all_dosages_A: List[int], all_dosages_B: List[int],
                              vesm_per_variant: List[float],
                              n_perm: int = 1000,
                              rng_seed: int = 1) -> float:
    """One-sided p-value for |burden_A - burden_B| under label permutation.

    Inputs are per-gene flat lists: for each variant in this gene, the list of
    per-chromosome dosages on each arrangement, plus the vesm score per variant.

    For speed, we shuffle the arrangement labels of the chromosomes and recompute.
    """
    import random
    rng = random.Random(rng_seed)
    obs = abs(burden_A - burden_B)
    if not vesm_per_variant or not all_dosages_A or not all_dosages_B:
        return float('nan')
    n_A = len(all_dosages_A) // len(vesm_per_variant) if vesm_per_variant else 0
    n_B = len(all_dosages_B) // len(vesm_per_variant) if vesm_per_variant else 0
    if n_A == 0 or n_B == 0:
        return float('nan')
    # Combined per-chromosome × per-variant matrix
    combined_A = [all_dosages_A[i*n_A:(i+1)*n_A] for i in range(len(vesm_per_variant))]
    combined_B = [all_dosages_B[i*n_B:(i+1)*n_B] for i in range(len(vesm_per_variant))]
    per_var = [a + b for a, b in zip(combined_A, combined_B)]
    n_total = n_A + n_B
    n_ge = 0
    for _ in range(n_perm):
        # Shuffle labels: pick n_A indices to be "A"
        idx = list(range(n_total))
        rng.shuffle(idx)
        a_idx = set(idx[:n_A])
        bA = bB = 0.0
        for vi, llr in enumerate(vesm_per_variant):
            row = per_var[vi]
            for ci, d in enumerate(row):
                if ci in a_idx:
                    bA += llr * d
                else:
                    bB += llr * d
        if abs(bA - bB) >= obs:
            n_ge += 1
    return (n_ge + 1) / (n_perm + 1)


# =============================================================================
# Main per-candidate driver
# =============================================================================

def process_candidate(cid: str, chrom: str, start_bp: int, end_bp: int,
                       reg) -> None:
    """Compute and register all 6C Level 1 outputs for one candidate."""
    g_REF = f'inv_{cid}_HOM_REF'
    g_INV = f'inv_{cid}_HOM_INV'

    samples_REF = load_sample_group(g_REF)
    samples_INV = load_sample_group(g_INV)
    n_REF, n_INV = len(samples_REF), len(samples_INV)
    print(f"  [{cid}] HOM_REF={n_REF}  HOM_INV={n_INV}  span={chrom}:{start_bp}-{end_bp}")

    # ── Genotypes per arrangement ──
    print(f"  [{cid}] Extracting genotypes (HOM_REF)...")
    vk_A, gts_list_A = extract_genotypes(chrom, start_bp, end_bp, samples_REF)
    print(f"  [{cid}] Extracting genotypes (HOM_INV)...")
    vk_B, gts_list_B = extract_genotypes(chrom, start_bp, end_bp, samples_INV)

    # Convert to var_key → {sample: dosage} dicts for easy lookup
    gts_A = {vk: gt for vk, gt in zip(vk_A, gts_list_A)}
    gts_B = {vk: gt for vk, gt in zip(vk_B, gts_list_B)}

    # ── Annotations ──
    print(f"  [{cid}] Loading annotations...")
    variants = read_tsv_dicts(VARIANT_MASTER)
    vesm = load_vesm()
    ann = index_annotations(variants, vesm)

    # ── Per-arrangement per-gene stats ──
    summ_A, sfs_A = per_arrangement_stats(vk_A, gts_list_A, ann, n_REF)
    summ_B, sfs_B = per_arrangement_stats(vk_B, gts_list_B, ann, n_INV)

    # ── Write per-arrangement summary tables to registry ──
    out_dir_cand = CARGO_BURDEN_DIR / cid
    out_dir_cand.mkdir(parents=True, exist_ok=True)

    # Table 1: per-gene per-arrangement burden + counts
    def write_arrangement_table(group_id: str, summ: Dict[str, Dict[str, float]],
                                 stat_name: str, fname: str):
        rows = []
        for gid in sorted(summ.keys()):
            d = summ[gid]
            rows.append({
                'gene_id': gid,
                'n_NS': int(d['n_NS']),
                'n_SS': int(d['n_SS']),
                'n_LoF': int(d['n_LoF']),
                'burden_sum': f"{d['burden_sum']:.4f}",
                'burden_mean': f"{d['burden_mean']:.4f}",
                'n_carriers_total': int(d['n_carriers_total']),
            })
        if not rows:
            return
        cols = list(rows[0].keys())
        # Write via registry put_interval_summary
        try:
            reg.results.put_interval_summary(
                chrom=chrom, start_bp=start_bp, end_bp=end_bp,
                group_1=group_id, stat=stat_name,
                rows=rows, fieldnames=cols,
                source_script=SCRIPT_NAME,
                upstream_files=[str(VARIANT_MASTER), str(VESM_SCORES)],
            )
        except ValueError as e:
            # If the stat enum doesn't include cargo_*, fall back to direct write
            print(f"    [WARN] registry rejected stat={stat_name}: {e}", file=sys.stderr)
            print(f"    Falling back to direct write under {out_dir_cand}", file=sys.stderr)
            fp = out_dir_cand / f'{fname}__{group_id}.tsv'
            with open(fp, 'w') as out:
                out.write('\t'.join(cols) + '\n')
                for r in rows:
                    out.write('\t'.join(str(r[c]) for c in cols) + '\n')

    write_arrangement_table(g_REF, summ_A, 'cargo_burden', 'per_gene_summary')
    write_arrangement_table(g_INV, summ_B, 'cargo_burden', 'per_gene_summary')

    # Table 2: joint missense SFS (one row per arrangement, frequency bins)
    def write_sfs(group_id: str, sfs_dict: Dict[str, List[float]]):
        # 11-bin folded SFS: [0, 0.1), [0.1, 0.2), ..., [0.9, 1.0]
        bins = [(i / 10, (i + 1) / 10) for i in range(10)]
        rows = []
        for cls in ('NS', 'SS', 'LoF'):
            freqs = sfs_dict.get(f'{cls}_freqs', [])
            for lo, hi in bins:
                if hi == 1.0:
                    n = sum(1 for f in freqs if lo <= f <= hi)
                else:
                    n = sum(1 for f in freqs if lo <= f < hi)
                rows.append({
                    'variant_class': cls,
                    'freq_bin_lo': f"{lo:.2f}",
                    'freq_bin_hi': f"{hi:.2f}",
                    'n_sites': n,
                })
        try:
            reg.results.put_interval_summary(
                chrom=chrom, start_bp=start_bp, end_bp=end_bp,
                group_1=group_id, stat='cargo_sfs',
                rows=rows, fieldnames=['variant_class', 'freq_bin_lo', 'freq_bin_hi', 'n_sites'],
                source_script=SCRIPT_NAME,
            )
        except ValueError as e:
            print(f"    [WARN] registry rejected cargo_sfs: {e}", file=sys.stderr)
            fp = out_dir_cand / f'sfs__{group_id}.tsv'
            with open(fp, 'w') as out:
                out.write('variant_class\tfreq_bin_lo\tfreq_bin_hi\tn_sites\n')
                for r in rows:
                    out.write(f"{r['variant_class']}\t{r['freq_bin_lo']}\t{r['freq_bin_hi']}\t{r['n_sites']}\n")

    write_sfs(g_REF, sfs_A)
    write_sfs(g_INV, sfs_B)

    # Table 3: per-site FST (pairwise A vs B)
    fst_rows = per_site_fst_table(
        list(set(vk_A) | set(vk_B)),
        gts_A, gts_B, ann, n_REF, n_INV)
    if fst_rows:
        cols = list(fst_rows[0].keys())
        try:
            reg.results.put_pairwise(
                chrom=chrom, group_1=g_REF, group_2=g_INV,
                stat='cargo_fst_per_site',
                rows=fst_rows, fieldnames=cols,
                source_script=SCRIPT_NAME,
                upstream_files=[str(VARIANT_MASTER), str(VESM_SCORES)],
            )
        except (ValueError, Exception) as e:
            print(f"    [WARN] registry rejected cargo_fst_per_site: {e}", file=sys.stderr)
            fp = out_dir_cand / 'fst_per_site.tsv'
            with open(fp, 'w') as out:
                out.write('\t'.join(cols) + '\n')
                for r in fst_rows:
                    out.write('\t'.join(str(r[c]) for c in cols) + '\n')

    # Table 4: per-gene burden difference + permutation p
    common_genes = set(summ_A) & set(summ_B)
    diff_rows = []
    for gid in sorted(common_genes | set(summ_A) | set(summ_B)):
        bA = summ_A.get(gid, {}).get('burden_sum', 0.0)
        bB = summ_B.get(gid, {}).get('burden_sum', 0.0)
        nA = summ_A.get(gid, {}).get('n_NS', 0)
        nB = summ_B.get(gid, {}).get('n_NS', 0)
        # Lightweight permutation: only run when both arrangements have >0 NS
        # and a meaningful burden difference
        p_val = float('nan')
        if nA > 0 and nB > 0 and abs(bA - bB) > 1e-6:
            # Build per-variant per-chromosome dosage matrices for this gene
            gene_vk = [vk for vk in (set(vk_A) | set(vk_B))
                       if ann.get(vk, {}).get('gene_id') == gid
                       and ann.get(vk, {}).get('class') == 'missense']
            if gene_vk:
                vesm_list = []
                A_flat = []
                B_flat = []
                for vk in gene_vk:
                    try:
                        v = float(ann[vk]['vesm_llr'])
                    except (ValueError, KeyError):
                        continue
                    vesm_list.append(v)
                    gA = gts_A.get(vk, {})
                    gB = gts_B.get(vk, {})
                    # Per-chromosome dosages: for each sample, two chromosomes
                    # Approximation: use diploid dosage / 2 expansion is overkill.
                    # Treat each sample as single observation (sum of 2 chr dosages).
                    A_flat.extend([gA.get(s, 0) for s in samples_REF])
                    B_flat.extend([gB.get(s, 0) for s in samples_INV])
                if vesm_list:
                    p_val = burden_permutation_test(
                        bA, bB, A_flat, B_flat, vesm_list,
                        n_perm=CARGO_PERM_N)
        diff_rows.append({
            'gene_id': gid,
            'burden_sum_HOM_REF': f"{bA:.4f}",
            'burden_sum_HOM_INV': f"{bB:.4f}",
            'delta_burden_INV_minus_REF': f"{bB - bA:.4f}",
            'n_NS_HOM_REF': nA,
            'n_NS_HOM_INV': nB,
            'n_LoF_HOM_REF': summ_A.get(gid, {}).get('n_LoF', 0),
            'n_LoF_HOM_INV': summ_B.get(gid, {}).get('n_LoF', 0),
            'permutation_p': '' if p_val != p_val else f"{p_val:.6f}",
        })
    if diff_rows:
        cols = list(diff_rows[0].keys())
        try:
            reg.results.put_pairwise(
                chrom=chrom, group_1=g_REF, group_2=g_INV,
                stat='cargo_burden_diff',
                rows=diff_rows, fieldnames=cols,
                source_script=SCRIPT_NAME,
            )
        except (ValueError, Exception) as e:
            print(f"    [WARN] registry rejected cargo_burden_diff: {e}", file=sys.stderr)
            fp = out_dir_cand / 'burden_diff.tsv'
            with open(fp, 'w') as out:
                out.write('\t'.join(cols) + '\n')
                for r in diff_rows:
                    out.write('\t'.join(str(r[c]) for c in cols) + '\n')

    print(f"  [{cid}] Done — {len(summ_A)} genes scored on HOM_REF, {len(summ_B)} on HOM_INV, {len(fst_rows)} sites with FST")


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--candidate-id')
    g.add_argument('--all', action='store_true',
                   help='Process all candidates marked level_1_ok in diagnostic_table.tsv')
    args = ap.parse_args()

    reg = load_registry()

    if args.all:
        diag = read_tsv_dicts(CARGO_INVENTORY_DIR / 'diagnostic_table.tsv')
        targets = [(r['candidate_id'], r['chrom'],
                    int(r['start_bp']), int(r['end_bp']))
                   for r in diag if r.get('level_1_ok') == '1']
    else:
        diag = read_tsv_dicts(CARGO_INVENTORY_DIR / 'diagnostic_table.tsv')
        rows = [r for r in diag if r['candidate_id'] == args.candidate_id]
        if not rows:
            print(f"[ERROR] candidate {args.candidate_id} not in diagnostic_table.tsv — run STEP_C60 first")
            sys.exit(1)
        r = rows[0]
        if r.get('level_1_ok') != '1':
            print(f"[WARN] {args.candidate_id} fails level_1 gate: {r.get('notes')}")
            sys.exit(0)
        targets = [(r['candidate_id'], r['chrom'],
                    int(r['start_bp']), int(r['end_bp']))]

    print(f"[STEP_C61] Processing {len(targets)} candidate(s)")
    for cid, chrom, s, e in targets:
        try:
            process_candidate(cid, chrom, s, e, reg)
        except Exception as exc:
            print(f"  [ERROR] {cid}: {exc}", file=sys.stderr)
            import traceback; traceback.print_exc()


if __name__ == '__main__':
    main()
