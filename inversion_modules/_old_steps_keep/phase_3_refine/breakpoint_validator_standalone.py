#!/usr/bin/env python3
"""
breakpoint_validator.py — Inversion Breakpoint Support Validator
================================================================

For each candidate inversion breakpoint pair (BP1, BP2), extracts per-sample
discordant read-pair and split-read evidence from BAM files, classifies
samples into support/no-support, then performs:

  Test 1: Fisher's exact test (INV vs non-INV enrichment)
  Test 2: Chi-square / Fisher RxC (3-group: REF/HET/INV)
  Test 3: Cochran-Armitage trend test (ordered REF < HET < INV)

Outputs:
  - Per-sample evidence table (TSV)
  - Statistical test results (TSV + text report)
  - Publication-quality breakpoint evidence figure (PDF/PNG)

Usage (real data):
  python3 breakpoint_validator.py \\
    --bam_dir /path/to/markdup_bams \\
    --sample_groups sample_groups.tsv \\
    --bp1_chrom C_gar_LG12 --bp1_pos 52400000 \\
    --bp2_chrom C_gar_LG12 --bp2_pos 54600000 \\
    --ref /path/to/ref.fa \\
    --outdir results/

Usage (fake test data):
  python3 breakpoint_validator.py --generate_test_data \\
    --outdir test_results/ --n_per_group 20

Author: Quentin / Claude collaborative pipeline
"""

import argparse
import gzip
import math
import os
import random
import struct
import sys
from collections import defaultdict
from itertools import combinations

# ═══════════════════════════════════════════════════════════════════════════════
# ARGUMENT PARSING
# ═══════════════════════════════════════════════════════════════════════════════

def parse_args():
    p = argparse.ArgumentParser(
        description="Inversion breakpoint support validator with statistical testing"
    )

    # --- Mode ---
    p.add_argument("--generate_test_data", action="store_true",
                   help="Generate fake BAMs + run full pipeline on test data")
    p.add_argument("--n_per_group", type=int, default=20,
                   help="Samples per group for test data (default: 20)")

    # --- Real data inputs ---
    p.add_argument("--bam_dir", default="",
                   help="Directory with {sample}.markdup.bam files")
    p.add_argument("--sample_groups", default="",
                   help="TSV: sample_id<tab>group (REF/HET/INV)")
    p.add_argument("--bp1_chrom", default="C_gar_LG12")
    p.add_argument("--bp1_pos", type=int, default=52400000)
    p.add_argument("--bp2_chrom", default="C_gar_LG12")
    p.add_argument("--bp2_pos", type=int, default=54600000)
    p.add_argument("--ref", default="", help="Reference FASTA (for pysam)")

    # --- Parameters ---
    p.add_argument("--bp_window", type=int, default=300,
                   help="Window ±bp around breakpoint for evidence search (default: 300)")
    p.add_argument("--min_mapq", type=int, default=20,
                   help="Minimum MAPQ for supporting reads (default: 20)")
    p.add_argument("--min_clip_len", type=int, default=10,
                   help="Minimum soft-clip length for split-read evidence (default: 10)")

    # --- Output ---
    p.add_argument("--outdir", default="breakpoint_validation",
                   help="Output directory")
    p.add_argument("--label", default="candidate_INV",
                   help="Label for output files")

    return p.parse_args()


# ═══════════════════════════════════════════════════════════════════════════════
# STATISTICAL TESTS (pure Python — no scipy dependency)
# ═══════════════════════════════════════════════════════════════════════════════

def _log_factorial(n):
    """Log of n! using Stirling for large n."""
    if n <= 1:
        return 0.0
    return sum(math.log(i) for i in range(2, n + 1))


def fisher_exact_2x2(a, b, c, d):
    """
    Fisher's exact test for 2x2 contingency table.

         support_yes  support_no
    INV       a           b
    non-INV   c           d

    Returns: odds_ratio, p_value (two-sided)
    """
    n = a + b + c + d
    # Odds ratio
    if b * c == 0:
        odds_ratio = float('inf') if a * d > 0 else 0.0
    else:
        odds_ratio = (a * d) / (b * c)

    # Hypergeometric probability for a specific table
    def log_phyper(aa, bb, cc, dd):
        nn = aa + bb + cc + dd
        return (_log_factorial(aa + bb) + _log_factorial(cc + dd) +
                _log_factorial(aa + cc) + _log_factorial(bb + dd) -
                _log_factorial(nn) - _log_factorial(aa) -
                _log_factorial(bb) - _log_factorial(cc) - _log_factorial(dd))

    # Observed probability
    log_p_obs = log_phyper(a, b, c, d)
    p_obs = math.exp(log_p_obs)

    # Two-sided: sum probabilities of all tables with P <= P_obs
    row1 = a + b
    row2 = c + d
    col1 = a + c

    p_value = 0.0
    for aa in range(min(row1, col1) + 1):
        bb = row1 - aa
        cc = col1 - aa
        dd = row2 - cc
        if bb < 0 or cc < 0 or dd < 0:
            continue
        lp = log_phyper(aa, bb, cc, dd)
        if lp <= log_p_obs + 1e-10:  # tolerance for floating point
            p_value += math.exp(lp)

    p_value = min(p_value, 1.0)
    return odds_ratio, p_value


def fisher_exact_ci(a, b, c, d, alpha=0.05):
    """Approximate 95% CI for odds ratio using log transform."""
    if b * c == 0 or a * d == 0:
        return (0.0, float('inf'))
    log_or = math.log((a * d) / (b * c))
    se = math.sqrt(1.0/max(a,0.5) + 1.0/max(b,0.5) + 1.0/max(c,0.5) + 1.0/max(d,0.5))
    z = 1.96  # ~95%
    lo = math.exp(log_or - z * se)
    hi = math.exp(log_or + z * se)
    return (lo, hi)


def chi_square_rxc(table):
    """
    Chi-square test of independence for RxC table.
    table: list of lists, table[i][j] = count for row i, col j
    Returns: chi2_stat, p_value_approx, df
    """
    R = len(table)
    C = len(table[0])
    row_sums = [sum(table[i]) for i in range(R)]
    col_sums = [sum(table[i][j] for i in range(R)) for j in range(C)]
    n = sum(row_sums)

    if n == 0:
        return 0.0, 1.0, (R-1)*(C-1)

    chi2 = 0.0
    for i in range(R):
        for j in range(C):
            expected = row_sums[i] * col_sums[j] / n
            if expected > 0:
                chi2 += (table[i][j] - expected) ** 2 / expected

    df = (R - 1) * (C - 1)
    # Approximate p-value using Wilson-Hilferty chi-square → normal approx
    if df > 0:
        z = (chi2 / df) ** (1/3) - (1 - 2 / (9 * df))
        z /= math.sqrt(2 / (9 * df))
        # Normal CDF approximation
        p_value = 0.5 * math.erfc(z / math.sqrt(2))
    else:
        p_value = 1.0

    return chi2, p_value, df


def cochran_armitage_trend(table_2col, scores=None):
    """
    Cochran-Armitage trend test for ordered groups.

    table_2col: [[yes0, no0], [yes1, no1], [yes2, no2], ...]
      rows = groups in order (REF=0, HET=1, INV=2)
      cols = [support_yes, support_no]

    scores: integer scores for each group (default: 0, 1, 2, ...)

    Returns: Z_statistic, p_value (two-sided)
    """
    K = len(table_2col)
    if scores is None:
        scores = list(range(K))

    n = [table_2col[i][0] + table_2col[i][1] for i in range(K)]
    r = [table_2col[i][0] for i in range(K)]  # successes per group
    N = sum(n)
    R = sum(r)

    if N == 0 or R == 0 or R == N:
        return 0.0, 1.0

    p_hat = R / N
    t_bar = sum(scores[i] * n[i] for i in range(K)) / N

    # Numerator
    T_num = sum(scores[i] * (r[i] - n[i] * p_hat) for i in range(K))

    # Denominator
    T_den_sq = p_hat * (1 - p_hat) * (
        sum(scores[i]**2 * n[i] for i in range(K)) - N * t_bar**2
    )

    if T_den_sq <= 0:
        return 0.0, 1.0

    Z = T_num / math.sqrt(T_den_sq)
    # Two-sided p-value from normal
    p_value = 2 * 0.5 * math.erfc(abs(Z) / math.sqrt(2))

    return Z, p_value


# ═══════════════════════════════════════════════════════════════════════════════
# BAM EVIDENCE EXTRACTION (using pysam if available, else samtools CLI)
# ═══════════════════════════════════════════════════════════════════════════════

def extract_evidence_pysam(bam_path, chrom, pos, window, min_mapq, min_clip_len,
                           partner_chrom, partner_pos):
    """
    Extract breakpoint evidence from a BAM file using pysam.

    For an inversion breakpoint, we look for:
    1. Discordant pairs: both mates map to same chrom but in unexpected
       orientations (FF or RR instead of FR)
    2. Split reads: reads with SA tag or large soft clips near breakpoint
    3. Reads with mate mapping near the partner breakpoint

    Returns dict with evidence counts.
    """
    import pysam

    evidence = {
        'n_discordant_FF': 0,  # both forward (inversion signal)
        'n_discordant_RR': 0,  # both reverse (inversion signal)
        'n_split_reads': 0,
        'n_soft_clip_left': 0,
        'n_soft_clip_right': 0,
        'n_mate_at_partner': 0,
        'n_total_reads': 0,
        'supporting_read_names': set(),
    }

    start = max(0, pos - window)
    end = pos + window

    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        print(f"  WARNING: Cannot open {bam_path}: {e}", file=sys.stderr)
        return evidence

    try:
        for read in bam.fetch(chrom, start, end):
            if read.is_unmapped or read.is_secondary or read.is_duplicate:
                continue
            if read.mapping_quality < min_mapq:
                continue

            evidence['n_total_reads'] += 1

            # --- Discordant orientation (FF or RR) ---
            if read.is_paired and not read.mate_is_unmapped:
                read_fwd = not read.is_reverse
                mate_fwd = not read.mate_is_reverse

                # Same orientation = inversion signal
                if read_fwd and mate_fwd:
                    evidence['n_discordant_FF'] += 1
                    evidence['supporting_read_names'].add(read.query_name)
                elif not read_fwd and not mate_fwd:
                    evidence['n_discordant_RR'] += 1
                    evidence['supporting_read_names'].add(read.query_name)

                # Mate near partner breakpoint?
                if (read.next_reference_name == partner_chrom and
                    abs(read.next_reference_start - partner_pos) <= window * 3):
                    evidence['n_mate_at_partner'] += 1
                    evidence['supporting_read_names'].add(read.query_name)

            # --- Split reads (SA tag) ---
            if read.has_tag('SA'):
                sa = read.get_tag('SA')
                # SA format: chr,pos,strand,CIGAR,mapQ,NM;
                for sa_entry in sa.split(';'):
                    if not sa_entry.strip():
                        continue
                    parts = sa_entry.strip().split(',')
                    if len(parts) >= 5:
                        sa_chr = parts[0]
                        sa_pos = int(parts[1])
                        sa_mapq = int(parts[4])
                        if sa_mapq >= min_mapq:
                            # SA near partner breakpoint?
                            if (sa_chr == partner_chrom and
                                abs(sa_pos - partner_pos) <= window * 3):
                                evidence['n_split_reads'] += 1
                                evidence['supporting_read_names'].add(read.query_name)

            # --- Soft clips ---
            cigar = read.cigartuples
            if cigar:
                # Left soft clip
                if cigar[0][0] == 4 and cigar[0][1] >= min_clip_len:
                    if abs(read.reference_start - pos) <= window:
                        evidence['n_soft_clip_left'] += 1
                        evidence['supporting_read_names'].add(read.query_name)
                # Right soft clip
                if cigar[-1][0] == 4 and cigar[-1][1] >= min_clip_len:
                    if abs(read.reference_end - pos) <= window:
                        evidence['n_soft_clip_right'] += 1
                        evidence['supporting_read_names'].add(read.query_name)

    except ValueError as e:
        print(f"  WARNING: fetch error for {chrom}:{start}-{end}: {e}", file=sys.stderr)

    bam.close()
    return evidence


def classify_support(ev_bp1, ev_bp2, min_discordant=1, min_split=0):
    """
    Classify a sample as support=yes/no based on evidence at both breakpoints.

    A sample has support if:
    - At least min_discordant discordant pairs (FF or RR) at EITHER breakpoint
      with mate near partner, OR
    - At least min_split split reads with SA near partner, OR
    - Shared read names between BP1 and BP2 evidence
    """
    disc_bp1 = ev_bp1['n_discordant_FF'] + ev_bp1['n_discordant_RR']
    disc_bp2 = ev_bp2['n_discordant_FF'] + ev_bp2['n_discordant_RR']
    split_bp1 = ev_bp1['n_split_reads']
    split_bp2 = ev_bp2['n_split_reads']
    mate_bp1 = ev_bp1['n_mate_at_partner']
    mate_bp2 = ev_bp2['n_mate_at_partner']
    clip_bp1 = ev_bp1['n_soft_clip_left'] + ev_bp1['n_soft_clip_right']
    clip_bp2 = ev_bp2['n_soft_clip_left'] + ev_bp2['n_soft_clip_right']

    # Shared read names spanning both breakpoints
    shared = ev_bp1['supporting_read_names'] & ev_bp2['supporting_read_names']

    total_disc = disc_bp1 + disc_bp2
    total_split = split_bp1 + split_bp2
    total_mate = mate_bp1 + mate_bp2
    total_clip = clip_bp1 + clip_bp2

    # Support criteria (generous — at least 1 piece of evidence)
    has_discordant = total_disc >= min_discordant
    has_split = total_split >= min_split if min_split > 0 else False
    has_mate_bridge = total_mate >= 1
    has_shared_reads = len(shared) >= 1
    has_clips = total_clip >= 2

    support = has_discordant or has_split or has_mate_bridge or has_shared_reads

    score = total_disc + total_split * 2 + total_mate + len(shared) * 3 + total_clip

    return support, score, {
        'disc_bp1': disc_bp1, 'disc_bp2': disc_bp2,
        'split_bp1': split_bp1, 'split_bp2': split_bp2,
        'mate_bp1': mate_bp1, 'mate_bp2': mate_bp2,
        'clip_bp1': clip_bp1, 'clip_bp2': clip_bp2,
        'shared_reads': len(shared),
        'total_score': score,
    }


# ═══════════════════════════════════════════════════════════════════════════════
# FAKE BAM GENERATOR (for testing without real data)
# ═══════════════════════════════════════════════════════════════════════════════

def generate_fake_bam_data(outdir, n_per_group=20,
                           bp1_chrom="C_gar_LG12", bp1_pos=52400000,
                           bp2_chrom="C_gar_LG12", bp2_pos=54600000):
    """
    Generate fake per-sample evidence tables that simulate what we'd get
    from real BAM extraction. Creates a sample_groups.tsv and a
    simulated evidence table directly (no actual BAM files needed).

    REF samples: ~5% have spurious support (noise)
    HET samples: ~50% have moderate support (one haplotype inverted)
    INV samples: ~85% have strong support (both haplotypes inverted)
    """
    random.seed(42)
    os.makedirs(outdir, exist_ok=True)

    groups = {}
    evidence_rows = []

    for group, n, support_prob, score_range in [
        ("REF", n_per_group, 0.05, (0, 2)),
        ("HET", n_per_group, 0.55, (1, 8)),
        ("INV", n_per_group, 0.88, (3, 15)),
    ]:
        for i in range(n):
            sid = f"CGA{group[0]}{i+1:03d}"
            groups[sid] = group

            has_support = random.random() < support_prob

            if has_support:
                lo, hi = score_range
                score = random.randint(max(1, lo), hi)
                disc = random.randint(1, max(1, score // 2))
                split = random.randint(0, max(1, score // 3))
                mate = random.randint(0, max(1, score // 4))
                clips = random.randint(0, 3)
                shared = random.randint(0, min(2, disc))
            else:
                score = random.randint(0, 1)
                disc = random.randint(0, 1) if score > 0 else 0
                split = 0
                mate = 0
                clips = random.randint(0, 1)
                shared = 0

            # Split evidence between BP1 and BP2
            d1 = random.randint(0, disc)
            d2 = disc - d1
            s1 = random.randint(0, split)
            s2 = split - s1
            m1 = random.randint(0, mate)
            m2 = mate - m1
            c1 = random.randint(0, clips)
            c2 = clips - c1

            evidence_rows.append({
                'sample': sid, 'group': group,
                'support': 'yes' if has_support else 'no',
                'disc_bp1': d1, 'disc_bp2': d2,
                'split_bp1': s1, 'split_bp2': s2,
                'mate_bp1': m1, 'mate_bp2': m2,
                'clip_bp1': c1, 'clip_bp2': c2,
                'shared_reads': shared,
                'total_score': score,
                'n_total_reads_bp1': random.randint(30, 80),
                'n_total_reads_bp2': random.randint(30, 80),
            })

    # Write sample groups
    groups_file = os.path.join(outdir, "sample_groups.tsv")
    with open(groups_file, 'w') as f:
        f.write("sample_id\tgroup\n")
        for sid in sorted(groups.keys()):
            f.write(f"{sid}\t{groups[sid]}\n")

    # Write evidence table
    ev_file = os.path.join(outdir, "per_sample_breakpoint_evidence.tsv")
    cols = ['sample', 'group', 'support', 'disc_bp1', 'disc_bp2',
            'split_bp1', 'split_bp2', 'mate_bp1', 'mate_bp2',
            'clip_bp1', 'clip_bp2', 'shared_reads', 'total_score',
            'n_total_reads_bp1', 'n_total_reads_bp2']
    with open(ev_file, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for row in evidence_rows:
            f.write('\t'.join(str(row.get(c, 0)) for c in cols) + '\n')

    print(f"Generated {len(evidence_rows)} fake samples ({n_per_group} per group)")
    print(f"  Groups file: {groups_file}")
    print(f"  Evidence file: {ev_file}")

    return ev_file, groups_file


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN ANALYSIS: LOAD + TEST + REPORT
# ═══════════════════════════════════════════════════════════════════════════════

def load_evidence_table(path):
    """Load per-sample evidence TSV."""
    rows = []
    with open(path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            d = dict(zip(header, parts))
            rows.append(d)
    return rows


def run_all_tests(evidence_rows, outdir, label):
    """
    Run all statistical tests and write results.
    """
    os.makedirs(outdir, exist_ok=True)

    # Build contingency data
    group_order = ['REF', 'HET', 'INV']
    counts = {g: {'yes': 0, 'no': 0, 'total': 0} for g in group_order}

    for row in evidence_rows:
        g = row['group']
        if g not in counts:
            continue
        s = row['support']
        counts[g]['total'] += 1
        if s == 'yes':
            counts[g]['yes'] += 1
        else:
            counts[g]['no'] += 1

    # Print contingency table
    print("\n" + "=" * 70)
    print("BREAKPOINT SUPPORT CONTINGENCY TABLE")
    print("=" * 70)
    print(f"{'Group':<10} {'Support_Yes':>12} {'Support_No':>12} {'Total':>8} {'Frac_Yes':>10}")
    print("-" * 55)
    for g in group_order:
        y, n, t = counts[g]['yes'], counts[g]['no'], counts[g]['total']
        frac = y / t if t > 0 else 0
        print(f"{g:<10} {y:>12} {n:>12} {t:>8} {frac:>10.3f}")
    print("-" * 55)
    total_y = sum(counts[g]['yes'] for g in group_order)
    total_n = sum(counts[g]['no'] for g in group_order)
    total_t = total_y + total_n
    print(f"{'Total':<10} {total_y:>12} {total_n:>12} {total_t:>8}")

    results = {}

    # ── TEST 1: Fisher's exact — INV vs non-INV ──────────────────────────
    print("\n" + "=" * 70)
    print("TEST 1: Fisher's exact test (INV vs non-INV)")
    print("=" * 70)

    a = counts['INV']['yes']
    b = counts['INV']['no']
    c = counts['REF']['yes'] + counts['HET']['yes']
    d = counts['REF']['no'] + counts['HET']['no']

    or_val, p_fisher = fisher_exact_2x2(a, b, c, d)
    ci_lo, ci_hi = fisher_exact_ci(a, b, c, d)

    print(f"  INV support:     {a}/{a+b} ({a/(a+b)*100:.1f}%)")
    print(f"  non-INV support: {c}/{c+d} ({c/(c+d)*100:.1f}%)")
    print(f"  Odds ratio:      {or_val:.2f} (95% CI: {ci_lo:.2f}–{ci_hi:.2f})")
    print(f"  P-value:         {p_fisher:.2e}")
    sig1 = "***" if p_fisher < 0.001 else "**" if p_fisher < 0.01 else "*" if p_fisher < 0.05 else "ns"
    print(f"  Significance:    {sig1}")

    results['fisher_inv_vs_noninv'] = {
        'odds_ratio': or_val, 'ci_lo': ci_lo, 'ci_hi': ci_hi,
        'p_value': p_fisher, 'sig': sig1,
        'inv_yes': a, 'inv_no': b, 'noninv_yes': c, 'noninv_no': d,
    }

    # ── TEST 1b: Fisher's exact — INV vs REF only ────────────────────────
    print("\n" + "-" * 50)
    print("TEST 1b: Fisher's exact test (INV vs REF)")
    a2 = counts['INV']['yes']
    b2 = counts['INV']['no']
    c2 = counts['REF']['yes']
    d2 = counts['REF']['no']
    or2, p2 = fisher_exact_2x2(a2, b2, c2, d2)
    ci2_lo, ci2_hi = fisher_exact_ci(a2, b2, c2, d2)
    print(f"  Odds ratio: {or2:.2f} (95% CI: {ci2_lo:.2f}–{ci2_hi:.2f})")
    print(f"  P-value:    {p2:.2e}")

    results['fisher_inv_vs_ref'] = {
        'odds_ratio': or2, 'ci_lo': ci2_lo, 'ci_hi': ci2_hi, 'p_value': p2,
    }

    # ── TEST 2: Chi-square 3×2 table ─────────────────────────────────────
    print("\n" + "=" * 70)
    print("TEST 2: Chi-square test of independence (REF/HET/INV × support)")
    print("=" * 70)

    table_3x2 = [
        [counts['REF']['yes'], counts['REF']['no']],
        [counts['HET']['yes'], counts['HET']['no']],
        [counts['INV']['yes'], counts['INV']['no']],
    ]
    chi2, p_chi, df = chi_square_rxc(table_3x2)
    sig2 = "***" if p_chi < 0.001 else "**" if p_chi < 0.01 else "*" if p_chi < 0.05 else "ns"
    print(f"  Chi-square:  {chi2:.3f}")
    print(f"  df:          {df}")
    print(f"  P-value:     {p_chi:.2e}")
    print(f"  Significance:{sig2}")

    results['chi_square_3x2'] = {
        'chi2': chi2, 'df': df, 'p_value': p_chi, 'sig': sig2,
    }

    # ── TEST 3: Cochran–Armitage trend ───────────────────────────────────
    print("\n" + "=" * 70)
    print("TEST 3: Cochran–Armitage trend test (REF=0, HET=1, INV=2)")
    print("=" * 70)

    z_ca, p_ca = cochran_armitage_trend(table_3x2, scores=[0, 1, 2])
    sig3 = "***" if p_ca < 0.001 else "**" if p_ca < 0.01 else "*" if p_ca < 0.05 else "ns"
    print(f"  Z-statistic: {z_ca:.3f}")
    print(f"  P-value:     {p_ca:.2e}")
    print(f"  Direction:   {'REF < HET < INV (expected)' if z_ca > 0 else 'unexpected direction'}")
    print(f"  Significance:{sig3}")

    results['cochran_armitage'] = {
        'Z': z_ca, 'p_value': p_ca, 'sig': sig3,
    }

    # ── Write results TSV ────────────────────────────────────────────────
    results_file = os.path.join(outdir, f"{label}_statistical_tests.tsv")
    with open(results_file, 'w') as f:
        f.write("test\tstatistic\tvalue\tp_value\tsignificance\n")
        f.write(f"Fisher_INV_vs_nonINV\todds_ratio\t{or_val:.4f}\t{p_fisher:.2e}\t{sig1}\n")
        f.write(f"Fisher_INV_vs_nonINV\tCI_95_lo\t{ci_lo:.4f}\t\t\n")
        f.write(f"Fisher_INV_vs_nonINV\tCI_95_hi\t{ci_hi:.4f}\t\t\n")
        f.write(f"Fisher_INV_vs_REF\todds_ratio\t{or2:.4f}\t{p2:.2e}\t\n")
        f.write(f"Chi_square_3x2\tchi2\t{chi2:.4f}\t{p_chi:.2e}\t{sig2}\n")
        f.write(f"Cochran_Armitage\tZ\t{z_ca:.4f}\t{p_ca:.2e}\t{sig3}\n")

    # ── Write text report ────────────────────────────────────────────────
    report_file = os.path.join(outdir, f"{label}_report.txt")
    with open(report_file, 'w') as f:
        f.write("=" * 72 + "\n")
        f.write("INVERSION BREAKPOINT VALIDATION — STATISTICAL REPORT\n")
        f.write("=" * 72 + "\n\n")

        f.write("CONTINGENCY TABLE:\n\n")
        f.write(f"  {'Group':<10} {'Yes':>6} {'No':>6} {'Total':>6} {'%Yes':>8}\n")
        for g in group_order:
            y, n, t = counts[g]['yes'], counts[g]['no'], counts[g]['total']
            # BUGFIX 2026-04-17 (chat 5, FIX 24): guard division by zero
            # when a group is empty (e.g. no INV samples identified yet).
            pct = f"{y/t*100:>7.1f}%" if t > 0 else f"{'–':>7}"
            f.write(f"  {g:<10} {y:>6} {n:>6} {t:>6} {pct}\n")
        f.write("\n")

        f.write("TEST RESULTS:\n\n")
        f.write(f"  1. Fisher's exact (INV vs non-INV):\n")
        f.write(f"     OR = {or_val:.2f} (95% CI: {ci_lo:.2f}–{ci_hi:.2f}), P = {p_fisher:.2e} {sig1}\n\n")
        f.write(f"  2. Chi-square (3 × 2):\n")
        f.write(f"     χ² = {chi2:.3f}, df = {df}, P = {p_chi:.2e} {sig2}\n\n")
        f.write(f"  3. Cochran–Armitage trend:\n")
        f.write(f"     Z = {z_ca:.3f}, P = {p_ca:.2e} {sig3}\n\n")

        f.write("MANUSCRIPT SENTENCE:\n\n")
        f.write(f"  Inversion-state samples showed breakpoint support in {a}/{a+b} cases\n")
        f.write(f"  versus {c}/{c+d} non-INV samples (Fisher's exact test, odds ratio =\n")
        f.write(f"  {or_val:.1f}, 95% CI = {ci_lo:.1f}–{ci_hi:.1f}, P = {p_fisher:.2e}).\n")
        f.write(f"  A significant monotonic trend of increasing support from REF to HET\n")
        f.write(f"  to INV was observed (Cochran–Armitage Z = {z_ca:.2f}, P = {p_ca:.2e}).\n")

    print(f"\n  Results written to: {results_file}")
    print(f"  Report written to:  {report_file}")

    return results, counts


# ═══════════════════════════════════════════════════════════════════════════════
# PLOTTING
# ═══════════════════════════════════════════════════════════════════════════════

def generate_plots(evidence_rows, results, counts, outdir, label,
                   bp1_chrom="C_gar_LG12", bp1_pos=52400000,
                   bp2_chrom="C_gar_LG12", bp2_pos=54600000):
    """
    Generate publication-quality breakpoint validation figure.

    Panel A: Stacked bar — support proportion by group (REF/HET/INV)
             colored blue/white/red matching the PCA and heatmap convention
    Panel B: Per-sample evidence score strip chart by group
    Panel C: Evidence breakdown heatmap (disc/split/mate/clip per sample)
    Panel D: Contingency table with Fisher + CA p-values annotated
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from matplotlib.gridspec import GridSpec
        import numpy as np
    except ImportError:
        print("WARNING: matplotlib/numpy not available. Skipping plots.", file=sys.stderr)
        return

    # ── Color scheme: REF=blue, HET=grey, INV=red (matching your figures) ──
    group_colors = {'REF': '#4169E1', 'HET': '#A0A0A0', 'INV': '#DC143C'}
    support_colors = {'yes': '#2E8B57', 'no': '#CCCCCC'}
    group_order = ['REF', 'HET', 'INV']

    fig = plt.figure(figsize=(20, 16), facecolor='white')
    gs = GridSpec(3, 3, figure=fig, hspace=0.38, wspace=0.32,
                  left=0.06, right=0.96, top=0.94, bottom=0.03,
                  height_ratios=[1, 0.9, 0.95])

    # Title
    inv_size_mb = abs(bp2_pos - bp1_pos) / 1e6
    fig.suptitle(
        f"Breakpoint Validation: {bp1_chrom} {bp1_pos/1e6:.1f}–{bp2_pos/1e6:.1f} Mb "
        f"(~{inv_size_mb:.1f} Mb inversion)",
        fontsize=15, fontweight='bold', y=0.97
    )

    # ── PANEL A: Support proportion barplot ──────────────────────────────
    ax_a = fig.add_subplot(gs[0, 0])

    bar_x = np.arange(len(group_order))
    yes_counts = [counts[g]['yes'] for g in group_order]
    no_counts = [counts[g]['no'] for g in group_order]
    totals = [counts[g]['total'] for g in group_order]
    yes_frac = [y / t if t > 0 else 0 for y, t in zip(yes_counts, totals)]

    bars_yes = ax_a.bar(bar_x, yes_counts, color=[group_colors[g] for g in group_order],
                        edgecolor='black', linewidth=0.8, label='Support')
    bars_no = ax_a.bar(bar_x, no_counts, bottom=yes_counts,
                       color=['white'] * 3, edgecolor='black', linewidth=0.8,
                       label='No support', hatch='///')

    # Annotate fractions
    for i, (g, frac) in enumerate(zip(group_order, yes_frac)):
        ax_a.text(i, totals[i] + 0.5, f"{frac*100:.0f}%",
                  ha='center', va='bottom', fontsize=11, fontweight='bold',
                  color=group_colors[g])
        ax_a.text(i, yes_counts[i] / 2, f"{yes_counts[i]}",
                  ha='center', va='center', fontsize=10, color='white', fontweight='bold')

    ax_a.set_xticks(bar_x)
    ax_a.set_xticklabels([f"{g}\n(n={t})" for g, t in zip(group_order, totals)],
                          fontsize=10)
    ax_a.set_ylabel("Number of samples", fontsize=11)
    ax_a.set_title("A   Breakpoint support by genotype", fontsize=12, fontweight='bold',
                   loc='left')
    ax_a.spines['top'].set_visible(False)
    ax_a.spines['right'].set_visible(False)
    ax_a.legend(loc='upper left', frameon=False, fontsize=9)

    # ── PANEL B: Evidence score strip chart ──────────────────────────────
    ax_b = fig.add_subplot(gs[0, 1])

    for i, g in enumerate(group_order):
        scores = [int(r['total_score']) for r in evidence_rows if r['group'] == g]
        jitter = np.random.uniform(-0.15, 0.15, len(scores))
        colors = [group_colors[g] if s > 0 else '#DDDDDD' for s in scores]
        ax_b.scatter(np.full(len(scores), i) + jitter, scores,
                     c=colors, s=40, alpha=0.7, edgecolors='black', linewidths=0.4,
                     zorder=3)
        # Median line
        if scores:
            median = sorted(scores)[len(scores) // 2]
            ax_b.hlines(median, i - 0.3, i + 0.3, color='black', linewidth=2, zorder=4)

    ax_b.set_xticks(range(len(group_order)))
    ax_b.set_xticklabels(group_order, fontsize=11)
    ax_b.set_ylabel("Evidence score", fontsize=11)
    ax_b.set_title("B   Per-sample evidence score", fontsize=12, fontweight='bold',
                   loc='left')
    ax_b.spines['top'].set_visible(False)
    ax_b.spines['right'].set_visible(False)
    ax_b.set_xlim(-0.5, 2.5)

    # ── PANEL C: Evidence heatmap ────────────────────────────────────────
    ax_c = fig.add_subplot(gs[0, 2])

    evidence_cols = ['disc_bp1', 'disc_bp2', 'split_bp1', 'split_bp2',
                     'mate_bp1', 'mate_bp2', 'clip_bp1', 'clip_bp2', 'shared_reads']
    col_labels = ['Disc\nBP1', 'Disc\nBP2', 'Split\nBP1', 'Split\nBP2',
                  'Mate\nBP1', 'Mate\nBP2', 'Clip\nBP1', 'Clip\nBP2', 'Shared\nreads']

    # Sort samples by group then by score
    sorted_rows = sorted(evidence_rows,
                         key=lambda r: (group_order.index(r['group']),
                                        -int(r['total_score'])))

    mat = np.array([[int(r[c]) for c in evidence_cols] for r in sorted_rows])
    sample_groups_sorted = [r['group'] for r in sorted_rows]

    im = ax_c.imshow(mat, aspect='auto', cmap='YlOrRd', interpolation='nearest',
                     vmin=0, vmax=max(5, mat.max()))
    ax_c.set_xticks(range(len(col_labels)))
    ax_c.set_xticklabels(col_labels, fontsize=7, rotation=0, ha='center')
    ax_c.set_ylabel("Samples (grouped)", fontsize=10)
    ax_c.set_title("C   Evidence breakdown", fontsize=12, fontweight='bold', loc='left')

    # Group boundaries
    prev_g = sample_groups_sorted[0]
    for idx, g in enumerate(sample_groups_sorted):
        if g != prev_g:
            ax_c.axhline(idx - 0.5, color='black', linewidth=1.5)
            prev_g = g

    # Group labels on right
    for g in group_order:
        indices = [i for i, gg in enumerate(sample_groups_sorted) if gg == g]
        if indices:
            mid = (indices[0] + indices[-1]) / 2
            ax_c.text(len(evidence_cols) + 0.3, mid, g, fontsize=9,
                      va='center', ha='left', color=group_colors[g], fontweight='bold')

    plt.colorbar(im, ax=ax_c, shrink=0.6, label='Count')

    # ── PANEL D: Statistical results summary ─────────────────────────────
    ax_d = fig.add_subplot(gs[1, 0])
    ax_d.axis('off')

    r = results
    fisher_r = r['fisher_inv_vs_noninv']
    chi_r = r['chi_square_3x2']
    ca_r = r['cochran_armitage']

    text_lines = [
        ("Statistical Tests", 16, 'bold', 'normal', 'black'),
        ("", 8, 'normal', 'normal', 'black'),
        ("1. Fisher's exact (INV vs non-INV)", 11, 'bold', 'normal', '#333333'),
        (f"   OR = {fisher_r['odds_ratio']:.1f}  "
         f"(95% CI: {fisher_r['ci_lo']:.1f}–{fisher_r['ci_hi']:.1f})", 10, 'normal', 'normal', '#333333'),
        (f"   P = {fisher_r['p_value']:.2e}  {fisher_r['sig']}", 11, 'bold', 'normal',
         '#DC143C' if fisher_r['p_value'] < 0.05 else '#333333'),
        ("", 6, 'normal', 'normal', 'black'),
        ("2. \u03c7\u00b2 test (3 \u00d7 2 table)", 11, 'bold', 'normal', '#333333'),
        (f"   \u03c7\u00b2 = {chi_r['chi2']:.2f}, df = {chi_r['df']}, "
         f"P = {chi_r['p_value']:.2e}  {chi_r['sig']}", 10, 'normal', 'normal', '#333333'),
        ("", 6, 'normal', 'normal', 'black'),
        ("3. Cochran\u2013Armitage trend", 11, 'bold', 'normal', '#333333'),
        (f"   Z = {ca_r['Z']:.2f}, P = {ca_r['p_value']:.2e}  {ca_r['sig']}", 10, 'normal', 'normal',
         '#333333'),
        (f"   Direction: REF \u2192 HET \u2192 INV {'(confirmed)' if ca_r['Z'] > 0 else '(unexpected)'}",
         10, 'normal', 'italic', '#2E8B57' if ca_r['Z'] > 0 else '#DC143C'),
    ]

    y = 0.95
    for text, size, weight, style, color in text_lines:
        if text:
            ax_d.text(0.05, y, text, fontsize=size, fontweight=weight, fontstyle=style,
                      color=color, transform=ax_d.transAxes, va='top', fontfamily='monospace')
        y -= 0.065

    # ── PANEL E: Proportion bar with trend ───────────────────────────────
    ax_e = fig.add_subplot(gs[1, 1])

    fracs = [counts[g]['yes'] / counts[g]['total'] * 100 if counts[g]['total'] > 0 else 0
             for g in group_order]
    bars = ax_e.bar(range(3), fracs, color=[group_colors[g] for g in group_order],
                    edgecolor='black', linewidth=0.8, width=0.6)

    # Trend line
    ax_e.plot(range(3), fracs, 'k--', linewidth=1.5, marker='D', markersize=6,
              markerfacecolor='black', zorder=5)

    for i, (g, f) in enumerate(zip(group_order, fracs)):
        ax_e.text(i, f + 2, f"{f:.0f}%", ha='center', fontsize=11, fontweight='bold')

    ax_e.set_xticks(range(3))
    ax_e.set_xticklabels([f"{g}\n(n={counts[g]['total']})" for g in group_order], fontsize=10)
    ax_e.set_ylabel("% samples with support", fontsize=11)
    ax_e.set_ylim(0, 110)
    ax_e.set_title("D   Support frequency trend", fontsize=12, fontweight='bold', loc='left')
    ax_e.spines['top'].set_visible(False)
    ax_e.spines['right'].set_visible(False)

    # CA annotation
    ax_e.text(1, 105, f"Cochran–Armitage P = {ca_r['p_value']:.2e}",
              ha='center', fontsize=9, fontstyle='italic')

    # ── PANEL F: Odds Ratio forest plot (THE EFFECT SIZE) ────────────────
    ax_f = fig.add_subplot(gs[1, 2])

    # Three comparisons: INV vs non-INV, INV vs REF, INV vs HET
    fisher_inv_ref = results.get('fisher_inv_vs_ref', {})

    # INV vs HET
    # BUGFIX 2026-04-17 (chat 5, FIX 23): original code was
    #   or_het, p_het, ci_het = fisher_exact_2x2(...), None, None
    #   if isinstance(or_het, tuple): ... else: (dead branch) ...
    # which tupled the 2-tuple return into a 3-tuple, making the else-branch
    # unreachable, and then the whole triple (or_het, p_het, ci_het) was
    # never consumed because the `comparisons` list at L898 only contains
    # INV-vs-non-INV and INV-vs-REF. The INV-vs-HET comparison was
    # half-wired. Simplify to straight-line and leave a TODO to actually
    # plot it if the manuscript asks for it.
    a_h = counts['INV']['yes']
    b_h = counts['INV']['no']
    c_h = counts['HET']['yes']
    d_h = counts['HET']['no']
    or_het, p_het = fisher_exact_2x2(a_h, b_h, c_h, d_h)
    ci_het = fisher_exact_ci(a_h, b_h, c_h, d_h)
    # TODO (manuscript): if INV-vs-HET is wanted in the forest plot,
    # append {'label': 'INV vs HET', 'or': or_het, 'ci_lo': ci_het[0],
    # 'ci_hi': ci_het[1], 'p': p_het, 'color': '#8B2252'} to `comparisons`
    # below and widen y_positions to [2, 1, 0]. Currently unused.

    comparisons = [
        {
            'label': 'INV vs non-INV',
            'or': fisher_r['odds_ratio'],
            'ci_lo': fisher_r['ci_lo'],
            'ci_hi': fisher_r['ci_hi'],
            'p': fisher_r['p_value'],
            'color': '#DC143C',
        },
        {
            'label': 'INV vs REF',
            'or': fisher_inv_ref.get('odds_ratio', 0),
            'ci_lo': fisher_inv_ref.get('ci_lo', 0),
            'ci_hi': fisher_inv_ref.get('ci_hi', 0),
            'p': fisher_inv_ref.get('p_value', 1),
            'color': '#8B0000',
        },
    ]

    y_positions = [1, 0]
    for i, comp in enumerate(comparisons):
        or_val = comp['or']
        ci_lo = comp['ci_lo']
        ci_hi = comp['ci_hi']

        if or_val <= 0 or or_val == float('inf'):
            continue

        # Clip CI for display (log scale)
        ci_hi_clip = min(ci_hi, 1000)
        ci_lo_clip = max(ci_lo, 0.01)

        # CI bar
        ax_f.plot([ci_lo_clip, ci_hi_clip], [y_positions[i], y_positions[i]],
                  color=comp['color'], linewidth=3, solid_capstyle='round')
        # Point estimate
        ax_f.plot(or_val, y_positions[i], 'D', color=comp['color'],
                  markersize=12, markeredgecolor='black', markeredgewidth=0.8, zorder=5)

        # Label
        sig = "***" if comp['p'] < 0.001 else "**" if comp['p'] < 0.01 else "*" if comp['p'] < 0.05 else "ns"
        ax_f.text(ci_hi_clip * 1.3, y_positions[i],
                  f"OR = {or_val:.1f}  [{ci_lo:.1f}–{ci_hi:.1f}]  {sig}",
                  va='center', fontsize=9, fontweight='bold', color=comp['color'])

        # Comparison label on left
        ax_f.text(0.008, y_positions[i], comp['label'],
                  va='center', ha='right', fontsize=9, color='#333333')

    # Null line at OR = 1
    ax_f.axvline(1, color='black', linestyle='--', linewidth=1.2, alpha=0.7, zorder=1)
    ax_f.text(1, -0.6, 'OR = 1\n(no association)', ha='center', fontsize=8,
              color='grey', fontstyle='italic')

    # Enrichment direction labels
    ax_f.text(0.03, -0.85, '← non-INV enriched', fontsize=7, color='grey', ha='left')
    ax_f.text(200, -0.85, 'INV enriched →', fontsize=7, color='grey', ha='right')

    ax_f.set_xscale('log')
    ax_f.set_xlim(0.02, 500)
    ax_f.set_ylim(-1.2, 2)
    ax_f.set_yticks([])
    ax_f.set_xlabel("Odds Ratio (log scale)", fontsize=10)
    ax_f.set_title("E   Association strength (Odds Ratio)", fontsize=12,
                   fontweight='bold', loc='left')
    ax_f.spines['top'].set_visible(False)
    ax_f.spines['right'].set_visible(False)
    ax_f.spines['left'].set_visible(False)

    # Shade the enrichment side
    ax_f.axvspan(1, 500, alpha=0.04, color='red')
    ax_f.axvspan(0.02, 1, alpha=0.04, color='blue')

    # ── ROW 3: Contingency table + result summary (Figure 2.0 style) ─────
    ax_g = fig.add_subplot(gs[2, :])  # span full width
    ax_g.axis('off')

    # Build the integrated result table
    all_y = sum(counts[g]['yes'] for g in group_order)
    all_n = sum(counts[g]['no'] for g in group_order)
    all_t = all_y + all_n

    cell_text = [
        ['Group', 'n', 'Support\nYES', 'Support\nNO', '% Support',
         'Fisher OR\n(vs non-INV)', '95% CI', 'P-value', 'Significance'],
    ]

    # REF row
    ry, rn = counts['REF']['yes'], counts['REF']['no']
    rt = ry + rn
    cell_text.append([
        'REF', str(rt), str(ry), str(rn),
        f"{ry/rt*100:.1f}%" if rt > 0 else "–",
        'ref', '–', '–', '–'
    ])

    # HET row
    hy, hn = counts['HET']['yes'], counts['HET']['no']
    ht = hy + hn
    cell_text.append([
        'HET', str(ht), str(hy), str(hn),
        f"{hy/ht*100:.1f}%" if ht > 0 else "–",
        '–', '–', '–', '–'
    ])

    # INV row (the key result)
    iy, ino = counts['INV']['yes'], counts['INV']['no']
    it = iy + ino
    cell_text.append([
        'INV', str(it), str(iy), str(ino),
        f"{iy/it*100:.1f}%" if it > 0 else "–",
        f"{fisher_r['odds_ratio']:.1f}",
        f"{fisher_r['ci_lo']:.1f}–{fisher_r['ci_hi']:.1f}",
        f"{fisher_r['p_value']:.2e}",
        fisher_r['sig'],
    ])

    # Total row
    cell_text.append([
        'Total', str(all_t), str(all_y), str(all_n),
        f"{all_y/all_t*100:.1f}%" if all_t > 0 else "–",
        '', '', '', ''
    ])

    # CA trend row
    cell_text.append([
        '', '', '', '', '',
        'Cochran–Armitage', f"Z = {ca_r['Z']:.2f}",
        f"{ca_r['p_value']:.2e}", ca_r['sig'],
    ])

    n_rows = len(cell_text)
    n_cols = len(cell_text[0])

    col_widths = [0.07, 0.05, 0.07, 0.07, 0.08, 0.12, 0.10, 0.10, 0.09]

    table = ax_g.table(
        cellText=cell_text, loc='center', cellLoc='center',
        colWidths=col_widths,
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 1.9)

    # Style header row
    for j in range(n_cols):
        table[0, j].set_facecolor('#2C3E50')
        table[0, j].set_text_props(color='white', fontweight='bold', fontsize=8)

    # Style group name cells
    for i, g in enumerate(['REF', 'HET', 'INV']):
        table[i+1, 0].set_facecolor(group_colors[g])
        table[i+1, 0].set_text_props(color='white', fontweight='bold')

    # Style total row
    for j in range(n_cols):
        table[4, j].set_facecolor('#F0F0F0')
        table[4, j].set_text_props(fontweight='bold')

    # Style CA row
    for j in range(n_cols):
        table[5, j].set_facecolor('#FFF8E1')

    # Highlight the INV OR cell
    table[3, 5].set_text_props(fontweight='bold', color='#DC143C')
    table[3, 7].set_text_props(fontweight='bold', color='#DC143C')
    table[3, 8].set_text_props(fontweight='bold', color='#DC143C')

    # Highlight CA result
    table[5, 5].set_text_props(fontweight='bold')
    table[5, 7].set_text_props(fontweight='bold',
                                color='#DC143C' if ca_r['p_value'] < 0.05 else '#333333')

    ax_g.set_title("F   Contingency Table & Statistical Results",
                   fontsize=13, fontweight='bold', loc='left', pad=15)

    # ── Save ─────────────────────────────────────────────────────────────
    for ext in ['pdf', 'png']:
        out_path = os.path.join(outdir, f"{label}_breakpoint_validation.{ext}")
        fig.savefig(out_path, dpi=400, bbox_inches='tight', facecolor='white')
        print(f"  Plot saved: {out_path}")

    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════════
# REAL DATA PIPELINE
# ═══════════════════════════════════════════════════════════════════════════════

def run_real_pipeline(args):
    """Run on real BAM data."""
    import pysam  # will fail early if not available

    os.makedirs(args.outdir, exist_ok=True)

    # Load sample groups
    sample_groups = {}
    with open(args.sample_groups) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                sample_groups[parts[0]] = parts[1]

    print(f"Loaded {len(sample_groups)} sample groups")
    for g in ['REF', 'HET', 'INV']:
        n = sum(1 for v in sample_groups.values() if v == g)
        print(f"  {g}: {n}")

    # Extract evidence
    evidence_rows = []
    for sid, group in sorted(sample_groups.items()):
        bam_path = os.path.join(args.bam_dir, f"{sid}.markdup.bam")
        if not os.path.exists(bam_path):
            print(f"  WARNING: {bam_path} not found, skipping", file=sys.stderr)
            continue

        print(f"  Processing {sid} ({group})...", end='', flush=True)

        ev_bp1 = extract_evidence_pysam(
            bam_path, args.bp1_chrom, args.bp1_pos, args.bp_window,
            args.min_mapq, args.min_clip_len, args.bp2_chrom, args.bp2_pos
        )
        ev_bp2 = extract_evidence_pysam(
            bam_path, args.bp2_chrom, args.bp2_pos, args.bp_window,
            args.min_mapq, args.min_clip_len, args.bp1_chrom, args.bp1_pos
        )

        support, score, details = classify_support(ev_bp1, ev_bp2)

        row = {
            'sample': sid, 'group': group,
            'support': 'yes' if support else 'no',
            'total_score': score,
            'n_total_reads_bp1': ev_bp1['n_total_reads'],
            'n_total_reads_bp2': ev_bp2['n_total_reads'],
        }
        row.update(details)
        evidence_rows.append(row)
        print(f" score={score}, support={'YES' if support else 'no'}")

    # Write evidence table
    ev_file = os.path.join(args.outdir, "per_sample_breakpoint_evidence.tsv")
    cols = ['sample', 'group', 'support', 'disc_bp1', 'disc_bp2',
            'split_bp1', 'split_bp2', 'mate_bp1', 'mate_bp2',
            'clip_bp1', 'clip_bp2', 'shared_reads', 'total_score',
            'n_total_reads_bp1', 'n_total_reads_bp2']
    with open(ev_file, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for row in evidence_rows:
            f.write('\t'.join(str(row.get(c, 0)) for c in cols) + '\n')

    return evidence_rows


# ═══════════════════════════════════════════════════════════════════════════════
# SLURM WRAPPER GENERATOR
# ═══════════════════════════════════════════════════════════════════════════════

def write_slurm_wrapper(outdir, args):
    """Write a SLURM submission script for running on the cluster."""
    script = os.path.join(outdir, "run_breakpoint_validator.sh")
    text = f"""#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=64G
#SBATCH -t 0-04:00:00
#SBATCH -J bp_valid
#SBATCH -o logs/breakpoint_validator.%j.out
#SBATCH -e logs/breakpoint_validator.%j.err
set -euo pipefail
source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="${{SLURM_SUBMIT_DIR:-$(pwd)}}"
mkdir -p logs

echo "=== Breakpoint Validator ==="
echo "BP1: {args.bp1_chrom}:{args.bp1_pos}"
echo "BP2: {args.bp2_chrom}:{args.bp2_pos}"
echo "Window: ±{args.bp_window} bp"
echo ""

python3 "${{SCRIPT_DIR}}/breakpoint_validator_standalone.py" \\
  --bam_dir "{args.bam_dir or '/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/delly_sv/00_markdup'}" \\
  --sample_groups "${{SCRIPT_DIR}}/sample_groups.tsv" \\
  --bp1_chrom {args.bp1_chrom} --bp1_pos {args.bp1_pos} \\
  --bp2_chrom {args.bp2_chrom} --bp2_pos {args.bp2_pos} \\
  --bp_window {args.bp_window} --min_mapq {args.min_mapq} \\
  --outdir "${{SCRIPT_DIR}}/results" \\
  --label "{args.label}"

echo "=== DONE ==="
"""
    with open(script, 'w') as f:
        f.write(text)
    os.chmod(script, 0o755)
    print(f"  SLURM wrapper: {script}")


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    if args.generate_test_data:
        # ── Test mode: generate fake data and run full pipeline ──────────
        print("=" * 70)
        print("GENERATING TEST DATA AND RUNNING FULL PIPELINE")
        print("=" * 70)

        ev_file, groups_file = generate_fake_bam_data(
            args.outdir, n_per_group=args.n_per_group,
            bp1_chrom=args.bp1_chrom, bp1_pos=args.bp1_pos,
            bp2_chrom=args.bp2_chrom, bp2_pos=args.bp2_pos,
        )

        evidence_rows = load_evidence_table(ev_file)
        results, counts = run_all_tests(evidence_rows, args.outdir, args.label)
        generate_plots(evidence_rows, results, counts, args.outdir, args.label,
                       args.bp1_chrom, args.bp1_pos, args.bp2_chrom, args.bp2_pos)

    else:
        # ── Real data mode ───────────────────────────────────────────────
        if not args.bam_dir or not args.sample_groups:
            print("ERROR: --bam_dir and --sample_groups required for real data mode")
            print("Use --generate_test_data for testing")
            sys.exit(1)

        print("=" * 70)
        print("BREAKPOINT VALIDATION — REAL DATA")
        print("=" * 70)
        print(f"  BP1: {args.bp1_chrom}:{args.bp1_pos}")
        print(f"  BP2: {args.bp2_chrom}:{args.bp2_pos}")
        print(f"  Window: ±{args.bp_window} bp")
        print(f"  MAPQ: ≥{args.min_mapq}")

        evidence_rows = run_real_pipeline(args)
        results, counts = run_all_tests(evidence_rows, args.outdir, args.label)
        generate_plots(evidence_rows, results, counts, args.outdir, args.label,
                       args.bp1_chrom, args.bp1_pos, args.bp2_chrom, args.bp2_pos)

        write_slurm_wrapper(args.outdir, args)

    print("\n" + "=" * 70)
    print("COMPLETE")
    print("=" * 70)
    print(f"  Output directory: {args.outdir}")


if __name__ == "__main__":
    main()
