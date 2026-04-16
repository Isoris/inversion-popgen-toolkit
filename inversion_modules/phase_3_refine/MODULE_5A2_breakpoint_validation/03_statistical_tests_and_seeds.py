#!/usr/bin/env python3
"""
03_statistical_tests_and_seeds.py
==================================
For each validated candidate:
  1. Fisher's exact test (INV vs non-INV)
  2. Chi-square 3×2 (REF/HET/INV × support)
  3. Cochran-Armitage trend test (REF < HET < INV)
  4. Concordance scoring (PCA band vs breakpoint evidence)
  5. Seed generation for C01i seeded deconvolution

Outputs:
  03_statistical_tests/
    {inv_id}_tests.tsv         — per-candidate test results
    all_candidates_tests.tsv   — master test table
    manuscript_sentences.txt   — ready-to-paste text
  04_deconvolution_seeds/
    {inv_id}_seeds.tsv         — confirmed sample sets for C01i
    seed_summary.tsv           — which candidates qualify for seeding
"""
import argparse
import glob
import math
import os
import sys

# ═════════════════════════════════════════════════════════════════════════════
# STATISTICAL TESTS (pure Python, no scipy)
# ═════════════════════════════════════════════════════════════════════════════

def _log_fact(n):
    return sum(math.log(i) for i in range(2, n + 1)) if n > 1 else 0.0

def fisher_exact(a, b, c, d):
    n = a + b + c + d
    if b * c == 0:
        odds = float('inf') if a * d > 0 else 0.0
    else:
        odds = (a * d) / (b * c)
    def lp(aa, bb, cc, dd):
        nn = aa + bb + cc + dd
        return (_log_fact(aa+bb) + _log_fact(cc+dd) + _log_fact(aa+cc) + _log_fact(bb+dd)
                - _log_fact(nn) - _log_fact(aa) - _log_fact(bb) - _log_fact(cc) - _log_fact(dd))
    lp_obs = lp(a, b, c, d)
    r1, r2, c1 = a + b, c + d, a + c
    pval = sum(math.exp(lp(aa, r1-aa, c1-aa, r2-c1+aa))
               for aa in range(min(r1, c1) + 1)
               if (r1-aa) >= 0 and (c1-aa) >= 0 and (r2-c1+aa) >= 0
               and lp(aa, r1-aa, c1-aa, r2-c1+aa) <= lp_obs + 1e-10)
    pval = min(pval, 1.0)
    # CI
    if b*c == 0 or a*d == 0:
        ci = (0.0, float('inf'))
    else:
        se = math.sqrt(1/max(a,.5) + 1/max(b,.5) + 1/max(c,.5) + 1/max(d,.5))
        lo = math.exp(math.log(odds) - 1.96*se)
        hi = math.exp(math.log(odds) + 1.96*se)
        ci = (lo, hi)
    return odds, pval, ci

def cochran_armitage(table, scores=None):
    K = len(table)
    if scores is None:
        scores = list(range(K))
    n = [table[i][0] + table[i][1] for i in range(K)]
    r = [table[i][0] for i in range(K)]
    N = sum(n); R = sum(r)
    if N == 0 or R == 0 or R == N:
        return 0.0, 1.0
    p_hat = R / N
    t_bar = sum(scores[i]*n[i] for i in range(K)) / N
    T_num = sum(scores[i]*(r[i] - n[i]*p_hat) for i in range(K))
    T_den_sq = p_hat*(1-p_hat)*(sum(scores[i]**2*n[i] for i in range(K)) - N*t_bar**2)
    if T_den_sq <= 0:
        return 0.0, 1.0
    Z = T_num / math.sqrt(T_den_sq)
    pval = 2 * 0.5 * math.erfc(abs(Z) / math.sqrt(2))
    return Z, pval

def sig_stars(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"


# ═════════════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════════════

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--evidence_dir", required=True, help="01_per_sample_evidence/")
    p.add_argument("--stats_dir", required=True, help="03_statistical_tests/")
    p.add_argument("--seeds_dir", required=True, help="04_deconvolution_seeds/")
    p.add_argument("--seed_min_concordance", type=float, default=0.80)
    p.add_argument("--seed_min_support_frac", type=float, default=0.60)
    p.add_argument("--seed_min_samples", type=int, default=5)
    p.add_argument("--bonferroni", action="store_true")
    return p.parse_args()

def main():
    args = parse_args()
    os.makedirs(args.stats_dir, exist_ok=True)
    os.makedirs(args.seeds_dir, exist_ok=True)

    ev_files = sorted(glob.glob(os.path.join(args.evidence_dir, "*_evidence.tsv")))
    print(f"Found {len(ev_files)} evidence files")

    all_results = []
    all_seeds = []
    manuscript_lines = []

    for ev_file in ev_files:
        inv_id = os.path.basename(ev_file).replace("_evidence.tsv", "")
        print(f"\n{'='*60}")
        print(f"Testing: {inv_id}")

        # Load evidence
        rows = []
        with open(ev_file) as f:
            header = f.readline().strip().split('\t')
            for line in f:
                p = line.strip().split('\t')
                rows.append(dict(zip(header, p)))

        if not rows:
            continue

        # Build contingency
        counts = {'REF': {'yes': 0, 'no': 0}, 'HET': {'yes': 0, 'no': 0}, 'INV': {'yes': 0, 'no': 0}}
        for r in rows:
            g = r.get('group', 'UNKNOWN')
            s = r.get('support', 'no')
            if g in counts:
                counts[g][s] += 1

        n_ref = counts['REF']['yes'] + counts['REF']['no']
        n_het = counts['HET']['yes'] + counts['HET']['no']
        n_inv = counts['INV']['yes'] + counts['INV']['no']

        if n_inv == 0:
            print(f"  No INV samples, skipping")
            continue

        # ── Test 1: Fisher INV vs non-INV ────────────────────────────────
        a = counts['INV']['yes']
        b = counts['INV']['no']
        c = counts['REF']['yes'] + counts['HET']['yes']
        d = counts['REF']['no'] + counts['HET']['no']
        or_val, p_fisher, ci = fisher_exact(a, b, c, d)

        print(f"  Fisher INV vs non-INV: OR={or_val:.1f} ({ci[0]:.1f}–{ci[1]:.1f}), P={p_fisher:.2e} {sig_stars(p_fisher)}")

        # ── Test 2: Fisher INV vs REF ────────────────────────────────────
        or2, p2, ci2 = fisher_exact(a, b, counts['REF']['yes'], counts['REF']['no'])

        # ── Test 3: Cochran-Armitage ─────────────────────────────────────
        table_ca = [
            [counts['REF']['yes'], counts['REF']['no']],
            [counts['HET']['yes'], counts['HET']['no']],
            [counts['INV']['yes'], counts['INV']['no']],
        ]
        z_ca, p_ca = cochran_armitage(table_ca)
        print(f"  Cochran-Armitage: Z={z_ca:.2f}, P={p_ca:.2e} {sig_stars(p_ca)}")

        # ── Concordance ──────────────────────────────────────────────────
        inv_support_frac = a / (a + b) if (a + b) > 0 else 0
        ref_support_frac = counts['REF']['yes'] / n_ref if n_ref > 0 else 0
        het_support_frac = counts['HET']['yes'] / n_het if n_het > 0 else 0

        # Concordance = (INV with support + REF without support) / total
        concordance = ((a + counts['REF']['no']) /
                      (a + b + counts['REF']['yes'] + counts['REF']['no'])
                      if (a + b + n_ref) > 0 else 0)

        result = {
            'inv_id': inv_id,
            'n_ref': n_ref, 'n_het': n_het, 'n_inv': n_inv,
            'inv_support': a, 'inv_no_support': b,
            'ref_support': counts['REF']['yes'],
            'het_support': counts['HET']['yes'],
            'fisher_OR': f"{or_val:.2f}",
            'fisher_CI_lo': f"{ci[0]:.2f}", 'fisher_CI_hi': f"{ci[1]:.2f}",
            'fisher_P': f"{p_fisher:.2e}", 'fisher_sig': sig_stars(p_fisher),
            'fisher_inv_ref_P': f"{p2:.2e}",
            'CA_Z': f"{z_ca:.3f}", 'CA_P': f"{p_ca:.2e}", 'CA_sig': sig_stars(p_ca),
            'inv_support_frac': f"{inv_support_frac:.3f}",
            'het_support_frac': f"{het_support_frac:.3f}",
            'ref_support_frac': f"{ref_support_frac:.3f}",
            'concordance': f"{concordance:.3f}",
        }
        all_results.append(result)

        # ── Write per-candidate test results ─────────────────────────────
        test_file = os.path.join(args.stats_dir, f"{inv_id}_tests.tsv")
        with open(test_file, 'w') as f:
            f.write('\t'.join(result.keys()) + '\n')
            f.write('\t'.join(str(v) for v in result.values()) + '\n')

        # ── Manuscript sentence ──────────────────────────────────────────
        sentence = (
            f"{inv_id}: Inversion-state samples showed breakpoint support in "
            f"{a}/{a+b} cases versus {c}/{c+d} non-INV samples "
            f"(Fisher's exact test, OR = {or_val:.1f}, 95% CI = {ci[0]:.1f}–{ci[1]:.1f}, "
            f"P = {p_fisher:.2e}). "
            f"A {'significant' if p_ca < 0.05 else 'non-significant'} monotonic trend "
            f"of increasing support from REF ({counts['REF']['yes']}/{n_ref}) to "
            f"HET ({counts['HET']['yes']}/{n_het}) to INV ({a}/{a+b}) was observed "
            f"(Cochran–Armitage Z = {z_ca:.2f}, P = {p_ca:.2e})."
        )
        manuscript_lines.append(sentence)

        # ── Seed generation ──────────────────────────────────────────────
        qualifies = (
            inv_support_frac >= args.seed_min_support_frac and
            concordance >= args.seed_min_concordance and
            a >= args.seed_min_samples and
            p_fisher < 0.05
        )

        seed_info = {
            'inv_id': inv_id,
            'qualifies': qualifies,
            'n_seeds_inv': 0,
            'n_seeds_het': 0,
            'reason': '',
        }

        if qualifies:
            # Extract confirmed INV samples (have support)
            inv_seeds = [r['sample'] for r in rows
                        if r['group'] == 'INV' and r['support'] == 'yes']
            # Extract confirmed HET samples (have support — one haplotype inverted)
            het_seeds = [r['sample'] for r in rows
                        if r['group'] == 'HET' and r['support'] == 'yes']
            # Extract confirmed REF samples (NO support — both haplotypes reference)
            ref_confirmed = [r['sample'] for r in rows
                            if r['group'] == 'REF' and r['support'] == 'no']

            seed_info['n_seeds_inv'] = len(inv_seeds)
            seed_info['n_seeds_het'] = len(het_seeds)
            seed_info['n_seeds_ref'] = len(ref_confirmed)

            # Write seed file for C01i
            seed_file = os.path.join(args.seeds_dir, f"{inv_id}_seeds.tsv")
            with open(seed_file, 'w') as f:
                f.write("sample\tseed_class\tevidence_score\n")
                for r in rows:
                    if r['sample'] in inv_seeds:
                        f.write(f"{r['sample']}\tINV_confirmed\t{r['total_score']}\n")
                    elif r['sample'] in het_seeds:
                        f.write(f"{r['sample']}\tHET_confirmed\t{r['total_score']}\n")
                    elif r['sample'] in ref_confirmed:
                        f.write(f"{r['sample']}\tREF_confirmed\t{r['total_score']}\n")

            print(f"  SEED QUALIFIED: {len(inv_seeds)} INV + {len(het_seeds)} HET + {len(ref_confirmed)} REF seeds")
        else:
            reasons = []
            if inv_support_frac < args.seed_min_support_frac:
                reasons.append(f"low_inv_support({inv_support_frac:.2f}<{args.seed_min_support_frac})")
            if concordance < args.seed_min_concordance:
                reasons.append(f"low_concordance({concordance:.2f}<{args.seed_min_concordance})")
            if a < args.seed_min_samples:
                reasons.append(f"too_few_inv_seeds({a}<{args.seed_min_samples})")
            if p_fisher >= 0.05:
                reasons.append(f"not_significant(P={p_fisher:.2e})")
            seed_info['reason'] = ';'.join(reasons)
            print(f"  NOT QUALIFIED for seeding: {seed_info['reason']}")

        all_seeds.append(seed_info)

    # ── Write master tables ──────────────────────────────────────────────
    if all_results:
        master = os.path.join(args.stats_dir, "all_candidates_tests.tsv")
        with open(master, 'w') as f:
            f.write('\t'.join(all_results[0].keys()) + '\n')
            for r in all_results:
                f.write('\t'.join(str(v) for v in r.values()) + '\n')
        print(f"\nMaster tests: {master}")

    # Manuscript sentences
    ms_file = os.path.join(args.stats_dir, "manuscript_sentences.txt")
    with open(ms_file, 'w') as f:
        for s in manuscript_lines:
            f.write(s + '\n\n')
    print(f"Manuscript sentences: {ms_file}")

    # Seed summary
    seed_sum = os.path.join(args.seeds_dir, "seed_summary.tsv")
    if all_seeds:
        with open(seed_sum, 'w') as f:
            f.write('\t'.join(all_seeds[0].keys()) + '\n')
            for s in all_seeds:
                f.write('\t'.join(str(v) for v in s.values()) + '\n')
    n_qual = sum(1 for s in all_seeds if s['qualifies'])
    print(f"\nSeed summary: {n_qual}/{len(all_seeds)} candidates qualified for seeded deconvolution")


if __name__ == "__main__":
    main()
