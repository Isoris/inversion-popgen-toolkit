#!/usr/bin/env python3
"""
population_regenotype.py — Bayesian population-prior SV regenotyping
=====================================================================
Rescues per-sample genotype dropout caused by stochastic coverage at
low-coverage WGS (~9×). Reads raw FORMAT fields from merged DELLY2 or
Manta VCFs and applies a cohort frequency prior to recompute genotype
posteriors for low-evidence samples.

PROBLEM:
  At 9× mean coverage, ~34% of true HET inversion carriers have ≤3
  supporting reads at a breakpoint. DELLY/Manta genotype these as 0/0
  independently per sample, with no population prior. This propagates
  false-REF calls into all downstream analyses (Fisher OR, decomposition
  seeds, carrier concordance for BND rescue).

SOLUTION:
  1. For each SV site, identify "confident" samples (high evidence) and
     estimate cohort allele frequency from them.
  2. For "low-evidence" samples, recompute P(GT | data, freq) using:
       P(GT=0/1 | data) ∝ P(data | GT=0/1) × P(GT=0/1 | freq)
     where the likelihood uses a binomial model on alt/(alt+ref) reads.
  3. Output a corrected GT matrix with both hard calls and posterior
     probabilities.

CALLER SUPPORT:
  DELLY2:  FORMAT/DR,DV (spanning ref/alt), FORMAT/RR,RV (junction ref/alt)
  Manta:   FORMAT/PR (ref,alt), FORMAT/SR (ref,alt)

USAGE:
  python3 population_regenotype.py \\
    --vcf cohort_226.INV.raw.vcf.gz \\
    --caller delly \\
    --outdir regenotyped/ \\
    [--sv_type INV|BND|ALL] \\
    [--min_confident_reads 4] \\
    [--min_confident_gq 10] \\
    [--prior_floor 0.005] \\
    [--posterior_threshold 0.5] \\
    [--error_rate 0.01]

OUTPUT:
  regenotyped/
    population_regenotype.GT_matrix.tsv     — corrected GT matrix
    population_regenotype.posteriors.tsv.gz  — full posterior probabilities
    population_regenotype.rescue_summary.tsv — per-site rescue statistics
    population_regenotype.log               — detailed log

PIPELINE INTEGRATION:
  Slots in as STEP_A03b after bcftools merge (STEP 5 in MODULE_4D/4E),
  before strict filtering (STEP 8). The corrected GT matrix replaces
  the raw GT matrix for all downstream use.

Author: Claude (for Quentin, Kasetsart University)
Date: 2026-04-16
"""
import argparse
import gzip
import math
import os
import sys
from collections import defaultdict, Counter

# =============================================================================
# CONSTANTS
# =============================================================================
DELLY_FORMAT_FIELDS = {
    'spanning_ref': 'DR',  # spanning ref reads
    'spanning_alt': 'DV',  # spanning alt reads (discordant pairs)
    'junction_ref': 'RR',  # junction ref reads
    'junction_alt': 'RV',  # junction alt reads (split reads)
}
MANTA_FORMAT_FIELDS = {
    'paired_ref_alt': 'PR',  # ref,alt paired reads
    'split_ref_alt':  'SR',  # ref,alt split reads
}

# =============================================================================
# MATH HELPERS
# =============================================================================

def log_binom_pmf(k, n, p):
    """Log probability of binomial(k; n, p). Handles edge cases."""
    if p <= 0:
        return 0.0 if k == 0 else float('-inf')
    if p >= 1:
        return 0.0 if k == n else float('-inf')
    if n == 0:
        return 0.0 if k == 0 else float('-inf')
    # Use log-gamma for binomial coefficient
    try:
        log_coeff = (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1))
        return log_coeff + k * math.log(p) + (n - k) * math.log(1 - p)
    except (ValueError, OverflowError):
        return float('-inf')


def genotype_posteriors(n_alt, n_ref, freq, error_rate=0.01):
    """
    Compute posterior P(GT | data, freq) for GT in {0/0, 0/1, 1/1}.

    Model:
      GT=0/0: expect alt reads ~ Binom(n_total, error_rate)
      GT=0/1: expect alt reads ~ Binom(n_total, 0.5)
      GT=1/1: expect alt reads ~ Binom(n_total, 1 - error_rate)

    Prior from population frequency (Hardy-Weinberg as approximation,
    though we know F1 hybrids violate HWE — the prior is still
    directionally correct and conservative):
      P(0/0) = (1-freq)^2
      P(0/1) = 2*freq*(1-freq)
      P(1/1) = freq^2

    For the F1 hybrid context where HWE doesn't hold, the prior
    essentially captures "what fraction of the cohort carries this
    variant" — even if the genotype frequencies don't follow HWE,
    the prior still helps because a common variant is more likely
    to be present in any given sample.
    """
    n_total = n_alt + n_ref

    if n_total == 0:
        # No data — return prior only
        p00 = (1 - freq) ** 2
        p01 = 2 * freq * (1 - freq)
        p11 = freq ** 2
        total = p00 + p01 + p11
        return p00/total, p01/total, p11/total

    # Log-likelihoods
    ll_00 = log_binom_pmf(n_alt, n_total, error_rate)
    ll_01 = log_binom_pmf(n_alt, n_total, 0.5)
    ll_11 = log_binom_pmf(n_alt, n_total, 1.0 - error_rate)

    # Log-priors (with floor to avoid log(0))
    freq_safe = max(freq, 1e-10)
    freq_safe = min(freq_safe, 1 - 1e-10)
    lp_00 = 2 * math.log(1 - freq_safe)
    lp_01 = math.log(2) + math.log(freq_safe) + math.log(1 - freq_safe)
    lp_11 = 2 * math.log(freq_safe)

    # Log-posteriors (unnormalized)
    log_post = [ll_00 + lp_00, ll_01 + lp_01, ll_11 + lp_11]

    # Normalize in log-space
    max_lp = max(log_post)
    if max_lp == float('-inf'):
        return 1/3, 1/3, 1/3

    exp_posts = [math.exp(lp - max_lp) for lp in log_post]
    total = sum(exp_posts)
    if total == 0:
        return 1/3, 1/3, 1/3

    return exp_posts[0]/total, exp_posts[1]/total, exp_posts[2]/total


# =============================================================================
# VCF PARSING
# =============================================================================

def parse_vcf_header(header_lines):
    """Extract sample names from VCF header."""
    for line in header_lines:
        if line.startswith('#CHROM'):
            fields = line.strip().split('\t')
            return fields[9:]  # sample names
    return []


def extract_delly_evidence(format_str, sample_str):
    """
    Extract alt/ref read counts from DELLY FORMAT fields.
    Returns (n_alt, n_ref, is_confident_carrier, original_gt).

    DELLY FORMAT: GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV
    We want DR (spanning ref), DV (spanning alt), RR (junction ref), RV (junction alt).
    Total alt = DV + RV, total ref = DR + RR.
    """
    fmt_keys = format_str.split(':')
    vals = sample_str.split(':')

    # Build dict
    d = {}
    for i, key in enumerate(fmt_keys):
        d[key] = vals[i] if i < len(vals) else '.'

    gt = d.get('GT', './.')

    # Extract read counts
    def safe_int(s):
        try:
            return int(s)
        except (ValueError, TypeError):
            return 0

    dr = safe_int(d.get('DR', '0'))
    dv = safe_int(d.get('DV', '0'))
    rr = safe_int(d.get('RR', '0'))
    rv = safe_int(d.get('RV', '0'))

    n_alt = dv + rv
    n_ref = dr + rr

    # GQ for confidence
    gq = safe_int(d.get('GQ', '0'))

    return n_alt, n_ref, gq, gt


def extract_manta_evidence(format_str, sample_str):
    """
    Extract alt/ref read counts from Manta FORMAT fields.
    Returns (n_alt, n_ref, gq, original_gt).

    Manta FORMAT: GT:FT:GQ:PL:PR:SR
    PR = ref,alt (paired reads)
    SR = ref,alt (split reads)
    """
    fmt_keys = format_str.split(':')
    vals = sample_str.split(':')

    d = {}
    for i, key in enumerate(fmt_keys):
        d[key] = vals[i] if i < len(vals) else '.'

    gt = d.get('GT', './.')

    def parse_pair(s):
        """Parse 'ref,alt' format."""
        try:
            parts = s.split(',')
            return int(parts[0]), int(parts[1])
        except (ValueError, IndexError, TypeError):
            return 0, 0

    pr_ref, pr_alt = parse_pair(d.get('PR', '0,0'))
    sr_ref, sr_alt = parse_pair(d.get('SR', '0,0'))

    n_alt = pr_alt + sr_alt
    n_ref = pr_ref + sr_ref

    gq = 0
    try:
        gq = int(d.get('GQ', '0'))
    except (ValueError, TypeError):
        pass

    return n_alt, n_ref, gq, gt


# =============================================================================
# MAIN REGENOTYPING
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description="Population-prior SV regenotyping for low-coverage WGS",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    p.add_argument("--vcf", required=True,
                   help="Input merged VCF/BCF (raw, before strict filtering)")
    p.add_argument("--caller", required=True, choices=["delly", "manta"],
                   help="SV caller (determines FORMAT field parsing)")
    p.add_argument("--outdir", required=True,
                   help="Output directory")
    p.add_argument("--sv_type", default="ALL",
                   help="Filter to this SV type (INV, BND, DEL, ALL) [ALL]")

    # Evidence thresholds
    p.add_argument("--min_confident_reads", type=int, default=4,
                   help="Min alt reads for a sample to be 'confident carrier' [4]")
    p.add_argument("--min_confident_gq", type=int, default=10,
                   help="Min GQ for a sample to contribute to frequency estimate [10]")

    # Prior parameters
    p.add_argument("--prior_floor", type=float, default=0.005,
                   help="Minimum allele frequency prior (avoids zero prior) [0.005]")
    p.add_argument("--posterior_threshold", type=float, default=0.5,
                   help="P(carrier) threshold for hard call rescue [0.5]")
    p.add_argument("--error_rate", type=float, default=0.01,
                   help="Expected alt read rate in true REF samples [0.01]")

    # Output control
    p.add_argument("--write_posteriors", action="store_true", default=True,
                   help="Write full posterior probabilities (large file)")

    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    log_path = os.path.join(args.outdir, "population_regenotype.log")
    logf = open(log_path, 'w')

    def log(msg):
        line = f"[REGENO] {msg}"
        print(line)
        logf.write(line + '\n')

    log(f"Input VCF: {args.vcf}")
    log(f"Caller: {args.caller}")
    log(f"SV type filter: {args.sv_type}")
    log(f"Min confident reads: {args.min_confident_reads}")
    log(f"Prior floor: {args.prior_floor}")
    log(f"Posterior threshold: {args.posterior_threshold}")
    log(f"Error rate: {args.error_rate}")

    # Choose evidence extractor
    if args.caller == "delly":
        extract_fn = extract_delly_evidence
    else:
        extract_fn = extract_manta_evidence

    # ── Pass 1: Read VCF, compute per-site frequency, regenotype ──────────

    opener = gzip.open if args.vcf.endswith('.gz') else open
    header_lines = []
    sample_names = []
    sites = []

    # Collect all data first
    log("Pass 1: Reading VCF and extracting per-sample evidence...")

    n_records = 0
    n_skipped_type = 0

    with opener(args.vcf, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                header_lines.append(line)
                if line.startswith('#CHROM'):
                    sample_names = line.strip().split('\t')[9:]
                    log(f"  Samples: {len(sample_names)}")
                continue

            n_records += 1
            fields = line.strip().split('\t')
            chrom, pos, sv_id = fields[0], int(fields[1]), fields[2]
            qual = fields[5]
            filt = fields[6]
            info_str = fields[7]
            format_str = fields[8]

            # Parse INFO
            info = {}
            for kv in info_str.split(';'):
                if '=' in kv:
                    k, v = kv.split('=', 1)
                    info[k] = v
                else:
                    info[kv] = True

            svtype = info.get('SVTYPE', '.')
            end = info.get('END', pos)
            ct = info.get('CT', '.')

            # Filter by SV type
            if args.sv_type != "ALL" and svtype != args.sv_type:
                n_skipped_type += 1
                continue

            # Extract per-sample evidence
            sample_data = []
            for si, sample_str in enumerate(fields[9:]):
                n_alt, n_ref, gq, orig_gt = extract_fn(format_str, sample_str)
                sample_data.append({
                    'n_alt': n_alt,
                    'n_ref': n_ref,
                    'n_total': n_alt + n_ref,
                    'gq': gq,
                    'orig_gt': orig_gt,
                    'raw_sample_str': sample_str,
                })

            # ── Estimate cohort frequency from confident calls ────────
            # "Confident carrier" = alt reads >= threshold AND original GT != 0/0
            # "Confident non-carrier" = alt reads == 0 AND total reads >= threshold
            n_confident_carrier = 0
            n_confident_ref = 0
            n_ambiguous = 0

            for sd in sample_data:
                gt = sd['orig_gt']
                if gt in ('./.', '.'):
                    n_ambiguous += 1
                    continue

                if sd['n_alt'] >= args.min_confident_reads and sd['gq'] >= args.min_confident_gq:
                    if gt in ('0/1', '0|1', '1|0', '1/1', '1|1'):
                        n_confident_carrier += 1
                    elif gt == '0/0' or gt == '0|0':
                        # Has alt reads but called REF — suspicious, count as carrier
                        # (DELLY's GL model may have been wrong)
                        n_confident_carrier += 1
                elif sd['n_alt'] == 0 and sd['n_total'] >= args.min_confident_reads:
                    n_confident_ref += 1
                else:
                    n_ambiguous += 1

            n_confident_total = n_confident_carrier + n_confident_ref
            if n_confident_total > 0:
                # Carrier frequency (not allele frequency — we don't know ploidy state)
                # For diploid: carrier_freq ≈ 2*allele_freq for rare variants
                carrier_freq = n_confident_carrier / n_confident_total
                # Convert to allele frequency (approximate)
                # For HWE: carrier_freq = 2pq + q² = 1 - (1-q)² → q = 1 - sqrt(1-cf)
                # But simpler: if most carriers are HET, allele_freq ≈ carrier_freq / 2
                # We use carrier_freq directly as our prior parameter since the
                # posterior model parameterizes by allele freq
                allele_freq = max(args.prior_floor,
                                  min(0.5, 1.0 - math.sqrt(max(0, 1.0 - carrier_freq))))
            else:
                allele_freq = args.prior_floor

            # ── Regenotype each sample ────────────────────────────────
            rescued_count = 0
            new_gts = []

            for sd in sample_data:
                p00, p01, p11 = genotype_posteriors(
                    sd['n_alt'], sd['n_ref'], allele_freq, args.error_rate
                )
                p_carrier = p01 + p11

                # Decision logic
                orig_is_ref = sd['orig_gt'] in ('0/0', '0|0')
                orig_is_carrier = sd['orig_gt'] in ('0/1', '0|1', '1|0', '1/1', '1|1')
                orig_is_missing = sd['orig_gt'] in ('./.', '.')

                if orig_is_ref and p_carrier >= args.posterior_threshold:
                    # RESCUE: was called REF but posterior says carrier
                    if p11 > p01:
                        new_gt = '1/1'
                    else:
                        new_gt = '0/1'
                    rescued_count += 1
                    rescue_flag = 'RESCUED'
                elif orig_is_missing and p_carrier >= args.posterior_threshold:
                    if p11 > p01:
                        new_gt = '1/1'
                    else:
                        new_gt = '0/1'
                    rescued_count += 1
                    rescue_flag = 'RESCUED_FROM_MISSING'
                elif orig_is_carrier and p_carrier < 0.1:
                    # DEMOTE: was called carrier but posterior says very unlikely
                    # (rare — only happens if a single noise read in a very rare variant)
                    new_gt = '0/0'
                    rescue_flag = 'DEMOTED'
                else:
                    new_gt = sd['orig_gt']
                    rescue_flag = 'KEPT'

                sd['new_gt'] = new_gt
                sd['p00'] = p00
                sd['p01'] = p01
                sd['p11'] = p11
                sd['p_carrier'] = p_carrier
                sd['rescue_flag'] = rescue_flag
                new_gts.append(new_gt)

            # Count new genotypes
            new_carrier_count = sum(1 for g in new_gts if g in ('0/1', '0|1', '1|0', '1/1', '1|1'))
            orig_carrier_count = sum(1 for sd in sample_data
                                     if sd['orig_gt'] in ('0/1', '0|1', '1|0', '1/1', '1|1'))

            sites.append({
                'chrom': chrom,
                'pos': pos,
                'sv_id': sv_id,
                'svtype': svtype,
                'end': end,
                'ct': ct,
                'qual': qual,
                'filter': filt,
                'allele_freq_est': allele_freq,
                'n_confident_carrier': n_confident_carrier,
                'n_confident_ref': n_confident_ref,
                'n_ambiguous': n_ambiguous,
                'orig_carrier_count': orig_carrier_count,
                'new_carrier_count': new_carrier_count,
                'rescued_count': rescued_count,
                'sample_data': sample_data,
                'new_gts': new_gts,
            })

            if n_records % 500 == 0:
                log(f"  Processed {n_records} records...")

    log(f"  Total records: {n_records}")
    log(f"  Skipped (type filter): {n_skipped_type}")
    log(f"  Processed: {len(sites)}")

    # ── Write outputs ─────────────────────────────────────────────────────

    # 1. Corrected GT matrix
    gt_path = os.path.join(args.outdir, "population_regenotype.GT_matrix.tsv")
    log(f"Writing corrected GT matrix: {gt_path}")
    with open(gt_path, 'w') as out:
        header = ['CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'CT',
                  'FREQ_EST', 'ORIG_CARRIERS', 'NEW_CARRIERS', 'RESCUED']
        header.extend(sample_names)
        out.write('\t'.join(header) + '\n')

        for site in sites:
            row = [
                site['chrom'], str(site['pos']), str(site['end']),
                site['sv_id'], site['svtype'], site['ct'],
                f"{site['allele_freq_est']:.4f}",
                str(site['orig_carrier_count']),
                str(site['new_carrier_count']),
                str(site['rescued_count']),
            ]
            row.extend(site['new_gts'])
            out.write('\t'.join(row) + '\n')

    # 2. Full posteriors (gzipped — can be large)
    if args.write_posteriors:
        post_path = os.path.join(args.outdir, "population_regenotype.posteriors.tsv.gz")
        log(f"Writing posteriors: {post_path}")
        with gzip.open(post_path, 'wt') as out:
            header = ['sv_id', 'sample', 'orig_gt', 'new_gt', 'rescue_flag',
                      'n_alt', 'n_ref', 'n_total', 'gq',
                      'p00', 'p01', 'p11', 'p_carrier', 'allele_freq_est']
            out.write('\t'.join(header) + '\n')

            for site in sites:
                for si, sd in enumerate(site['sample_data']):
                    row = [
                        site['sv_id'],
                        sample_names[si],
                        sd['orig_gt'],
                        sd['new_gt'],
                        sd['rescue_flag'],
                        str(sd['n_alt']),
                        str(sd['n_ref']),
                        str(sd['n_total']),
                        str(sd['gq']),
                        f"{sd['p00']:.4f}",
                        f"{sd['p01']:.4f}",
                        f"{sd['p11']:.4f}",
                        f"{sd['p_carrier']:.4f}",
                        f"{site['allele_freq_est']:.4f}",
                    ]
                    out.write('\t'.join(row) + '\n')

    # 3. Rescue summary
    summary_path = os.path.join(args.outdir, "population_regenotype.rescue_summary.tsv")
    log(f"Writing rescue summary: {summary_path}")
    with open(summary_path, 'w') as out:
        header = ['sv_id', 'chrom', 'pos', 'end', 'svtype', 'ct', 'qual', 'filter',
                  'allele_freq_est', 'n_confident_carrier', 'n_confident_ref',
                  'n_ambiguous', 'orig_carriers', 'new_carriers', 'rescued',
                  'pct_increase']
        out.write('\t'.join(header) + '\n')

        for site in sites:
            pct = ((site['new_carrier_count'] - site['orig_carrier_count'])
                   / max(1, site['orig_carrier_count']) * 100)
            row = [
                site['sv_id'], site['chrom'], str(site['pos']),
                str(site['end']), site['svtype'], site['ct'],
                site['qual'], site['filter'],
                f"{site['allele_freq_est']:.4f}",
                str(site['n_confident_carrier']),
                str(site['n_confident_ref']),
                str(site['n_ambiguous']),
                str(site['orig_carrier_count']),
                str(site['new_carrier_count']),
                str(site['rescued_count']),
                f"{pct:.1f}",
            ]
            out.write('\t'.join(row) + '\n')

    # ── Summary statistics ────────────────────────────────────────────────

    total_rescued = sum(s['rescued_count'] for s in sites)
    total_orig_carriers = sum(s['orig_carrier_count'] for s in sites)
    total_new_carriers = sum(s['new_carrier_count'] for s in sites)
    sites_with_rescue = sum(1 for s in sites if s['rescued_count'] > 0)

    log("")
    log("=" * 60)
    log("REGENOTYPING SUMMARY")
    log("=" * 60)
    log(f"  Sites processed:          {len(sites)}")
    log(f"  Sites with ≥1 rescue:     {sites_with_rescue}")
    log(f"  Total original carriers:  {total_orig_carriers}")
    log(f"  Total new carriers:       {total_new_carriers}")
    log(f"  Total rescued samples:    {total_rescued}")
    if total_orig_carriers > 0:
        log(f"  Overall increase:         {(total_new_carriers - total_orig_carriers) / total_orig_carriers * 100:.1f}%")

    # Per-site rescue distribution
    rescue_counts = Counter(s['rescued_count'] for s in sites)
    log("")
    log("  Rescue distribution (rescued_per_site: n_sites):")
    for k in sorted(rescue_counts.keys()):
        if k <= 5 or rescue_counts[k] > 1:
            log(f"    {k}: {rescue_counts[k]}")
    if max(rescue_counts.keys(), default=0) > 5:
        log(f"    >5: {sum(v for k, v in rescue_counts.items() if k > 5)}")

    # Frequency distribution of allele freq estimates
    freq_bins = [0, 0.01, 0.05, 0.10, 0.20, 0.50, 1.01]
    freq_labels = ['<1%', '1-5%', '5-10%', '10-20%', '20-50%', '>50%']
    freq_counts = [0] * len(freq_labels)
    for s in sites:
        for bi in range(len(freq_bins) - 1):
            if freq_bins[bi] <= s['allele_freq_est'] < freq_bins[bi + 1]:
                freq_counts[bi] += 1
                break
    log("")
    log("  Allele frequency distribution:")
    for label, count in zip(freq_labels, freq_counts):
        log(f"    {label}: {count}")

    log("")
    log(f"Output: {args.outdir}/")
    log("=" * 60)

    logf.close()


if __name__ == "__main__":
    main()
