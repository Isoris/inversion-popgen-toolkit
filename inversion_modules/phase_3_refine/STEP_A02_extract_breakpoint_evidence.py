#!/usr/bin/env python3
"""
STEP_A02_extract_breakpoint_evidence.py
=======================================
For each DELLY2 INV candidate, extract per-sample breakpoint evidence
from markdup BAMs and assign sample groups from two sources:

  Source A: Snake 2 band assignments (majority vote across overlapping windows)
  Source B: Local PCA on the INV region (direct k=3 on PC1 of dosage)

Outputs per candidate:
  {inv_id}_evidence.tsv        — per-sample evidence + group + support
  {inv_id}_group_assignments.tsv — sample groups from both sources

Usage:
  python3 STEP_A02_extract_breakpoint_evidence.py \\
    --candidates matched_inv_candidates.tsv \\
    --bam_dir /path/to/00_markdup \\
    --samples samples_all_226.txt \\
    [--snake2_bands snake2_band_assignments.tsv.gz] \\
    [--dosage_dir /path/to/04_dosage_by_chr] \\
    --outdir 01_per_sample_evidence
"""
import argparse
import gzip
import math
import os
import sys
from collections import Counter, defaultdict

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--candidates", required=True,
                   help="matched_inv_candidates.tsv from STEP_A01")
    p.add_argument("--bam_dir", required=True,
                   help="Directory with {sample}.markdup.bam files")
    p.add_argument("--samples", required=True,
                   help="samples_all_226.txt")
    p.add_argument("--snake2_bands", default="",
                   help="snake2_band_assignments.tsv.gz (optional)")
    p.add_argument("--dosage_dir", default="",
                   help="Dosage directory for local PCA fallback (optional)")
    p.add_argument("--outdir", required=True)
    p.add_argument("--bp_window", type=int, default=300)
    p.add_argument("--min_mapq", type=int, default=20)
    p.add_argument("--min_clip_len", type=int, default=10)
    p.add_argument("--threads", type=int, default=8)
    return p.parse_args()


def extract_evidence_pysam(bam_path, chrom, pos, window, min_mapq, min_clip_len,
                           partner_chrom, partner_pos):
    """Extract breakpoint evidence from BAM. Returns dict."""
    import pysam

    ev = {
        'n_discordant_FF': 0, 'n_discordant_RR': 0,
        'n_split_reads': 0, 'n_soft_clip': 0,
        'n_mate_at_partner': 0, 'n_total_reads': 0,
        'supporting_read_names': set(),
    }

    start = max(0, pos - window)
    end = pos + window

    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        return ev

    try:
        for read in bam.fetch(chrom, start, end):
            if read.is_unmapped or read.is_secondary or read.is_duplicate:
                continue
            if read.mapping_quality < min_mapq:
                continue
            ev['n_total_reads'] += 1

            if read.is_paired and not read.mate_is_unmapped:
                read_fwd = not read.is_reverse
                mate_fwd = not read.mate_is_reverse
                if read_fwd == mate_fwd:  # same orientation = inversion signal
                    key = 'n_discordant_FF' if read_fwd else 'n_discordant_RR'
                    ev[key] += 1
                    ev['supporting_read_names'].add(read.query_name)

                if (read.next_reference_name == partner_chrom and
                    abs(read.next_reference_start - partner_pos) <= window * 3):
                    ev['n_mate_at_partner'] += 1
                    ev['supporting_read_names'].add(read.query_name)

            if read.has_tag('SA'):
                for sa_entry in read.get_tag('SA').split(';'):
                    if not sa_entry.strip():
                        continue
                    parts = sa_entry.strip().split(',')
                    if len(parts) >= 5:
                        sa_chr, sa_pos, sa_mapq = parts[0], int(parts[1]), int(parts[4])
                        if (sa_mapq >= min_mapq and sa_chr == partner_chrom and
                            abs(sa_pos - partner_pos) <= window * 3):
                            ev['n_split_reads'] += 1
                            ev['supporting_read_names'].add(read.query_name)

            cigar = read.cigartuples
            if cigar:
                if cigar[0][0] == 4 and cigar[0][1] >= min_clip_len:
                    ev['n_soft_clip'] += 1
                    ev['supporting_read_names'].add(read.query_name)
                if cigar[-1][0] == 4 and cigar[-1][1] >= min_clip_len:
                    ev['n_soft_clip'] += 1
    except:
        pass

    bam.close()
    return ev


def assign_groups_snake2(snake2_file, chrom, bp1, bp2, samples):
    """
    Assign REF/HET/INV from Snake 2 band assignments.
    Majority vote across windows overlapping [bp1, bp2].
    """
    groups = {}
    if not snake2_file or not os.path.exists(snake2_file):
        return groups

    opener = gzip.open if snake2_file.endswith('.gz') else open
    band_counts = defaultdict(lambda: Counter())  # sample -> Counter of bands

    with opener(snake2_file, 'rt') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            p = line.strip().split('\t')
            d = dict(zip(header, p))
            wc = d.get('chrom', '')
            if wc != chrom:
                continue
            ws = int(d.get('start_bp', d.get('window_start', 0)))
            we = int(d.get('end_bp', d.get('window_end', 0)))
            # Window overlaps candidate?
            if we < bp1 or ws > bp2:
                continue
            # Parse sample assignments
            for col in header:
                if col.startswith('band_') or col.startswith('sample_band_'):
                    continue
            # Look for band assignment columns
            band_col = None
            for c in ['band_assignments', 'sample_bands', 'band_indices']:
                if c in d:
                    band_col = c
                    break

            if band_col and d[band_col]:
                # Format: sample1:0,sample2:1,sample3:2 or similar
                for pair in d[band_col].split(','):
                    if ':' in pair:
                        sid, band = pair.split(':', 1)
                        band_counts[sid][band] += 1

    # Majority vote
    for sid in samples:
        if sid in band_counts:
            most_common = band_counts[sid].most_common(1)[0][0]
            # Map band index to REF/HET/INV
            # Convention: band 0/2 = tails (REF/INV), band 1 = middle (HET)
            # But which tail is which depends on PC1 polarity
            if most_common in ('0', 'tailA', 'REF'):
                groups[sid] = 'REF'
            elif most_common in ('1', 'middle', 'HET'):
                groups[sid] = 'HET'
            elif most_common in ('2', 'tailB', 'INV'):
                groups[sid] = 'INV'
            else:
                groups[sid] = 'UNKNOWN'

    return groups


def assign_groups_local_pca(dosage_dir, chrom, bp1, bp2, samples):
    """
    Fallback: do local PCA on dosage in the INV region, k-means k=3 on PC1.
    Returns {sample: REF/HET/INV}.
    """
    import numpy as np

    groups = {}
    if not dosage_dir:
        return groups

    # Find dosage file for this chromosome
    dosage_file = None
    for ext in ['.dosage.tsv.gz', '.dosage.tsv', '.beagle.gz']:
        candidate = os.path.join(dosage_dir, f"{chrom}{ext}")
        if os.path.exists(candidate):
            dosage_file = candidate
            break

    if not dosage_file:
        print(f"  WARNING: No dosage file for {chrom} in {dosage_dir}")
        return groups

    # Load dosage for region
    opener = gzip.open if dosage_file.endswith('.gz') else open
    sample_indices = []
    dosage_matrix = []

    with opener(dosage_file, 'rt') as f:
        header = f.readline().strip().split('\t')
        # Find sample columns (typically after marker/pos/ref/alt columns)
        meta_cols = 0
        for i, h in enumerate(header):
            if h in samples:
                meta_cols = i
                break
        sample_names = header[meta_cols:]

        for line in f:
            p = line.strip().split('\t')
            pos = int(p[1]) if len(p) > 1 else 0
            if pos < bp1 or pos > bp2:
                continue
            vals = []
            for v in p[meta_cols:]:
                try:
                    vals.append(float(v))
                except:
                    vals.append(0.0)
            dosage_matrix.append(vals)

    if len(dosage_matrix) < 10:
        return groups

    mat = np.array(dosage_matrix).T  # samples × markers
    # Remove zero-variance markers
    var = mat.var(axis=0)
    mat = mat[:, var > 0]
    if mat.shape[1] < 5:
        return groups

    # PCA
    mat_c = mat - mat.mean(axis=0)
    U, S, Vt = np.linalg.svd(mat_c, full_matrices=False)
    pc1 = U[:, 0] * S[0]

    # k-means k=3 on PC1
    from numpy import array
    centers = np.array([np.percentile(pc1, 15), np.median(pc1), np.percentile(pc1, 85)])
    for _ in range(20):
        labels = np.argmin(np.abs(pc1[:, None] - centers[None, :]), axis=1)
        for k in range(3):
            mask = labels == k
            if mask.sum() > 0:
                centers[k] = pc1[mask].mean()

    # Sort centers to assign REF (lowest) / HET (middle) / INV (highest)
    order = np.argsort(centers)
    label_map = {order[0]: 'REF', order[1]: 'HET', order[2]: 'INV'}

    for i, sid in enumerate(sample_names):
        if i < len(labels):
            groups[sid] = label_map.get(labels[i], 'UNKNOWN')

    return groups


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    import pysam  # early import check

    # Load samples
    with open(args.samples) as f:
        samples = [l.strip() for l in f if l.strip()]
    print(f"Samples: {len(samples)}")

    # Load candidates
    candidates = []
    with open(args.candidates) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            p = line.strip().split('\t')
            d = dict(zip(header, p))
            candidates.append(d)
    print(f"Candidates: {len(candidates)}")

    # Process each candidate
    summary_rows = []

    for ci, cand in enumerate(candidates):
        inv_id = cand['inv_id']
        chrom = cand['chrom']
        bp1 = int(cand['bp1_pos'])
        bp2 = int(cand['bp2_pos'])
        svlen = int(cand['svlen_bp'])

        # Parse CIPOS/CIEND to widen breakpoint search window
        # CIPOS format: "-152,153" meaning bp1 is uncertain in [bp1-152, bp1+153]
        # CIEND format: "-171,171" meaning bp2 is uncertain in [bp2-171, bp2+171]
        bp1_window = args.bp_window  # default fallback
        bp2_window = args.bp_window
        cipos_str = cand.get('cipos', '.')
        ciend_str = cand.get('ciend', '.')
        if cipos_str and cipos_str != '.':
            try:
                ci_lo, ci_hi = cipos_str.split(',')
                bp1_window = max(args.bp_window, abs(int(ci_lo)) + abs(int(ci_hi)) + 50)
            except:
                pass
        if ciend_str and ciend_str != '.':
            try:
                ci_lo, ci_hi = ciend_str.split(',')
                bp2_window = max(args.bp_window, abs(int(ci_lo)) + abs(int(ci_hi)) + 50)
            except:
                pass

        print(f"\n[{ci+1}/{len(candidates)}] {inv_id}: {chrom}:{bp1}-{bp2} ({svlen/1e6:.2f} Mb)"
              f" bp1_window={bp1_window} bp2_window={bp2_window}"
              f" CIPOS={cipos_str} CIEND={ciend_str}")

        # ── Assign groups ────────────────────────────────────────────────
        groups_snake2 = assign_groups_snake2(
            args.snake2_bands, chrom, bp1, bp2, samples
        )
        groups_pca = assign_groups_local_pca(
            args.dosage_dir, chrom, bp1, bp2, samples
        )

        # Use Snake2 if available, else PCA, else skip
        if len(groups_snake2) >= len(samples) * 0.5:
            groups = groups_snake2
            group_source = "snake2"
        elif len(groups_pca) >= len(samples) * 0.5:
            groups = groups_pca
            group_source = "local_pca"
        else:
            print(f"  WARNING: Could not assign groups for {inv_id}, skipping")
            continue

        n_ref = sum(1 for v in groups.values() if v == 'REF')
        n_het = sum(1 for v in groups.values() if v == 'HET')
        n_inv = sum(1 for v in groups.values() if v == 'INV')
        print(f"  Groups ({group_source}): REF={n_ref}, HET={n_het}, INV={n_inv}")

        # ── Extract evidence ─────────────────────────────────────────────
        evidence_rows = []
        for si, sid in enumerate(samples):
            bam_path = os.path.join(args.bam_dir, f"{sid}.markdup.bam")
            if not os.path.exists(bam_path):
                continue

            ev_bp1 = extract_evidence_pysam(
                bam_path, chrom, bp1, bp1_window, args.min_mapq,
                args.min_clip_len, chrom, bp2
            )
            ev_bp2 = extract_evidence_pysam(
                bam_path, chrom, bp2, bp2_window, args.min_mapq,
                args.min_clip_len, chrom, bp1
            )

            disc = (ev_bp1['n_discordant_FF'] + ev_bp1['n_discordant_RR'] +
                    ev_bp2['n_discordant_FF'] + ev_bp2['n_discordant_RR'])
            split = ev_bp1['n_split_reads'] + ev_bp2['n_split_reads']
            mate = ev_bp1['n_mate_at_partner'] + ev_bp2['n_mate_at_partner']
            clip = ev_bp1['n_soft_clip'] + ev_bp2['n_soft_clip']
            shared = len(ev_bp1['supporting_read_names'] & ev_bp2['supporting_read_names'])
            score = disc + split * 2 + mate + shared * 3

            support = disc >= 1 or split >= 1 or mate >= 1 or shared >= 1
            group = groups.get(sid, 'UNASSIGNED')

            evidence_rows.append({
                'sample': sid,
                'group': group,
                'group_source': group_source,
                'support': 'yes' if support else 'no',
                'disc_bp1': ev_bp1['n_discordant_FF'] + ev_bp1['n_discordant_RR'],
                'disc_bp2': ev_bp2['n_discordant_FF'] + ev_bp2['n_discordant_RR'],
                'split_bp1': ev_bp1['n_split_reads'],
                'split_bp2': ev_bp2['n_split_reads'],
                'mate_bp1': ev_bp1['n_mate_at_partner'],
                'mate_bp2': ev_bp2['n_mate_at_partner'],
                'clip_bp1': ev_bp1['n_soft_clip'],
                'clip_bp2': ev_bp2['n_soft_clip'],
                'shared_reads': shared,
                'total_score': score,
                'reads_bp1': ev_bp1['n_total_reads'],
                'reads_bp2': ev_bp2['n_total_reads'],
            })

            if (si + 1) % 50 == 0:
                print(f"  ... {si+1}/{len(samples)} samples processed")

        # ── Write evidence table ─────────────────────────────────────────
        ev_file = os.path.join(args.outdir, f"{inv_id}_evidence.tsv")
        cols = list(evidence_rows[0].keys()) if evidence_rows else []
        with open(ev_file, 'w') as f:
            f.write('\t'.join(cols) + '\n')
            for row in evidence_rows:
                f.write('\t'.join(str(row[c]) for c in cols) + '\n')

        # ── Write group assignments (both sources) ───────────────────────
        grp_file = os.path.join(args.outdir, f"{inv_id}_group_assignments.tsv")
        with open(grp_file, 'w') as f:
            f.write("sample\tgroup_snake2\tgroup_pca\tgroup_used\tgroup_source\n")
            for sid in samples:
                gs2 = groups_snake2.get(sid, '.')
                gpc = groups_pca.get(sid, '.')
                gu = groups.get(sid, '.')
                f.write(f"{sid}\t{gs2}\t{gpc}\t{gu}\t{group_source}\n")

        # ── Summary ──────────────────────────────────────────────────────
        n_support = sum(1 for r in evidence_rows if r['support'] == 'yes')
        n_inv_support = sum(1 for r in evidence_rows
                          if r['group'] == 'INV' and r['support'] == 'yes')
        summary_rows.append({
            'inv_id': inv_id, 'chrom': chrom,
            'bp1': bp1, 'bp2': bp2, 'svlen': svlen,
            'n_samples': len(evidence_rows),
            'n_ref': n_ref, 'n_het': n_het, 'n_inv': n_inv,
            'n_support_total': n_support,
            'n_inv_with_support': n_inv_support,
            'frac_inv_support': n_inv_support / n_inv if n_inv > 0 else 0,
            'group_source': group_source,
        })
        print(f"  Support: {n_support}/{len(evidence_rows)} total, "
              f"{n_inv_support}/{n_inv} INV")

    # ── Write master summary ─────────────────────────────────────────────
    summary_file = os.path.join(args.outdir, "validation_candidates_summary.tsv")
    if summary_rows:
        cols = list(summary_rows[0].keys())
        with open(summary_file, 'w') as f:
            f.write('\t'.join(cols) + '\n')
            for row in summary_rows:
                f.write('\t'.join(str(row[c]) for c in cols) + '\n')
    print(f"\nSummary: {summary_file}")
    print(f"Evidence files: {args.outdir}/{{inv_id}}_evidence.tsv")


if __name__ == "__main__":
    main()
