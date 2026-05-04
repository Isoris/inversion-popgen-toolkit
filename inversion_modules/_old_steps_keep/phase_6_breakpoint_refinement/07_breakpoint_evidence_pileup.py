#!/usr/bin/env python3
# =============================================================================
# 07_breakpoint_evidence_pileup.py  (v1.0)
#
# PANEL D REPLICA — per-sample stacked read evidence at both breakpoints
#
# What this plots:
#
#   For each breakpoint (left + right), renders a multi-track panel:
#     - Track 1: split reads (soft clips)      → red tick marks per sample
#     - Track 2: discordant read pairs         → colored arcs per sample
#     - Track 3: coverage                       → grey filled area per sample
#     - Track 4: assembly alignment (optional) → blue bars per sample
#
#   Samples are stacked vertically, labeled on the left margin:
#     sample1
#     sample2
#     sample3
#     ...
#
# Input options (priority order):
#   A) Pre-extracted per-sample evidence TSV (fastest, if C01f already ran)
#   B) Extract fresh from BAMs via pysam (slower, runs MODULE_5A2 logic inline)
#
# Usage (mode A — pre-extracted):
#   python 07_breakpoint_evidence_pileup.py \
#     --evidence-tsv  results/c01f/candidate_LG28_evidence.tsv.gz \
#     --sample-list   results/c01f/candidate_LG28_samples.tsv \
#     --bp-left       15115243 \
#     --bp-right      18005891 \
#     --chrom         C_gar_LG28 \
#     --out           figures/LG28_panel_D.pdf
#
# Usage (mode B — extract from BAMs):
#   python 07_breakpoint_evidence_pileup.py \
#     --bam-dir       /data/bams \
#     --sample-list   results/c01f/candidate_LG28_samples.tsv \
#     --bp-left       15115243 \
#     --bp-right      18005891 \
#     --chrom         C_gar_LG28 \
#     --window        2000 \
#     --out           figures/LG28_panel_D.pdf
#
# The --sample-list TSV must have columns:
#   sample         (sample ID, used for BAM filename lookup)
#   group          (REF / HET / INV — for sort order and optional coloring)
#
# Samples are sorted in the output by group (INV on top, HET middle, REF
# bottom), then by supporting-read count within each group (most supportive
# first).
#
# =============================================================================

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Arc, Rectangle
from matplotlib.collections import LineCollection, PatchCollection

try:
    import pysam
    HAS_PYSAM = True
except ImportError:
    HAS_PYSAM = False


# -----------------------------------------------------------------------------
# Evidence extraction (mode B — port of MODULE_5A2 / C01f logic)
# -----------------------------------------------------------------------------
def extract_evidence_pysam(bam_path, chrom, pos, window, min_mapq,
                            min_clip_len, partner_chrom, partner_pos):
    """Port of the C01f extraction logic. Returns per-read records, not
    just counts — we need individual read positions for the plot."""

    records = {
        'split_reads': [],           # list of (read_start, clip_pos, strand, sa_pos)
        'discordant_pairs': [],      # list of (read_pos, mate_pos, orient)
        'soft_clips': [],            # list of (clip_pos, clip_len, side)
        'mate_at_partner': [],       # list of (read_pos, mate_pos)
        'coverage_bins': {},         # dict: position -> depth
    }

    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        return records

    start = max(0, pos - window)
    end = pos + window

    try:
        for read in bam.fetch(chrom, start, end):
            if read.is_unmapped or read.is_secondary or read.is_duplicate:
                continue
            if read.mapping_quality < min_mapq:
                continue

            # Discordant same-orientation pairs
            if read.is_paired and not read.mate_is_unmapped:
                read_fwd = not read.is_reverse
                mate_fwd = not read.mate_is_reverse
                if read_fwd == mate_fwd and read.next_reference_name == chrom:
                    orient = 'FF' if read_fwd else 'RR'
                    records['discordant_pairs'].append((
                        read.reference_start,
                        read.next_reference_start,
                        orient
                    ))

                # Mate at the partner breakpoint
                if (read.next_reference_name == partner_chrom and
                        abs(read.next_reference_start - partner_pos) <= window * 3):
                    records['mate_at_partner'].append((
                        read.reference_start,
                        read.next_reference_start
                    ))

            # Split reads (SA tag)
            if read.has_tag('SA'):
                for sa_entry in read.get_tag('SA').split(';'):
                    if not sa_entry.strip():
                        continue
                    parts = sa_entry.strip().split(',')
                    if len(parts) >= 5:
                        sa_chr, sa_pos = parts[0], int(parts[1])
                        sa_mapq = int(parts[4])
                        if (sa_mapq >= min_mapq and sa_chr == partner_chrom and
                                abs(sa_pos - partner_pos) <= window * 3):
                            records['split_reads'].append((
                                read.reference_start,
                                read.reference_end or read.reference_start,
                                '-' if read.is_reverse else '+',
                                sa_pos
                            ))

            # Soft clips
            cigar = read.cigartuples
            if cigar:
                if cigar[0][0] == 4 and cigar[0][1] >= min_clip_len:
                    records['soft_clips'].append((
                        read.reference_start, cigar[0][1], 'left'
                    ))
                if cigar[-1][0] == 4 and cigar[-1][1] >= min_clip_len:
                    records['soft_clips'].append((
                        read.reference_end or read.reference_start,
                        cigar[-1][1], 'right'
                    ))

        # Coverage in 50bp bins
        bin_size = 50
        bin_start = (start // bin_size) * bin_size
        for col in bam.pileup(chrom, start, end, truncate=True, min_mapping_quality=min_mapq):
            bin_pos = (col.reference_pos // bin_size) * bin_size
            records['coverage_bins'][bin_pos] = \
                records['coverage_bins'].get(bin_pos, 0) + col.nsegments
    except Exception:
        pass

    bam.close()
    return records


# -----------------------------------------------------------------------------
# Rendering — Panel D replica
# -----------------------------------------------------------------------------
def render_panel(ax_reads, ax_cov, ax_asm, per_sample_records,
                  bp_pos, window, sample_labels, title_str,
                  group_map=None, group_colors=None, draw_asm=True):
    """
    Render the 3 sub-tracks for one breakpoint.

    ax_reads  axes for combined split-reads + discordant arcs (one row per sample)
    ax_cov    axes for aggregated coverage (mean + IQR across samples)
    ax_asm    axes for assembly alignment (bottom)
    per_sample_records  dict[sample] -> records dict
    bp_pos    breakpoint position
    window    half-width in bp
    sample_labels  list of sample IDs in display order (top to bottom)
    """
    n_samples = len(sample_labels)
    x_lo, x_hi = bp_pos - window, bp_pos + window
    group_map = group_map or {}
    group_colors = group_colors or {}

    # ---- Stacked reads: one row per sample ----
    for row_idx, sid in enumerate(sample_labels):
        y = n_samples - row_idx - 1  # top = first sample
        recs = per_sample_records.get(sid, {})
        grp = group_map.get(sid, 'unknown')
        row_bg_color = group_colors.get(grp, '#eee')

        # Subtle row-background stripe so you can trace the sample
        ax_reads.axhspan(y, y + 1, color=row_bg_color, alpha=0.05,
                          linewidth=0)

        # Split-read primary alignments as red line segments
        for rs, re, strand, sa_pos in recs.get('split_reads', []):
            rs_c, re_c = max(rs, x_lo), min(re, x_hi)
            if rs_c < re_c:
                ax_reads.plot([rs_c, re_c], [y + 0.78, y + 0.78],
                               color='#c62828', linewidth=1.5,
                               solid_capstyle='butt', alpha=0.85)

        # Soft-clip tick marks (red)
        for clip_pos, clip_len, side in recs.get('soft_clips', []):
            if x_lo <= clip_pos <= x_hi:
                end_x = (clip_pos + clip_len) if side == 'right' else (clip_pos - clip_len)
                ax_reads.plot([clip_pos, end_x], [y + 0.78, y + 0.78],
                               color='#c62828', linewidth=2.2,
                               solid_capstyle='butt')

        # Discordant pair arcs — FF/RR in different colors
        for pos_a, pos_b, orient in recs.get('discordant_pairs', []):
            if not (x_lo <= pos_a <= x_hi or x_lo <= pos_b <= x_hi):
                continue
            center = (pos_a + pos_b) / 2
            width = abs(pos_b - pos_a)
            if width < 50:
                continue
            color = '#1f5d8e' if orient == 'FF' else '#b03a2e'
            arc = Arc((center, y + 0.15), width, 0.55, angle=0,
                       theta1=0, theta2=180,
                       color=color, linewidth=0.6, alpha=0.7)
            ax_reads.add_patch(arc)

        # Mate-at-partner arrows (small directional tick, subtle)
        for pos_a, pos_b in recs.get('mate_at_partner', []):
            if not (x_lo <= pos_a <= x_hi):
                continue
            direction = 1 if pos_b > pos_a else -1
            ax_reads.annotate(
                '', xy=(pos_a + direction * window * 0.05, y + 0.45),
                xytext=(pos_a, y + 0.45),
                arrowprops=dict(arrowstyle='->', color='#555',
                                  lw=0.35, alpha=0.5)
            )

    ax_reads.set_xlim(x_lo, x_hi)
    ax_reads.set_ylim(-0.5, n_samples + 0.5)
    ax_reads.axvline(bp_pos, color='#d32f2f', linestyle='--', linewidth=1.2, alpha=0.85)
    ax_reads.set_yticks([])
    ax_reads.set_xticks([])
    ax_reads.spines['top'].set_visible(False)
    ax_reads.spines['right'].set_visible(False)
    ax_reads.spines['bottom'].set_visible(False)
    ax_reads.set_title(title_str, fontsize=10, color='#d32f2f', loc='left')

    # ---- Coverage: mean + 10-90% ribbon across samples ----
    all_positions = set()
    for sid in sample_labels:
        for p in per_sample_records.get(sid, {}).get('coverage_bins', {}).keys():
            if x_lo <= p <= x_hi:
                all_positions.add(p)

    if all_positions:
        positions = sorted(all_positions)
        depths = np.zeros((n_samples, len(positions)))
        for row_idx, sid in enumerate(sample_labels):
            cov = per_sample_records.get(sid, {}).get('coverage_bins', {})
            for pi, p in enumerate(positions):
                depths[row_idx, pi] = cov.get(p, 0)
        mean_cov = depths.mean(axis=0)
        p10 = np.percentile(depths, 10, axis=0)
        p90 = np.percentile(depths, 90, axis=0)

        ax_cov.fill_between(positions, p10, p90, color='#bbbbbb', alpha=0.4,
                              linewidth=0)
        ax_cov.plot(positions, mean_cov, color='#333', linewidth=0.9)
        ax_cov.set_ylim(0, max(20, p90.max() * 1.1))

    ax_cov.set_xlim(x_lo, x_hi)
    ax_cov.axvline(bp_pos, color='#d32f2f', linestyle='--', linewidth=1.2, alpha=0.85)
    ax_cov.set_ylabel('Coverage', fontsize=7, rotation=90, labelpad=3)
    ax_cov.tick_params(axis='y', labelsize=6)
    ax_cov.set_xticks([])
    ax_cov.spines['top'].set_visible(False)
    ax_cov.spines['right'].set_visible(False)

    # ---- Assembly alignment ----
    if draw_asm:
        # Without real assembly BED, draw a schematic broken bar showing that
        # contigs align cleanly on both sides but break at the bp
        n_show = min(n_samples, 15)
        for row_idx in range(n_show):
            y = n_show - row_idx - 1
            for xl, xr in [(x_lo, bp_pos - 20), (bp_pos + 20, x_hi)]:
                ax_asm.add_patch(Rectangle(
                    (xl, y + 0.15), xr - xl, 0.7,
                    facecolor='#2c5aa0',
                    edgecolor='none', alpha=0.5
                ))
        ax_asm.set_ylim(-0.5, n_show + 0.5)

    ax_asm.set_xlim(x_lo, x_hi)
    ax_asm.axvline(bp_pos, color='#d32f2f', linestyle='--', linewidth=1.2, alpha=0.85)
    ax_asm.set_yticks([])
    ax_asm.set_ylabel('Assembly', fontsize=7, rotation=90, labelpad=3)
    ax_asm.tick_params(axis='x', labelsize=6)
    ax_asm.set_xlabel(f'Position (bp) — breakpoint at {bp_pos:,}', fontsize=7)
    ax_asm.spines['top'].set_visible(False)
    ax_asm.spines['right'].set_visible(False)


# -----------------------------------------------------------------------------
# Sample-label rendering on the left of the whole figure
# -----------------------------------------------------------------------------
def draw_sample_labels(ax, sample_labels, group_map, group_colors):
    """Small leftmost axis with sample names color-coded by group."""
    n = len(sample_labels)
    for i, sid in enumerate(sample_labels):
        y = n - i - 1
        grp = group_map.get(sid, 'unknown')
        color = group_colors.get(grp, '#333')
        ax.text(0.95, y + 0.5, sid, fontsize=5.5, ha='right', va='center',
                  color=color, family='monospace')
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.5, n + 0.5)
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ('top', 'right', 'bottom', 'left'):
        ax.spines[side].set_visible(False)


# -----------------------------------------------------------------------------
# Sample ordering — INV first, then HET, then REF; within group by support
# -----------------------------------------------------------------------------
def order_samples(sample_list_df, per_sample_records):
    """Sort by: group (INV/HET/REF), then by supporting-read count desc."""
    def support_score(sid):
        recs = per_sample_records.get(sid, {})
        return (len(recs.get('split_reads', [])) +
                  len(recs.get('discordant_pairs', [])) +
                  len(recs.get('mate_at_partner', [])))

    group_order = {'INV': 0, 'HET': 1, 'REF': 2, 'unknown': 3}
    sample_list_df = sample_list_df.copy()
    sample_list_df['grp_rank'] = sample_list_df['group'].map(
        lambda g: group_order.get(g, 3)
    )
    sample_list_df['support'] = sample_list_df['sample'].map(support_score)
    sample_list_df = sample_list_df.sort_values(
        ['grp_rank', 'support'], ascending=[True, False]
    )
    return sample_list_df['sample'].tolist()


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def _find_registry_loader_py():
    """Locate registries/api/python/registry_loader.py; return import dir or None."""
    candidates = []
    base = os.environ.get('BASE', '')
    if base:
        candidates += [
            os.path.join(base, 'registries', 'api', 'python'),
            os.path.join(base, 'inversion-popgen-toolkit', 'registries', 'api', 'python'),
        ]
    candidates += [
        'registries/api/python',
        '../registries/api/python',
    ]
    for d in candidates:
        if os.path.isfile(os.path.join(d, 'registry_loader.py')):
            return os.path.abspath(d)
    return None


def main():
    p = argparse.ArgumentParser(
        description='Per-sample breakpoint evidence pileup (Panel D replica)'
    )
    p.add_argument('--evidence-tsv', default=None,
                    help='Mode A: pre-extracted per-sample evidence TSV')
    p.add_argument('--bam-dir', default=None,
                    help='Mode B: directory containing <sample>.markdup.bam')
    p.add_argument('--bam-suffix', default='.markdup.bam',
                    help='Suffix appended to sample ID to form BAM filename')
    # chat-18: registry resolution. Given a candidate ID, looks up chrom,
    # bp coords (from boundary_refined if present, else interval_registry),
    # and the sample list (from reg.samples.get_groups_for_candidate). All
    # other flags become optional overrides.
    p.add_argument('--candidate', default=None,
                    help='Candidate ID (e.g. LG28_1). Auto-resolves chrom, '
                         'bp-left, bp-right, and sample-list via the registry.')
    p.add_argument('--sample-list', default=None,
                    help='TSV with columns: sample, group (optional if '
                         '--candidate is used).')
    p.add_argument('--chrom',        default=None)
    p.add_argument('--bp-left',      type=int, default=None)
    p.add_argument('--bp-right',     type=int, default=None)
    p.add_argument('--window',       type=int, default=2000,
                    help='Half-width around each breakpoint (bp)')
    p.add_argument('--min-mapq',     type=int, default=20)
    p.add_argument('--min-clip-len', type=int, default=10)
    p.add_argument('--max-samples',  type=int, default=30,
                    help='Max samples to show on plot (top N by support)')
    p.add_argument('--out',          default=None,
                    help='Output path (.pdf or .png). Defaults to the '
                         'evidence_registry figures dir when --candidate '
                         'is used.')
    p.add_argument('--anonymize',    action='store_true',
                    help='Replace sample IDs with sample1, sample2, ... in plot')
    args = p.parse_args()

    # Registry resolution when --candidate is given
    reg = None
    if args.candidate:
        loader_dir = _find_registry_loader_py()
        if loader_dir is None:
            sys.exit('[pileup] ERROR: --candidate requires '
                     'registries/api/python/registry_loader.py (not found). '
                     'Set $BASE or run from the repo root.')
        if loader_dir not in sys.path:
            sys.path.insert(0, loader_dir)
        from registry_loader import load_registry
        reg = load_registry()

        cand = reg.intervals.get_candidate(args.candidate)
        if cand is None:
            sys.exit(f'[pileup] ERROR: candidate {args.candidate} not in '
                     f'interval_registry')
        if args.chrom is None:
            args.chrom = cand['chrom']
        # Prefer boundary_refined_{side}.final_bp; fall back to the raw interval
        if args.bp_left is None:
            refined_L = reg.evidence.read_block(args.candidate,
                                                 'boundary_refined_left')
            if refined_L and refined_L.get('data', {}).get('final_bp'):
                args.bp_left = int(refined_L['data']['final_bp'])
            else:
                args.bp_left = int(cand['start_bp'])
        if args.bp_right is None:
            refined_R = reg.evidence.read_block(args.candidate,
                                                 'boundary_refined_right')
            if refined_R and refined_R.get('data', {}).get('final_bp'):
                args.bp_right = int(refined_R['data']['final_bp'])
            else:
                args.bp_right = int(cand['end_bp'])
        # Sample list from registry
        if args.sample_list is None:
            grs = reg.samples.get_groups_for_candidate(args.candidate)
            if grs is None or not any(grs.get(k) for k in
                                       ('HOM_REF', 'HET', 'HOM_INV')):
                sys.exit(f'[pileup] ERROR: no karyotype groups registered '
                         f'for {args.candidate}; pass --sample-list explicitly '
                         f'or register inv_{args.candidate}_HOM_REF/HET/HOM_INV.')
            reg_root = os.environ.get('REGISTRIES', '') or \
                       os.path.join(os.environ.get('BASE', '.'), 'registries')
            cand_raw = os.path.join(reg_root, 'data', 'evidence_registry',
                                     'per_candidate', args.candidate, 'raw')
            os.makedirs(cand_raw, exist_ok=True)
            auto_list = os.path.join(cand_raw, 'pileup_samples.tsv')
            rows = []
            for gname, label in [('HOM_REF', 'REF'),
                                  ('HET',     'HET'),
                                  ('HOM_INV', 'INV')]:
                for s in grs.get(gname, []) or []:
                    rows.append({'sample': s, 'group': label})
            pd.DataFrame(rows).to_csv(auto_list, sep='\t', index=False)
            args.sample_list = auto_list
            print(f'[pileup] resolved sample list: {len(rows)} samples -> {auto_list}')
        # Default output path under evidence figures/
        if args.out is None:
            reg_root = os.environ.get('REGISTRIES', '') or \
                       os.path.join(os.environ.get('BASE', '.'), 'registries')
            fig_dir = os.path.join(reg_root, 'data', 'evidence_registry',
                                    'per_candidate', args.candidate, 'figures')
            os.makedirs(fig_dir, exist_ok=True)
            args.out = os.path.join(fig_dir, f'panel_D_{args.candidate}.pdf')
        print(f'[pileup] --candidate {args.candidate} resolved: '
              f'chrom={args.chrom} bp_left={args.bp_left} bp_right={args.bp_right}')
        print(f'[pileup]   --out -> {args.out}')

    # Required-fields validation (when --candidate wasn't used)
    missing = [k for k in ('chrom', 'bp_left', 'bp_right', 'sample_list', 'out')
                if getattr(args, k) is None]
    if missing:
        sys.exit(f'[pileup] ERROR: missing required args (and --candidate '
                 f'not used to auto-resolve): {missing}')

    # Load sample list
    sample_df = pd.read_csv(args.sample_list, sep='\t')
    if 'sample' not in sample_df.columns:
        sys.exit('ERROR: --sample-list must have a "sample" column')
    if 'group' not in sample_df.columns:
        sample_df['group'] = 'unknown'

    # Gather per-sample records
    per_sample_records_L = {}
    per_sample_records_R = {}

    if args.evidence_tsv:
        # Mode A: pre-extracted evidence
        print(f'[pileup] mode A: reading {args.evidence_tsv}')
        ev = pd.read_csv(args.evidence_tsv, sep='\t')
        # Expected schema: sample, bp_side (L/R), read_type, read_pos, read_end, mate_pos, strand
        for sid, g in ev.groupby('sample'):
            for side, records_dict in [('L', per_sample_records_L),
                                          ('R', per_sample_records_R)]:
                sub = g[g['bp_side'] == side]
                recs = {'split_reads': [], 'discordant_pairs': [],
                          'soft_clips': [], 'mate_at_partner': [],
                          'coverage_bins': {}}
                for _, row in sub.iterrows():
                    rt = row.get('read_type', '')
                    if rt == 'split_read':
                        recs['split_reads'].append((
                            int(row['read_pos']), int(row.get('read_end', row['read_pos'])),
                            row.get('strand', '+'), int(row.get('mate_pos', 0))
                        ))
                    elif rt == 'discordant':
                        recs['discordant_pairs'].append((
                            int(row['read_pos']), int(row.get('mate_pos', 0)),
                            row.get('orient', 'FF')
                        ))
                    elif rt == 'soft_clip':
                        recs['soft_clips'].append((
                            int(row['read_pos']), int(row.get('clip_len', 10)),
                            row.get('side', 'left')
                        ))
                    elif rt == 'mate_at_partner':
                        recs['mate_at_partner'].append((
                            int(row['read_pos']), int(row.get('mate_pos', 0))
                        ))
                    elif rt == 'coverage':
                        recs['coverage_bins'][int(row['read_pos'])] = \
                            int(row.get('depth', 0))
                records_dict[sid] = recs

    elif args.bam_dir:
        # Mode B: extract fresh from BAMs
        if not HAS_PYSAM:
            sys.exit('ERROR: pysam required for --bam-dir mode. pip install pysam')
        print(f'[pileup] mode B: extracting from BAMs in {args.bam_dir}')

        for _, row in sample_df.iterrows():
            sid = row['sample']
            bam_path = os.path.join(args.bam_dir, sid + args.bam_suffix)
            if not os.path.exists(bam_path):
                print(f'[pileup]   skip {sid}: BAM not found at {bam_path}')
                continue
            per_sample_records_L[sid] = extract_evidence_pysam(
                bam_path, args.chrom, args.bp_left, args.window,
                args.min_mapq, args.min_clip_len, args.chrom, args.bp_right
            )
            per_sample_records_R[sid] = extract_evidence_pysam(
                bam_path, args.chrom, args.bp_right, args.window,
                args.min_mapq, args.min_clip_len, args.chrom, args.bp_left
            )
            n_ev = (len(per_sample_records_L[sid]['split_reads']) +
                    len(per_sample_records_L[sid]['discordant_pairs']) +
                    len(per_sample_records_R[sid]['split_reads']) +
                    len(per_sample_records_R[sid]['discordant_pairs']))
            print(f'[pileup]   {sid}: {n_ev} evidence reads')
    else:
        sys.exit('ERROR: must provide either --evidence-tsv (mode A) or '
                  '--bam-dir (mode B)')

    # Order + clip to max-samples
    samples_in_order = order_samples(sample_df, per_sample_records_L)
    samples_in_order = samples_in_order[:args.max_samples]

    # Optional anonymization of sample IDs
    if args.anonymize:
        display_labels = [f'sample{i + 1}' for i in range(len(samples_in_order))]
        label_map = dict(zip(samples_in_order, display_labels))
    else:
        display_labels = samples_in_order
        label_map = {s: s for s in samples_in_order}

    group_map = dict(zip(sample_df['sample'], sample_df['group']))
    group_colors = {'REF': '#1f5d8e', 'HET': '#8d8d8d', 'INV': '#b03a2e',
                      'unknown': '#333'}

    # ---- Layout ----
    # Each sample gets one row across split-reads + discordants (stacked).
    # Coverage and assembly are aggregated, shown below.
    #
    # Grid: [labels | BP1 tracks | BP2 tracks]
    #         rows:  [stacked_reads (split+disc) | coverage | assembly]
    fig = plt.figure(figsize=(14, max(7, 0.22 * len(samples_in_order) + 3)))
    gs = fig.add_gridspec(
        nrows=3, ncols=3,
        width_ratios=[0.09, 1.0, 1.0],
        height_ratios=[3.5, 0.9, 1.0],
        hspace=0.18, wspace=0.05,
        left=0.08, right=0.98, top=0.93, bottom=0.08
    )

    # Left column: one label per sample, aligned with the stacked-reads row
    ax_labels = fig.add_subplot(gs[0, 0])
    n = len(samples_in_order)
    for i, (orig, disp) in enumerate(zip(samples_in_order, display_labels)):
        y = n - i - 1
        grp = group_map.get(orig, 'unknown')
        color = group_colors.get(grp, '#333')
        ax_labels.text(0.98, y + 0.5, disp, fontsize=6, ha='right',
                        va='center', color=color, family='monospace')
    ax_labels.set_xlim(0, 1)
    ax_labels.set_ylim(-0.5, n + 0.5)
    ax_labels.set_xticks([]); ax_labels.set_yticks([])
    for side in ('top', 'right', 'bottom', 'left'):
        ax_labels.spines[side].set_visible(False)

    # BP1 — stacked reads, coverage, assembly
    ax_L_reads = fig.add_subplot(gs[0, 1])
    ax_L_cov   = fig.add_subplot(gs[1, 1], sharex=ax_L_reads)
    ax_L_asm   = fig.add_subplot(gs[2, 1], sharex=ax_L_reads)

    # BP2 — stacked reads, coverage, assembly
    ax_R_reads = fig.add_subplot(gs[0, 2])
    ax_R_cov   = fig.add_subplot(gs[1, 2], sharex=ax_R_reads)
    ax_R_asm   = fig.add_subplot(gs[2, 2], sharex=ax_R_reads)

    render_panel(
        ax_L_reads, ax_L_cov, ax_L_asm,
        per_sample_records_L, args.bp_left, args.window,
        samples_in_order,
        f'Left breakpoint (BP1, {args.bp_left / 1e6:.2f} Mb)',
        group_map=group_map, group_colors=group_colors
    )
    render_panel(
        ax_R_reads, ax_R_cov, ax_R_asm,
        per_sample_records_R, args.bp_right, args.window,
        samples_in_order,
        f'Right breakpoint (BP2, {args.bp_right / 1e6:.2f} Mb)',
        group_map=group_map, group_colors=group_colors
    )

    # Legend
    legend_handles = [
        mpatches.Patch(color='#1f5d8e', label='Discordant FF (left→right)'),
        mpatches.Patch(color='#b03a2e', label='Discordant RR (alt orient)'),
        mpatches.Patch(color='#c62828', label='Split read / soft clip'),
        mpatches.Patch(color='#bbbbbb', label='Coverage 10-90%'),
        mpatches.Patch(color='#2c5aa0', label='Assembly alignment'),
    ]
    fig.legend(handles=legend_handles, loc='lower center', ncol=5,
                 fontsize=7, frameon=False, bbox_to_anchor=(0.5, 0.01))

    fig.suptitle(
        f'Breakpoint evidence — {args.chrom}  '
        f'({len(samples_in_order)} samples, ordered INV→HET→REF, '
        f'within group by support)',
        fontsize=10, y=0.97
    )

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.suffix.lower() == '.png':
        fig.savefig(out_path, dpi=150, bbox_inches='tight')
    else:
        fig.savefig(out_path, bbox_inches='tight')
    # Also save PNG sibling for quick preview
    if out_path.suffix.lower() == '.pdf':
        fig.savefig(out_path.with_suffix('.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)

    print(f'[pileup] wrote {out_path}')

    # Register the figure path in evidence_registry if we were invoked via
    # --candidate (so reg is populated).
    if reg is not None and args.candidate:
        try:
            reg.evidence.add_evidence(args.candidate,
                                       'figure_panel_D_pdf_path',
                                       str(out_path))
            if out_path.suffix.lower() == '.pdf':
                reg.evidence.add_evidence(args.candidate,
                                           'figure_panel_D_png_path',
                                           str(out_path.with_suffix('.png')))
            print(f'[pileup] registered figure path under candidate {args.candidate}')
        except Exception as e:
            print(f'[pileup] could not register figure path: {e}')


if __name__ == '__main__':
    main()
