#!/usr/bin/env python3
"""
cheat20_extract_junctions.py — Breakpoint junction microhomology extraction

BIOLOGY:
  The actual DNA sequence at the breakpoint junction reveals the formation
  mechanism: microhomology (1-9 bp) → MMEJ/alt-EJ, blunt → classical NHEJ,
  insertion → template switching, duplication → staggered break.

USAGE:
  python cheat20_extract_junctions.py --bam <bam> --chr <chr> --bp <pos> \
    --ref <fasta> --out <tsv>

REQUIRES: pysam
"""

import argparse
import sys
from collections import Counter

try:
    import pysam
except ImportError:
    print("ERROR: pysam required. Install: pip install pysam", file=sys.stderr)
    sys.exit(1)


def extract_junction_reads(bam_path, chr, bp, window=500):
    """Collect split reads and soft-clipped reads at a breakpoint."""
    junctions = []
    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        print(f"[cheat20] Cannot open BAM: {e}", file=sys.stderr)
        return junctions

    region_start = max(0, bp - window)
    region_end = bp + window

    for read in bam.fetch(chr, region_start, region_end):
        if read.is_unmapped or read.is_duplicate or read.mapping_quality < 10:
            continue

        cigar = read.cigartuples
        if cigar is None:
            continue

        # Check for soft clips near breakpoint
        ref_pos = read.reference_start
        for op, length in cigar:
            if op == 4:  # soft clip
                clip_dist = abs(ref_pos - bp)
                if clip_dist <= window:
                    seq_at_clip = read.query_sequence
                    if seq_at_clip:
                        junctions.append({
                            'read_name': read.query_name,
                            'type': 'soft_clip',
                            'clip_pos': ref_pos,
                            'clip_length': length,
                            'distance_to_bp': clip_dist,
                            'mapq': read.mapping_quality,
                            'strand': '-' if read.is_reverse else '+',
                        })
            if op in (0, 1, 2, 3, 7, 8):  # consuming ref
                if op != 1:  # not insertion
                    ref_pos += length

        # Check for supplementary alignment (split read)
        if read.has_tag('SA'):
            sa_tag = read.get_tag('SA')
            for sa in sa_tag.strip(';').split(';'):
                parts = sa.split(',')
                if len(parts) >= 4:
                    sa_chr = parts[0]
                    sa_pos = int(parts[1])
                    sa_strand = parts[2]
                    if sa_chr == chr and abs(sa_pos - bp) <= window * 5:
                        junctions.append({
                            'read_name': read.query_name,
                            'type': 'split_read',
                            'clip_pos': read.reference_start,
                            'sa_pos': sa_pos,
                            'sa_strand': sa_strand,
                            'distance_to_bp': min(
                                abs(read.reference_start - bp),
                                abs(sa_pos - bp)),
                            'mapq': read.mapping_quality,
                            'strand': '-' if read.is_reverse else '+',
                        })

    bam.close()
    return junctions


def classify_junction(junctions, ref_fasta, chr, bp, window=50):
    """Classify junction type from collected reads."""
    if not junctions:
        return {
            'junction_class': 'NO_EVIDENCE',
            'microhomology_length': 0,
            'insertion_length': 0,
            'n_supporting_reads': 0,
            'confidence': 0.0,
        }

    n_reads = len(junctions)
    split_reads = [j for j in junctions if j['type'] == 'split_read']
    clip_reads = [j for j in junctions if j['type'] == 'soft_clip']

    # Extract reference sequence around breakpoint for microhomology check
    microhom_len = 0
    try:
        ref = pysam.FastaFile(ref_fasta)
        left_seq = ref.fetch(chr, max(0, bp - 20), bp).upper()
        right_seq = ref.fetch(chr, bp, bp + 20).upper()
        ref.close()

        # Check microhomology: shared bases at junction
        for i in range(min(len(left_seq), len(right_seq), 20)):
            if i < len(left_seq) and i < len(right_seq):
                if left_seq[-(i+1):] == right_seq[:i+1] if i > 0 else True:
                    # Simple overlap check
                    pass
        # Count matching bases from breakpoint outward on both sides
        for i in range(min(20, len(left_seq), len(right_seq))):
            li = len(left_seq) - 1 - i
            if li >= 0 and i < len(right_seq) and left_seq[li] == right_seq[i]:
                microhom_len += 1
            else:
                break
    except Exception:
        pass

    # Classify
    if microhom_len == 0:
        junction_class = 'BLUNT'
    elif microhom_len <= 3:
        junction_class = 'SHORT_MICROHOMOLOGY'
    elif microhom_len <= 9:
        junction_class = 'MICROHOMOLOGY'
    else:
        junction_class = 'LONG_HOMOLOGY'

    confidence = min(1.0, n_reads / 10.0)

    return {
        'junction_class': junction_class,
        'microhomology_length': microhom_len,
        'insertion_length': 0,
        'n_supporting_reads': n_reads,
        'n_split_reads': len(split_reads),
        'n_clip_reads': len(clip_reads),
        'confidence': round(confidence, 3),
    }


def main():
    parser = argparse.ArgumentParser(
        description='Extract and classify breakpoint junction reads')
    parser.add_argument('--bam', required=True, help='Markdup BAM file')
    parser.add_argument('--chr', required=True, help='Chromosome')
    parser.add_argument('--bp', required=True, type=int,
                        help='Breakpoint position')
    parser.add_argument('--ref', required=True, help='Reference FASTA')
    parser.add_argument('--out', required=True, help='Output TSV path')
    parser.add_argument('--window', type=int, default=500,
                        help='Window around breakpoint (default: 500)')
    args = parser.parse_args()

    print(f"[cheat20] Extracting junctions at {args.chr}:{args.bp} "
          f"from {args.bam}", file=sys.stderr)

    junctions = extract_junction_reads(args.bam, args.chr, args.bp, args.window)
    print(f"[cheat20] Found {len(junctions)} junction reads", file=sys.stderr)

    result = classify_junction(junctions, args.ref, args.chr, args.bp)
    print(f"[cheat20] Classification: {result['junction_class']} "
          f"(microhom={result['microhomology_length']}bp, "
          f"confidence={result['confidence']})", file=sys.stderr)

    # Write output
    with open(args.out, 'w') as f:
        f.write('\t'.join(result.keys()) + '\n')
        f.write('\t'.join(str(v) for v in result.values()) + '\n')

    print(f"[cheat20] Written → {args.out}", file=sys.stderr)


if __name__ == '__main__':
    main()
