#!/usr/bin/env python3
"""
STEP04 – Export weak filtered indel candidates with read-level evidence.

Identifies filtered indels that were NOT rescued in STEP03 (PENDING status)
and extracts read-level evidence using pysam.  These are candidates for
population-level rescue in later steps.

For each weak candidate, exports:
  - VCF-level fields
  - pysam-extracted read evidence: supporting reads, strand balance,
    CIGAR-derived indel counts
  - local event signature fields for cross-sample grouping
  - flanking reference sequence for motif/context matching

Outputs:
  weak_indel_candidates.tsv  – one row per weak filtered indel candidate

Usage:
  python STEP04_export_weak_indel_candidates.py \
      --classified  all_classified_step03.tsv \
      --bam         sample.bam \
      --ref         reference.fa \
      --outdir      <output_directory> \
      --sample_id   CGA009 \
      [--flank_bp   20] \
      [--region_pad 5]
"""

import os, sys, argparse
import numpy as np
import pandas as pd

try:
    import pysam
    HAS_PYSAM = True
except ImportError:
    HAS_PYSAM = False
    print("[WARN] pysam not available – read-level evidence will be NA", file=sys.stderr)


# ── event signature builder ──────────────────────────────────────────────────

def build_event_signature(chrom, pos, ref, alt, var_type, indel_len, left_flank, right_flank):
    """Build a composite event signature for cross-sample grouping.

    Signature = chrom:pos:type:len:norm_event:flank5mer
    where norm_event = normalized ref>alt
    and flank5mer = last 5bp of left flank + first 5bp of right flank
    """
    norm_event = f"{ref}>{alt}"
    fl5 = (left_flank[-5:] if left_flank else "") + "|" + (right_flank[:5] if right_flank else "")
    ilen = int(indel_len) if pd.notna(indel_len) else 0
    sig = f"{chrom}:{pos}:{var_type}:{ilen}:{norm_event}:{fl5}"
    return sig


def build_motif(ref, alt, var_type):
    """Extract the inserted or deleted sequence motif."""
    if var_type == "INS":
        # inserted bases = alt minus shared prefix
        shared = min(len(ref), len(alt))
        return alt[shared - 1:] if shared > 0 else alt
    elif var_type == "DEL":
        shared = min(len(ref), len(alt))
        return ref[shared - 1:] if shared > 0 else ref
    return ""


# ── pysam read-level extraction ─────────────────────────────────────────────

def extract_read_evidence(bam_path, chrom, pos, ref, alt, var_type, indel_len,
                          region_pad=5):
    """Use pysam to inspect reads around an indel candidate.

    Returns dict with:
      - n_reads_spanning: total reads overlapping the position
      - n_reads_support_indel: reads with a matching CIGAR indel at/near pos
      - support_read_names: comma-separated read names (for within-sample debugging)
      - strand_fwd / strand_rev: counts among supporting reads
      - mean_mapq_support: mean MAPQ of supporting reads
    """
    result = {
        "N_READS_SPANNING": np.nan,
        "N_READS_SUPPORT_INDEL": 0,
        "SUPPORT_READ_NAMES": "",
        "STRAND_FWD": 0,
        "STRAND_REV": 0,
        "MEAN_MAPQ_SUPPORT": np.nan,
    }

    if not HAS_PYSAM:
        return result

    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        print(f"[WARN] Cannot open BAM: {e}", file=sys.stderr)
        return result

    # 0-based coordinate for pysam
    pos0 = pos - 1
    start = max(0, pos0 - region_pad)
    end = pos0 + max(len(ref), 1) + region_pad

    try:
        reads = list(bam.fetch(chrom, start, end))
    except Exception:
        reads = []

    result["N_READS_SPANNING"] = len(reads)

    abs_ilen = abs(int(indel_len)) if pd.notna(indel_len) else 1
    target_op = 1 if var_type == "INS" else 2  # BAM_CINS=1, BAM_CDEL=2

    support_names = []
    mapqs = []
    fwd = rev = 0

    for read in reads:
        if read.is_unmapped or read.cigartuples is None:
            continue

        # Walk CIGAR to find indels near our target position
        ref_pos = read.reference_start  # 0-based
        for op, length in read.cigartuples:
            if op in (0, 7, 8):  # M, =, X → consumes ref
                ref_pos += length
            elif op == 1:  # INS → does NOT consume ref
                if op == target_op and abs(ref_pos - pos0) <= region_pad:
                    if abs(length - abs_ilen) <= max(1, abs_ilen // 3):
                        support_names.append(read.query_name)
                        mapqs.append(read.mapping_quality)
                        if read.is_reverse:
                            rev += 1
                        else:
                            fwd += 1
                        break
            elif op == 2:  # DEL → consumes ref
                if op == target_op and abs(ref_pos - pos0) <= region_pad:
                    if abs(length - abs_ilen) <= max(1, abs_ilen // 3):
                        support_names.append(read.query_name)
                        mapqs.append(read.mapping_quality)
                        if read.is_reverse:
                            rev += 1
                        else:
                            fwd += 1
                        ref_pos += length
                        break
                ref_pos += length
            elif op == 3:  # N → ref skip
                ref_pos += length
            elif op == 4:  # soft clip
                pass
            elif op == 5:  # hard clip
                pass

    result["N_READS_SUPPORT_INDEL"] = len(support_names)
    result["SUPPORT_READ_NAMES"] = ",".join(support_names[:50])  # cap at 50
    result["STRAND_FWD"] = fwd
    result["STRAND_REV"] = rev
    if mapqs:
        result["MEAN_MAPQ_SUPPORT"] = np.mean(mapqs)

    bam.close()
    return result


def get_flanking_seq(ref_fa, chrom, pos, flank_bp=20):
    """Extract flanking reference sequence around a position."""
    if ref_fa is None:
        return "", ""
    try:
        pos0 = pos - 1
        left_start = max(0, pos0 - flank_bp)
        left = ref_fa.fetch(chrom, left_start, pos0).upper()
        right_end = min(ref_fa.get_reference_length(chrom), pos0 + flank_bp)
        right = ref_fa.fetch(chrom, pos0, right_end).upper()
        return left, right
    except Exception:
        return "", ""


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description="STEP04: Export weak indel candidates with read evidence")
    ap.add_argument("--classified", required=True, help="all_classified_step03.tsv")
    ap.add_argument("--bam", required=True, help="BAM file for this sample")
    ap.add_argument("--ref", required=True, help="Reference FASTA")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--sample_id", required=True)
    ap.add_argument("--flank_bp", type=int, default=20)
    ap.add_argument("--region_pad", type=int, default=5)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.classified, sep="\t")

    # Select weak filtered indels = PENDING + IS_INDEL
    weak = df[(df["FINAL_CLASS"] == "PENDING") & (df["IS_INDEL"] == 1)].copy()
    print(f"[STEP04] Weak filtered indel candidates: {len(weak)}")

    if len(weak) == 0:
        print("[STEP04] No weak candidates to export.")
        out = os.path.join(args.outdir, "weak_indel_candidates.tsv")
        pd.DataFrame().to_csv(out, sep="\t", index=False)
        return

    # Open reference FASTA
    ref_fa = None
    if HAS_PYSAM:
        try:
            ref_fa = pysam.FastaFile(args.ref)
        except Exception as e:
            print(f"[WARN] Cannot open reference: {e}", file=sys.stderr)

    # Process each weak candidate
    evidence_rows = []
    for i, (idx, row) in enumerate(weak.iterrows()):
        if i > 0 and i % 500 == 0:
            print(f"[STEP04]   processed {i}/{len(weak)} candidates …")

        chrom = row["CHROM"]
        pos = int(row["POS"])
        ref_allele = str(row["REF"])
        alt_allele = str(row["ALT1"])
        var_type = row["VAR_TYPE"]
        indel_len = row["INDEL_LEN"]

        # flanking sequence
        left_flank, right_flank = get_flanking_seq(ref_fa, chrom, pos, args.flank_bp)

        # read evidence
        read_ev = extract_read_evidence(
            args.bam, chrom, pos, ref_allele, alt_allele,
            var_type, indel_len, region_pad=args.region_pad
        )

        # motif
        motif = build_motif(ref_allele, alt_allele, var_type)

        # event signature
        sig = build_event_signature(
            chrom, pos, ref_allele, alt_allele,
            var_type, indel_len, left_flank, right_flank
        )

        evidence_rows.append({
            "SAMPLE_ID": args.sample_id,
            "CHROM": chrom,
            "POS": pos,
            "END": int(row["END"]),
            "REF": ref_allele,
            "ALT": alt_allele,
            "VAR_TYPE": var_type,
            "INDEL_LEN": indel_len,
            "ABS_INDEL_LEN": row.get("ABS_INDEL_LEN", np.nan),
            "QUAL": row["QUAL"],
            "GQ": row["GQ"],
            "DP": row["DP"],
            "AD_REF": row["AD_REF"],
            "AD_ALT": row["AD_ALT"],
            "AF1": row["AF1"],
            "GT": row["GT"],
            "GT_CLASS": row["GT_CLASS"],
            "FILTER": row["FILTER"],
            "STATUS": row["STATUS"],
            "CALLER_SOURCE": row.get("CALLER_SOURCE", ""),
            "BLOCK_ID": row.get("BLOCK_ID", -1),
            "OVERLAP_CLUSTER_ID": row.get("OVERLAP_CLUSTER_ID", -1),
            "N_LOCAL_VAR_10BP": row.get("N_LOCAL_VAR_10BP", 0),
            "N_LOCAL_VAR_20BP": row.get("N_LOCAL_VAR_20BP", 0),
            # read evidence
            "N_READS_SPANNING": read_ev["N_READS_SPANNING"],
            "N_READS_SUPPORT_INDEL": read_ev["N_READS_SUPPORT_INDEL"],
            "SUPPORT_READ_NAMES": read_ev["SUPPORT_READ_NAMES"],
            "STRAND_FWD": read_ev["STRAND_FWD"],
            "STRAND_REV": read_ev["STRAND_REV"],
            "MEAN_MAPQ_SUPPORT": read_ev["MEAN_MAPQ_SUPPORT"],
            # motif / flank / signature
            "INDEL_MOTIF": motif,
            "LEFT_FLANK_20BP": left_flank,
            "RIGHT_FLANK_20BP": right_flank,
            "EVENT_SIGNATURE": sig,
            # class
            "FINAL_CLASS": "WEAK_SIGNATURE_CANDIDATE",
            "FINAL_LAYER": "weak_signature_candidate",
            "RESCUE_REASON": "weak_signature_candidate",
        })

    if ref_fa:
        ref_fa.close()

    out_df = pd.DataFrame(evidence_rows)
    out_path = os.path.join(args.outdir, "weak_indel_candidates.tsv")
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"[STEP04] Written: {out_path}  ({len(out_df)} rows)")

    # quick summary
    print(f"[STEP04] Support reads distribution:")
    if len(out_df) > 0:
        for n in [0, 1, 2, 3]:
            cnt = (out_df["N_READS_SUPPORT_INDEL"] == n).sum()
            print(f"         {n} supporting reads: {cnt}")
        cnt_gt3 = (out_df["N_READS_SUPPORT_INDEL"] > 3).sum()
        print(f"         >3 supporting reads: {cnt_gt3}")


if __name__ == "__main__":
    main()
