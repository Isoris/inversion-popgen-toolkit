#!/usr/bin/env python3
"""
STEP07A – Per-sample regenotyping (SLURM array worker)

Scans ONE sample's BAM against the shared regenotype catalog.
Writes one TSV per sample: regenotyped_rescued_indels.<SAMPLE>.tsv

This is the parallelizable part of STEP07.
After all samples complete, run STEP07B to merge + build cohort summary.

Usage:
  python STEP07A_regenotype_one_sample.py \\
      --catalog    rescued_indel_regenotype_catalog.tsv \\
      --sample_id  CGA097 \\
      --bam        /path/to/CGA097.bam \\
      --outdir     <population_output_dir> \\
      [--chrom C_gar_LG01] \\
      [--pad_bp 50]
"""

import os, sys, argparse, time
import numpy as np
import pandas as pd

try:
    import pysam
except ImportError:
    print("[ERROR] pysam is required for STEP07A", file=sys.stderr)
    sys.exit(1)


def _now(): return time.time()
def _elapsed(t0):
    dt = time.time() - t0
    return f"{dt:.1f}s" if dt < 60 else f"{dt/60:.1f}min"


# ── Candidate index (same as original STEP07) ───────────────────────────────

class CandidateIndex:
    def __init__(self, catalog_df, pad_bp=50):
        self.n = len(catalog_df)
        if self.n == 0:
            self.starts = np.array([], dtype=np.int64)
            self.ends = np.array([], dtype=np.int64)
            self.candidates = []
            return

        cat = catalog_df.sort_values("POS_REPRESENTATIVE").reset_index(drop=True)
        self.starts = np.empty(self.n, dtype=np.int64)
        self.ends = np.empty(self.n, dtype=np.int64)
        self.candidates = []

        for i, (_, row) in enumerate(cat.iterrows()):
            pos = int(row["POS_REPRESENTATIVE"])
            ref_allele = str(row["REF"])
            pos0 = pos - 1
            ref_len = max(len(ref_allele), 1)
            self.starts[i] = max(0, pos0 - pad_bp)
            self.ends[i] = pos0 + ref_len + pad_bp

            var_type = str(row["VAR_TYPE"])
            abs_ilen = int(row.get("ABS_INDEL_LEN", 1)) if pd.notna(row.get("ABS_INDEL_LEN", None)) else 1

            self.candidates.append({
                "idx": i,
                "regen_id": row["REGEN_ID"],
                "chrom": row["CHROM"],
                "pos": pos,
                "pos0": pos0,
                "ref": ref_allele,
                "alt": str(row["ALT"]),
                "var_type": var_type,
                "target_op": 1 if var_type == "INS" else 2,
                "abs_ilen": abs_ilen,
                "ilen_tol": max(1, abs_ilen // 3),
            })

    def find_overlapping(self, read_start, read_end):
        if self.n == 0:
            return []
        hi = np.searchsorted(self.starts, read_end, side='left')
        if hi == 0:
            return []
        mask = self.ends[:hi] > read_start
        return np.where(mask)[0].tolist()


# ── CIGAR evidence checker ───────────────────────────────────────────────────

REGION_PAD = 5

def check_read_evidence(read, cand):
    if read.cigartuples is None:
        return "none"
    pos0 = cand["pos0"]
    target_op = cand["target_op"]
    ilen = cand["abs_ilen"]
    ilen_tol = cand["ilen_tol"]
    ref_pos = read.reference_start
    found_indel = False
    found_ref_span = False

    for op, length in read.cigartuples:
        if op in (0, 7, 8):
            if ref_pos <= pos0 < ref_pos + length:
                found_ref_span = True
            ref_pos += length
        elif op == 1:
            if target_op == 1 and abs(ref_pos - pos0) <= REGION_PAD:
                if abs(length - ilen) <= ilen_tol:
                    found_indel = True
        elif op == 2:
            if target_op == 2 and abs(ref_pos - pos0) <= REGION_PAD:
                if abs(length - ilen) <= ilen_tol:
                    found_indel = True
            ref_pos += length
        elif op == 3:
            ref_pos += length
        elif op in (4, 5):
            pass

    if found_indel:
        return "support"
    elif found_ref_span:
        return "ref"
    return "none"


# ── Per-sample linear pass ───────────────────────────────────────────────────

def regenotype_sample_linear(bam_path, cand_index, chrom=None):
    accum = {}
    for c in cand_index.candidates:
        accum[c["idx"]] = {
            "n_spanning": 0, "n_support": 0, "n_ref": 0,
            "fwd": 0, "rev": 0,
            "mapqs_all": [], "mapqs_support": [],
        }

    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        print(f"[STEP07A] WARN: cannot open BAM: {e}")
        return accum

    try:
        iterator = bam.fetch(chrom) if chrom else bam.fetch()
        for read in iterator:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.reference_start is None or read.cigartuples is None:
                continue
            rstart = read.reference_start
            rend = read.reference_end
            if rend is None:
                continue

            hit_indices = cand_index.find_overlapping(rstart, rend)
            if not hit_indices:
                continue

            mq = read.mapping_quality
            is_rev = read.is_reverse

            for ci in hit_indices:
                cand = cand_index.candidates[ci]
                acc = accum[ci]
                acc["n_spanning"] += 1
                acc["mapqs_all"].append(mq)
                ev = check_read_evidence(read, cand)
                if ev == "support":
                    acc["n_support"] += 1
                    acc["mapqs_support"].append(mq)
                    if is_rev:
                        acc["rev"] += 1
                    else:
                        acc["fwd"] += 1
                elif ev == "ref":
                    acc["n_ref"] += 1
    except Exception as e:
        print(f"[STEP07A] WARN: BAM iteration error: {e}")
    finally:
        bam.close()

    return accum


def accum_to_result(acc):
    n_support = acc["n_support"]
    n_ref = acc["n_ref"]
    total_informative = n_support + n_ref

    result = {
        "N_READS_SPANNING": acc["n_spanning"],
        "N_READS_SUPPORT": n_support,
        "N_READS_REF": n_ref,
        "STRAND_FWD": acc["fwd"],
        "STRAND_REV": acc["rev"],
        "MEAN_MAPQ_ALL": round(float(np.mean(acc["mapqs_all"])), 1) if acc["mapqs_all"] else np.nan,
        "MEAN_MAPQ_SUPPORT": round(float(np.mean(acc["mapqs_support"])), 1) if acc["mapqs_support"] else np.nan,
        "SUPPORT_FRACTION": 0.0,
        "REGEN_GT": "./.",
        "REGEN_CONFIDENCE": "NO_DATA",
    }

    if total_informative > 0:
        sf = n_support / total_informative
        result["SUPPORT_FRACTION"] = round(sf, 4)
        if sf >= 0.8:
            result["REGEN_GT"] = "1/1"
            result["REGEN_CONFIDENCE"] = "HIGH" if n_support >= 3 else "LOW"
        elif sf >= 0.15:
            result["REGEN_GT"] = "0/1"
            result["REGEN_CONFIDENCE"] = "HIGH" if n_support >= 2 else "LOW"
        else:
            result["REGEN_GT"] = "0/0"
            result["REGEN_CONFIDENCE"] = "HIGH" if total_informative >= 3 else "LOW"

    return result


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    t0 = _now()

    ap = argparse.ArgumentParser(description="STEP07A: Regenotype one sample")
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--sample_id", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--pad_bp", type=int, default=50)
    ap.add_argument("--chrom", default=None)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Load catalog
    cat = pd.read_csv(args.catalog, sep="\t")
    if args.chrom:
        cat = cat[cat["CHROM"] == args.chrom].copy()

    print(f"[STEP07A] Sample={args.sample_id}  Candidates={len(cat)}  Chrom={args.chrom or 'all'}")

    if len(cat) == 0:
        print("[STEP07A] No candidates, writing empty output.")
        out_path = os.path.join(args.outdir, f"regenotyped_rescued_indels.{args.sample_id}.tsv")
        pd.DataFrame().to_csv(out_path, sep="\t", index=False)
        return

    # Build index
    cand_index = CandidateIndex(cat, pad_bp=args.pad_bp)
    print(f"[STEP07A] Index built: {cand_index.n} candidates ({_elapsed(t0)})")

    # Scan BAM
    t_scan = _now()
    accum = regenotype_sample_linear(args.bam, cand_index, chrom=args.chrom)
    print(f"[STEP07A] BAM scan done: {_elapsed(t_scan)}")

    # Build rows
    rows = []
    n_alt = 0
    for c in cand_index.candidates:
        result = accum_to_result(accum[c["idx"]])
        row = {
            "REGEN_ID": c["regen_id"],
            "CHROM": c["chrom"],
            "POS": c["pos"],
            "REF": c["ref"],
            "ALT": c["alt"],
            "VAR_TYPE": c["var_type"],
            "ABS_INDEL_LEN": c["abs_ilen"],
            "SAMPLE_ID": args.sample_id,
        }
        row.update(result)
        rows.append(row)
        if result["REGEN_GT"] in ("0/1", "1/1"):
            n_alt += 1

    # Write per-sample output
    out_path = os.path.join(args.outdir, f"regenotyped_rescued_indels.{args.sample_id}.tsv")
    pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)

    print(f"[STEP07A] {args.sample_id}: {n_alt} alt calls -> {out_path}")
    print(f"[STEP07A] Total time: {_elapsed(t0)}")


if __name__ == "__main__":
    main()
