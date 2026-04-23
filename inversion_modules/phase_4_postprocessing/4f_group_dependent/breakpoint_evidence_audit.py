#!/usr/bin/env python3
# =============================================================================
# breakpoint_evidence_audit.py
# =============================================================================
#
# POST-PROCESSING AUDIT: For each confirmed inversion candidate, counts
# read-level breakpoint evidence per sample WITHOUT modifying genotypes.
#
# PURPOSE:
#   DELLY's PE >= 4 filter controls false positives well, but it also
#   calls some true carriers as 0/0 due to stochastic dropout at 9x coverage.
#   We keep the DELLY filter as-is. This audit documents the evidence
#   landscape per inversion, showing HOW MUCH support the breakpoint has
#   across the population, including sub-threshold evidence.
#
# OUTPUT: Per-candidate audit, NOT per-sample genotype corrections.
#
# FOR EACH CANDIDATE:
#   - Count PCA carriers with strong SV support (GT=0/1 or 1/1)
#   - Count PCA carriers with weak support (GT=0/0 but 1-3 alt reads)
#   - Count PCA carriers with no support (GT=0/0 and 0 alt reads)
#   - Count non-carriers with unexpected support (rare)
#
# INPUT:
#   - merged regenotyped VCF (from DELLY regenotyping step, already done)
#   - PCA carrier lists per candidate (from decomposition output)
#   - Candidate coordinates
#
# OUTPUT:
#   - breakpoint_audit.tsv: one row per candidate with evidence distribution
#   - Registry keys (q7b_pca_carriers_strong_sv, etc.) added to each candidate
#
# Usage:
#   python breakpoint_evidence_audit.py \
#     --vcf catalog_226.INV.raw.vcf.gz \
#     --candidates candidate_regions.tsv \
#     --pca-carriers pca_carrier_assignments.tsv \
#     --outdir breakpoint_audit/
# =============================================================================

import argparse
import gzip
import os
import sys
from collections import defaultdict


def open_text(path):
    """Open .gz or plain text for reading."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_format_field(format_str, sample_str):
    """Parse VCF FORMAT/sample into dict."""
    fmt = format_str.split(":")
    vals = sample_str.split(":")
    return dict(zip(fmt, vals))


def extract_alt_evidence(sample_data, fmt_keys):
    """
    Extract alt-read count from FORMAT fields.
    Handles DELLY (DV, RV) and Manta (PR, SR) formats.

    Returns:
      n_alt: total alt-supporting reads (DV + RV for DELLY, PR_alt + SR_alt for Manta)
      n_ref: total ref-supporting reads
      gt:    genotype string (0/0, 0/1, 1/1, ./.)
    """
    parsed = parse_format_field(":".join(fmt_keys), sample_data)
    gt = parsed.get("GT", "./.")
    n_alt = 0
    n_ref = 0

    # DELLY format
    if "DV" in parsed and "DR" in parsed:
        try:
            n_alt += int(parsed["DV"])
            n_ref += int(parsed["DR"])
        except (ValueError, TypeError):
            pass
    if "RV" in parsed and "RR" in parsed:
        try:
            n_alt += int(parsed["RV"])
            n_ref += int(parsed["RR"])
        except (ValueError, TypeError):
            pass

    # Manta format
    if "PR" in parsed:
        try:
            pr = parsed["PR"].split(",")
            n_ref += int(pr[0])
            n_alt += int(pr[1])
        except (ValueError, TypeError, IndexError):
            pass
    if "SR" in parsed:
        try:
            sr = parsed["SR"].split(",")
            n_ref += int(sr[0])
            n_alt += int(sr[1])
        except (ValueError, TypeError, IndexError):
            pass

    return n_alt, n_ref, gt


def classify_evidence(gt, n_alt):
    """
    Classify a sample's breakpoint evidence:
      strong   = GT is 0/1 or 1/1 (DELLY called carrier)
      weak     = GT is 0/0 but n_alt >= 1 (dropout suspected)
      absent   = GT is 0/0 and n_alt = 0 (no evidence at all)
      missing  = GT is ./. (no data at this site)
    """
    if gt in ("./.", ".|.", "."):
        return "missing"
    if gt in ("0/1", "1/0", "1/1", "0|1", "1|0", "1|1"):
        return "strong"
    if gt in ("0/0", "0|0"):
        if n_alt >= 1:
            return "weak"
        else:
            return "absent"
    return "missing"


def load_pca_carriers(path):
    """
    Load PCA carrier assignments per candidate.
    Expected format: candidate_id<TAB>sample_id<TAB>class

    Returns: dict of {candidate_id: {sample_id: class}}
    """
    carriers = defaultdict(dict)
    with open_text(path) as f:
        header = f.readline().strip().split("\t")
        try:
            idx_cid = header.index("candidate_id")
            idx_sid = header.index("sample_id")
            idx_cls = header.index("class")
        except ValueError:
            sys.exit(f"ERROR: expected columns candidate_id/sample_id/class in {path}")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < max(idx_cid, idx_sid, idx_cls) + 1:
                continue
            carriers[parts[idx_cid]][parts[idx_sid]] = parts[idx_cls]
    return carriers


def load_candidates(path):
    """
    Load candidate table.
    Expected columns: candidate_id, chrom, start_bp, end_bp

    Returns: list of dicts
    """
    rows = []
    with open_text(path) as f:
        header = f.readline().strip().split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                continue
            rows.append({
                "candidate_id": parts[idx["candidate_id"]],
                "chrom": parts[idx["chrom"]],
                "start_bp": int(parts[idx.get("start_bp", idx.get("start"))]),
                "end_bp": int(parts[idx.get("end_bp", idx.get("end"))])
            })
    return rows


def audit_candidate(vcf_path, candidate, pca_carriers, bp_window=5000):
    """
    For one candidate, find SV sites within bp_window of breakpoints and
    aggregate evidence per sample.

    Returns: dict of audit results
    """
    chrom = candidate["chrom"]
    bp_left = candidate["start_bp"]
    bp_right = candidate["end_bp"]

    # Track best evidence per sample (best = strongest across all matching sites)
    # Precedence: strong > weak > absent > missing
    precedence = {"strong": 3, "weak": 2, "absent": 1, "missing": 0}
    sample_best_class = {}
    sample_best_n_alt = {}
    sample_ids_from_vcf = []
    n_sites_found = 0

    with open_text(vcf_path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                sample_ids_from_vcf = header[9:]
                continue

            parts = line.rstrip("\n").split("\t")
            if parts[0] != chrom:
                continue

            pos = int(parts[1])

            # Check if position is near one of the candidate breakpoints
            near_left = abs(pos - bp_left) <= bp_window
            near_right = abs(pos - bp_right) <= bp_window

            # Also check END (stored in INFO for DELLY)
            end_pos = None
            info = parts[7]
            for field in info.split(";"):
                if field.startswith("END="):
                    try:
                        end_pos = int(field.split("=")[1])
                    except ValueError:
                        pass
                    break

            end_near_left = end_pos and abs(end_pos - bp_left) <= bp_window
            end_near_right = end_pos and abs(end_pos - bp_right) <= bp_window

            if not (near_left or near_right or end_near_left or end_near_right):
                continue

            # Relevant site
            n_sites_found += 1
            fmt_keys = parts[8].split(":")

            for i, sdata in enumerate(parts[9:]):
                if i >= len(sample_ids_from_vcf):
                    break
                sid = sample_ids_from_vcf[i]
                n_alt, n_ref, gt = extract_alt_evidence(sdata, fmt_keys)
                cls = classify_evidence(gt, n_alt)

                prev_cls = sample_best_class.get(sid, "missing")
                if precedence[cls] > precedence[prev_cls]:
                    sample_best_class[sid] = cls
                    sample_best_n_alt[sid] = n_alt

    # Tabulate: cross PCA class with SV evidence class
    pca_assignments = pca_carriers.get(candidate["candidate_id"], {})

    counts = {
        "pca_carriers_strong_sv": 0,       # PCA carrier + SV strong
        "pca_carriers_weak_sv": 0,         # PCA carrier + SV weak
        "pca_carriers_no_sv": 0,           # PCA carrier + SV absent
        "pca_carriers_missing_sv": 0,      # PCA carrier + no SV site data
        "pca_ref_unexpected_sv": 0,        # PCA REF but SV shows evidence (rare)
        "pca_ref_absent_sv": 0,            # PCA REF + SV absent (expected)
        "pca_ref_missing_sv": 0,           # PCA REF + no SV data
        "n_pca_carriers": 0,
        "n_pca_ref": 0,
        "n_sv_sites_found": n_sites_found
    }

    for sid, pca_class in pca_assignments.items():
        sv_class = sample_best_class.get(sid, "missing")
        is_carrier = pca_class in ("HET", "HOM_INV", "RECOMBINANT")

        if is_carrier:
            counts["n_pca_carriers"] += 1
            if sv_class == "strong":
                counts["pca_carriers_strong_sv"] += 1
            elif sv_class == "weak":
                counts["pca_carriers_weak_sv"] += 1
            elif sv_class == "absent":
                counts["pca_carriers_no_sv"] += 1
            else:
                counts["pca_carriers_missing_sv"] += 1
        else:
            counts["n_pca_ref"] += 1
            if sv_class in ("strong", "weak"):
                counts["pca_ref_unexpected_sv"] += 1
            elif sv_class == "absent":
                counts["pca_ref_absent_sv"] += 1
            else:
                counts["pca_ref_missing_sv"] += 1

    # Compute summary fractions
    n_car = counts["n_pca_carriers"]
    if n_car > 0:
        counts["pca_carrier_sv_support_pct"] = counts["pca_carriers_strong_sv"] / n_car
        counts["dropout_suspected_fraction"] = counts["pca_carriers_weak_sv"] / n_car
        counts["true_absence_fraction"] = counts["pca_carriers_no_sv"] / n_car
    else:
        counts["pca_carrier_sv_support_pct"] = 0.0
        counts["dropout_suspected_fraction"] = 0.0
        counts["true_absence_fraction"] = 0.0

    counts["candidate_id"] = candidate["candidate_id"]
    counts["chrom"] = chrom
    counts["bp_left"] = bp_left
    counts["bp_right"] = bp_right
    return counts


def write_audit_table(audit_rows, outpath):
    """Write per-candidate audit as TSV."""
    columns = [
        "candidate_id", "chrom", "bp_left", "bp_right",
        "n_pca_carriers", "n_pca_ref", "n_sv_sites_found",
        "pca_carriers_strong_sv", "pca_carriers_weak_sv",
        "pca_carriers_no_sv", "pca_carriers_missing_sv",
        "pca_ref_unexpected_sv", "pca_ref_absent_sv", "pca_ref_missing_sv",
        "pca_carrier_sv_support_pct", "dropout_suspected_fraction",
        "true_absence_fraction"
    ]
    with open(outpath, "w") as f:
        f.write("\t".join(columns) + "\n")
        for row in audit_rows:
            f.write("\t".join(str(row.get(c, "")) for c in columns) + "\n")


def write_registry_keys(audit_rows, outpath):
    """
    Write registry-compatible keys (candidate_id, key, value).
    Each candidate gets ~7 q7b_* keys.
    """
    with open(outpath, "w") as f:
        f.write("candidate_id\tkey\tvalue\tsource\n")
        for row in audit_rows:
            cid = row["candidate_id"]
            key_map = {
                "q7b_pca_carriers_strong_sv": row["pca_carriers_strong_sv"],
                "q7b_pca_carriers_weak_sv": row["pca_carriers_weak_sv"],
                "q7b_pca_carriers_no_sv": row["pca_carriers_no_sv"],
                "q7b_pca_carrier_sv_support_pct": round(row["pca_carrier_sv_support_pct"], 3),
                "q7b_dropout_suspected_fraction": round(row["dropout_suspected_fraction"], 3),
                "q7b_true_absence_fraction": round(row["true_absence_fraction"], 3),
                "q7b_n_sv_sites_audited": row["n_sv_sites_found"]
            }
            for key, value in key_map.items():
                f.write(f"{cid}\t{key}\t{value}\tbreakpoint_evidence_audit.py\n")


def write_interpretation(audit_rows, outpath):
    """Write a human-readable interpretation for each candidate."""
    with open(outpath, "w") as f:
        f.write("=" * 78 + "\n")
        f.write("BREAKPOINT EVIDENCE AUDIT — PER-CANDIDATE INTERPRETATION\n")
        f.write("=" * 78 + "\n\n")
        for row in audit_rows:
            cid = row["candidate_id"]
            n_car = row["n_pca_carriers"]
            n_strong = row["pca_carriers_strong_sv"]
            n_weak = row["pca_carriers_weak_sv"]
            n_none = row["pca_carriers_no_sv"]
            pct = row["pca_carrier_sv_support_pct"] * 100

            f.write(f"Candidate: {cid}\n")
            f.write(f"  Region: {row['chrom']}:{row['bp_left']}-{row['bp_right']}\n")
            f.write(f"  PCA carriers: {n_car}\n")
            f.write(f"  SV sites audited within 5kb of breakpoints: {row['n_sv_sites_found']}\n")
            f.write(f"\n")
            f.write(f"  Breakpoint evidence distribution across PCA carriers:\n")
            f.write(f"    Strong (DELLY GT=0/1 or 1/1):       {n_strong:3d} ({n_strong/max(n_car,1)*100:.0f}%)\n")
            f.write(f"    Weak (GT=0/0 but 1-3 alt reads):    {n_weak:3d} ({n_weak/max(n_car,1)*100:.0f}%)  <- dropout\n")
            f.write(f"    Absent (GT=0/0 and 0 alt reads):    {n_none:3d} ({n_none/max(n_car,1)*100:.0f}%)  <- mask/miss\n")
            f.write(f"\n")
            if pct >= 70:
                f.write(f"  Interpretation: STRONG breakpoint support ({pct:.0f}% of PCA carriers).\n")
            elif pct >= 50:
                f.write(f"  Interpretation: MODERATE breakpoint support ({pct:.0f}%).\n")
            else:
                f.write(f"  Interpretation: WEAK breakpoint support ({pct:.0f}%). Check candidate quality.\n")

            if row["pca_ref_unexpected_sv"] > 0:
                f.write(f"  NOTE: {row['pca_ref_unexpected_sv']} PCA-REF samples show unexpected SV evidence.\n")
                f.write(f"        Possible recombinants, nested inversion carriers, or FP calls.\n")
            f.write("\n" + "-" * 78 + "\n\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, help="Regenotyped VCF (from DELLY merged output)")
    ap.add_argument("--candidates", required=True, help="Candidate regions TSV")
    ap.add_argument("--pca-carriers", required=True, help="PCA class assignments TSV")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--bp-window", type=int, default=5000,
                    help="Search window around breakpoints (bp, default 5000)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print(f"[audit] Loading candidates from {args.candidates}")
    candidates = load_candidates(args.candidates)
    print(f"[audit] Loaded {len(candidates)} candidates")

    print(f"[audit] Loading PCA carrier assignments from {args.pca_carriers}")
    pca_carriers = load_pca_carriers(args.pca_carriers)
    print(f"[audit] Loaded assignments for {len(pca_carriers)} candidates")

    print(f"[audit] Auditing breakpoint evidence from {args.vcf}")
    audit_rows = []
    for i, cand in enumerate(candidates, 1):
        if i % 10 == 0:
            print(f"  [audit] {i}/{len(candidates)}")
        audit_rows.append(audit_candidate(args.vcf, cand, pca_carriers, args.bp_window))

    # Write outputs
    audit_path = os.path.join(args.outdir, "breakpoint_audit.tsv")
    keys_path = os.path.join(args.outdir, "registry_q7b_audit_keys.tsv")
    interp_path = os.path.join(args.outdir, "breakpoint_audit_interpretation.txt")

    print(f"[audit] Writing {audit_path}")
    write_audit_table(audit_rows, audit_path)

    print(f"[audit] Writing {keys_path}")
    write_registry_keys(audit_rows, keys_path)

    print(f"[audit] Writing {interp_path}")
    write_interpretation(audit_rows, interp_path)

    # Summary
    print(f"\n[audit] SUMMARY")
    print(f"  Candidates audited: {len(audit_rows)}")
    n_strong = sum(1 for r in audit_rows if r["pca_carrier_sv_support_pct"] >= 0.7)
    n_moderate = sum(1 for r in audit_rows if 0.5 <= r["pca_carrier_sv_support_pct"] < 0.7)
    n_weak = sum(1 for r in audit_rows if r["pca_carrier_sv_support_pct"] < 0.5)
    print(f"  Strong breakpoint support (>=70%):  {n_strong}")
    print(f"  Moderate support (50-70%):          {n_moderate}")
    print(f"  Weak support (<50%):                {n_weak}")
    print(f"[audit] Done.")


if __name__ == "__main__":
    main()
