#!/usr/bin/env python3
"""
STEP_02_bnd_sided_support.py — single-sided BND support (STEP 2 of 4)

=============================================================================
PIPELINE POSITION
=============================================================================
  STEP 1  assembled-junction forensics       (q4_mechanism)
→ STEP 2  single-sided BND support           (q7_existence_audit)  ← THIS SCRIPT
  STEP 3A cross-species synteny bridge       (cross_species)
  STEP 3B bp-pipeline bridge                 (bp_bridge)
  STEP 4  structural-class assignment        (phase_9)

Renamed from bnd_sided_support.py in pass 18 (2026-04-24). Dispatched
by run_evidence_biology.sh after STEP 1. Output block is
bnd_sided_support.json; STEP 4 loads it (currently for future tier-
upgrade logic — the keys are extracted but not yet rules-consumed in v7).

=============================================================================
ROLE
=============================================================================
STEP06 (phase 3) pairs DELLY CT=3to3/CT=5to5 BNDs into orphan rescue
candidates. Single-sided BNDs with no partner are currently abandoned.
This script asks the DIFFERENT question: for every EXISTING Phase 2/3
candidate, is there a CT=3to3 BND near its left boundary and/or a
CT=5to5 BND near its right boundary, regardless of whether pairing was
possible?

Why this matters:
  The MODULE_4E WIKI §10 Q&A explicitly acknowledges this gap:
    "Single-sided orphans accumulate in the BND catalog but contribute
     no rescue candidates. They could in principle be cross-referenced
     with the population-signal regions... no module currently does this."

  At ~5-9x coverage, one breakpoint of an inversion frequently falls in
  clean sequence (strong SV support) and the other falls in a repeat
  (no SV support). The pipeline currently treats this as "no SV support"
  rather than "half-sided SV support" — which is strictly stronger than
  no support.

=============================================================================
KEYS WRITTEN (Q7B extension)
=============================================================================
  q7b_bnd_left_support           cat: precise | imprecise | absent
  q7b_bnd_right_support          cat: precise | imprecise | absent
  q7b_bnd_sided_support_class    cat: both_sides | left_only |
                                       right_only | none
  q7b_bnd_left_ct_count          int
  q7b_bnd_right_ct_count         int
  q7b_bnd_left_best_qual         numeric
  q7b_bnd_right_best_qual        numeric

=============================================================================
INPUTS
=============================================================================
  --candidates      TSV: candidate_id, chrom, left_bp, right_bp,
                         left_zone_bp (optional), right_zone_bp (optional)
  --delly_bnd_vcf   DELLY BND catalog (bgzipped, MODULE_4E strict catalog)
  --manta_raw_vcf   Manta raw pre-conversion VCF (bgzipped, optional)
                    — used to cross-confirm sided support across callers
  --outdir          Where to write per-candidate JSON
  --bp_window       Default 10000 bp: how close to a boundary must a BND
                    fall to count as support for that side

=============================================================================
OUTPUT
=============================================================================
  <outdir>/<candidate_id>/structured/bnd_sided_support.json

If --registries_root is provided, writes a block of the same name via
the registry API.

=============================================================================
"""
import argparse
import csv
import gzip
import json
import os
import sys
from pathlib import Path


def open_vcf(path):
    if not path or not os.path.exists(path):
        return None
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def parse_info(info_str):
    out = {}
    for field in info_str.split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            out[k] = v
        else:
            out[field] = True
    return out


def classify_alt_bracket(alt):
    """
    Classify Manta BND ALT bracket pattern to DELLY CT equivalent.
      ]p]t → 3to3
      t[p[ → 5to5
      t]p] → 3to5
      [p[t → 5to3
    """
    if alt.startswith("]") and alt.endswith("]"):
        return "3to3"
    if alt.startswith("[") and alt.endswith("["):
        return "5to3"
    # Starts with a base
    if not alt:
        return None
    # `t[p[`
    if "[" in alt and alt.index("[") > 0:
        return "5to5"
    if "]" in alt and alt.index("]") > 0:
        return "3to5"
    return None


def parse_bnd_ct(vcf_path, caller="delly"):
    """
    Yield BND records with inversion orientation (CT=3to3 or CT=5to5).
    Skip inter-chromosomal (CHR2 != CHROM).
    For DELLY: uses INFO/CT.
    For Manta raw: parses ALT bracket pattern.
    """
    if not vcf_path:
        return []
    f = open_vcf(vcf_path)
    if f is None:
        return []
    records = []
    try:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom = parts[0]
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            alt = parts[4]
            qual_str = parts[5]
            try:
                qual = float(qual_str) if qual_str not in (".", "") else None
            except ValueError:
                qual = None
            filt = parts[6]
            info = parse_info(parts[7])

            svtype = info.get("SVTYPE", "")
            if svtype != "BND":
                continue
            chr2 = info.get("CHR2", chrom)
            if chr2 != chrom:
                continue  # skip inter-chromosomal

            if caller == "delly":
                ct = info.get("CT", "")
            else:
                ct = classify_alt_bracket(alt)

            if ct not in ("3to3", "5to5"):
                continue

            precise = "PRECISE" in info
            records.append({
                "caller": caller,
                "chrom": chrom,
                "pos": pos,
                "ct": ct,
                "qual": qual,
                "filter": filt,
                "precise": precise,
            })
    finally:
        f.close()
    return records


def count_bnds_near_boundary(records, chrom, boundary_bp, ct_required,
                              bp_window):
    """
    Count + describe BNDs with CT==ct_required within bp_window of boundary_bp.
    Returns (count, best_record_precision, best_qual).
    """
    matches = [r for r in records
               if r["chrom"] == chrom and r["ct"] == ct_required and
               abs(r["pos"] - boundary_bp) <= bp_window]
    if not matches:
        return 0, "absent", None

    any_precise = any(m["precise"] for m in matches)
    precision = "precise" if any_precise else "imprecise"
    quals = [m["qual"] for m in matches if m["qual"] is not None]
    best_qual = max(quals) if quals else None
    return len(matches), precision, best_qual


def sided_support_class(left_precision, right_precision):
    """Combine left and right precision classes into overall sided support."""
    if left_precision in ("precise", "imprecise") and \
       right_precision in ("precise", "imprecise"):
        return "both_sides"
    if left_precision in ("precise", "imprecise"):
        return "left_only"
    if right_precision in ("precise", "imprecise"):
        return "right_only"
    return "none"


def load_candidates(path):
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def try_registry(registries_root):
    if not registries_root:
        return None
    api_path = os.path.join(registries_root, "api", "python")
    if not os.path.isdir(api_path):
        return None
    if api_path not in sys.path:
        sys.path.insert(0, api_path)
    try:
        from registry_loader import load_registry
        return load_registry()
    except Exception as e:
        print(f"[bnd_sided] registry unavailable: {e}", file=sys.stderr)
        return None


def main():
    p = argparse.ArgumentParser(description=__doc__.split("=" * 77)[0])
    p.add_argument("--candidates", required=True)
    p.add_argument("--delly_bnd_vcf", default="")
    p.add_argument("--manta_raw_vcf", default="")
    p.add_argument("--outdir", required=True)
    p.add_argument("--registries_root", default="")
    p.add_argument("--bp_window", type=int, default=10000)
    p.add_argument("--dry_run", action="store_true")
    args = p.parse_args()

    cand_rows = load_candidates(args.candidates)
    print(f"[bnd_sided] {len(cand_rows)} candidates loaded")

    delly_bnds = parse_bnd_ct(args.delly_bnd_vcf, caller="delly") \
        if args.delly_bnd_vcf else []
    manta_bnds = parse_bnd_ct(args.manta_raw_vcf, caller="manta") \
        if args.manta_raw_vcf else []

    all_bnds = delly_bnds + manta_bnds
    print(f"[bnd_sided] delly inversion-orientation BNDs: {len(delly_bnds)}")
    print(f"[bnd_sided] manta inversion-orientation BNDs: {len(manta_bnds)}")

    reg = try_registry(args.registries_root)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    counts_class = {"both_sides": 0, "left_only": 0,
                    "right_only": 0, "none": 0}
    for row in cand_rows:
        try:
            chrom = row["chrom"]
            left = int(row["left_bp"])
            right = int(row["right_bp"])
        except (KeyError, ValueError):
            continue

        # Left boundary gets CT=3to3 evidence; right gets CT=5to5
        left_n, left_prec, left_q = count_bnds_near_boundary(
            all_bnds, chrom, left, "3to3", args.bp_window
        )
        right_n, right_prec, right_q = count_bnds_near_boundary(
            all_bnds, chrom, right, "5to5", args.bp_window
        )

        sc = sided_support_class(left_prec, right_prec)
        counts_class[sc] += 1

        data = {
            "q7b_bnd_left_support": left_prec,
            "q7b_bnd_right_support": right_prec,
            "q7b_bnd_sided_support_class": sc,
            "q7b_bnd_left_ct_count": left_n,
            "q7b_bnd_right_ct_count": right_n,
            "q7b_bnd_left_best_qual": left_q,
            "q7b_bnd_right_best_qual": right_q,
        }

        cid = row["candidate_id"]
        cand_out = outdir / cid / "structured"
        cand_out.mkdir(parents=True, exist_ok=True)
        block = {
            "block_type": "bnd_sided_support",
            "candidate_id": cid,
            "source_script": "STEP_02_bnd_sided_support.py",
            "data": data,
        }
        if not args.dry_run:
            (cand_out / "bnd_sided_support.json").write_text(
                json.dumps(block, indent=2, default=str)
            )
            if reg is not None:
                try:
                    reg.evidence.write_block(
                        candidate_id=cid,
                        block_type="bnd_sided_support",
                        data=data,
                        source_script="STEP_02_bnd_sided_support.py",
                    )
                except Exception as e:
                    print(f"[bnd_sided] {cid}: registry write failed: {e}",
                          file=sys.stderr)

    print(f"[bnd_sided] sided-support distribution: {counts_class}")


if __name__ == "__main__":
    main()
