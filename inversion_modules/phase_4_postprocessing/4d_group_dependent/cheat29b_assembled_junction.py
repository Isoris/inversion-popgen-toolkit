#!/usr/bin/env python3
"""
cheat29b_assembled_junction.py — assembled-junction forensics

=============================================================================
ROLE
=============================================================================
Extracts junction evidence from ASSEMBLED contigs (DELLY CONSENSUS, Manta
PRECISE ALT) rather than from reference-genome context at the estimated
boundary_bp. Complements cheat29 (which reads reference sequence): this
script reads the actual junction sequence DELLY/Manta produced from split
reads.

Why this is different from cheat29:

  cheat29 does:
    samtools faidx REF <chr>:<bp-50>-<bp>   # LEFT of estimated boundary
    samtools faidx REF <chr>:<bp+1>-<bp+50> # RIGHT of estimated boundary
    → detect microhomology at the reference junction
    → If boundary_bp is off by kilobases (common at 9x for NAHR-like events),
      the result describes genomic context, not the rearrangement junction.

  cheat29b does:
    parse DELLY BND records with CT=3to3/5to5 that fall in the candidate's
    boundary zone → read INFO/CONSENSUS (assembled contig spanning the
    breakpoint, from split reads)
    parse Manta INV records with PRECISE + HOMLEN/HOMSEQ → read HOMLEN,
    HOMSEQ directly
    → detect microhomology on the assembled junction
    → This IS the rearrangement junction.

The two complement each other:
- When a PRECISE record is available, cheat29b is the real answer.
- When no PRECISE record is available, cheat29 (reference fallback) is all
  we have, and the output is flagged with ref_ prefix.

=============================================================================
KEYS WRITTEN (Q4b_asm sub-block of mechanism.schema.json v5)
=============================================================================
  q4b_asm_precise_record_available    bool
  q4b_asm_source                       cat: delly_inv | delly_bnd_pair |
                                            delly_bnd_single | manta_inv |
                                            manta_bnd_pair | none
  q4b_asm_homlen                       int   from INFO/HOMLEN
  q4b_asm_homseq                       str   from INFO/HOMSEQ
  q4b_asm_consensus_available          bool  from INFO/CONSENSUS presence
  q4b_asm_consensus_length             int
  q4b_asm_junction_class               cat: blunt | MH_short | MH_long |
                                            insertion | unclear
  q4b_asm_vs_ref_concordance           cat: agree | disagree |
                                            ref_unavailable | asm_unavailable

=============================================================================
INPUTS
=============================================================================
  --candidates      TSV with candidate_id, chrom, left_bp, right_bp, left_zone_bp, right_zone_bp
  --delly_inv_vcf   DELLY INV catalog (bgzipped)
  --delly_bnd_vcf   DELLY BND catalog (bgzipped) — the MODULE_4E output
  --manta_inv_vcf   Manta INV catalog (bgzipped) — post-convertInversion
  --manta_raw_vcf   Manta raw pre-conversion VCF (bgzipped) — optional
  --ref_fallback_json  cheat29 output per-candidate (for concordance check)
                       — optional
  --outdir          Where to write per-candidate structured/mechanism.json updates
  --registries_root Path to registries/ (enables registry API writes)
  --bp_window       Default 1000 bp: how close a BND must be to a boundary
                    to count as this candidate's junction

=============================================================================
OUTPUT
=============================================================================
  <outdir>/<candidate_id>/structured/mechanism_assembled.json
    (written independently from cheat29's reference-based mechanism.json;
     the synthesis step in 4e reads both and picks a preference)

If --registries_root is provided, also writes via
reg.evidence.write_block(candidate_id, "mechanism_assembled", data).

=============================================================================
"""
import argparse
import csv
import gzip
import json
import os
import sys
from collections import defaultdict
from pathlib import Path


def open_vcf(path):
    if not path or not os.path.exists(path):
        return None
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def parse_info(info_str):
    """Parse VCF INFO field to dict."""
    out = {}
    for field in info_str.split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            out[k] = v
        else:
            out[field] = True
    return out


def classify_junction(homlen, has_consensus, ins_len=0):
    """
    Classify from HOMLEN + presence of insertion.

    This is the ASSEMBLED-junction classifier (operates on data DELLY/Manta
    produced from split reads, not from reference-genome context).
    """
    if homlen is None:
        return "unclear"
    try:
        homlen = int(homlen)
    except (ValueError, TypeError):
        return "unclear"
    if ins_len and int(ins_len) > 0:
        return "insertion"
    if homlen == 0:
        return "blunt"
    if 1 <= homlen <= 7:
        return "MH_short"
    if homlen >= 8:
        return "MH_long"
    return "unclear"


def parse_delly_inv(path):
    """
    Parse DELLY INV catalog. Return list of records with relevant fields.
    Each record: dict with chrom, pos, end, homlen, homseq, consensus, precise.
    """
    if not path:
        return []
    out = []
    f = open_vcf(path)
    if f is None:
        return []
    try:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            info = parse_info(parts[7])
            # Only INV, only PRECISE
            svtype = info.get("SVTYPE", "")
            if svtype != "INV":
                continue
            precise = "PRECISE" in info  # DELLY flag without =value
            if not precise:
                continue
            try:
                end_bp = int(info.get("END", 0))
            except ValueError:
                end_bp = 0
            rec = {
                "source": "delly_inv",
                "chrom": parts[0],
                "pos": int(parts[1]),
                "end": end_bp,
                "homlen": info.get("HOMLEN"),
                "homseq": info.get("HOMSEQ"),
                "consensus": info.get("CONSENSUS"),
                "ct": info.get("CT", ""),
                "precise": True,
            }
            out.append(rec)
    finally:
        f.close()
    return out


def parse_delly_bnd(path):
    """
    Parse DELLY BND catalog. Only CT=3to3 and CT=5to5 (inversion orientation).
    Each record carries its boundary position and (if PRECISE) CONSENSUS.
    """
    if not path:
        return []
    out = []
    f = open_vcf(path)
    if f is None:
        return []
    try:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            info = parse_info(parts[7])
            ct = info.get("CT", "")
            if ct not in ("3to3", "5to5"):
                continue  # skip DEL/DUP orientations
            if "PRECISE" not in info:
                continue  # only PRECISE BNDs are useful for junction forensics
            # chr2 for intra-chromosomal check
            chr2 = info.get("CHR2", parts[0])
            if chr2 != parts[0]:
                continue  # translocation, not inversion
            try:
                end_bp = int(info.get("END", 0))
            except ValueError:
                end_bp = 0
            rec = {
                "source": "delly_bnd",
                "chrom": parts[0],
                "pos": int(parts[1]),
                "end": end_bp,
                "ct": ct,
                "homlen": info.get("HOMLEN"),
                "homseq": info.get("HOMSEQ"),
                "consensus": info.get("CONSENSUS"),
                "cipos": info.get("CIPOS"),
                "ciend": info.get("CIEND"),
                "precise": True,
            }
            out.append(rec)
    finally:
        f.close()
    return out


def parse_manta_inv(path):
    """
    Parse Manta INV catalog. Records come from convertInversion (INV3+INV5
    pair merged). PRECISE = assembly succeeded. HOMLEN/HOMSEQ available.
    """
    if not path:
        return []
    out = []
    f = open_vcf(path)
    if f is None:
        return []
    try:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            info = parse_info(parts[7])
            svtype = info.get("SVTYPE", "")
            if svtype != "INV":
                continue
            # Manta: PRECISE means full ALT assembled; IMPRECISE means
            # flanking-only. We want PRECISE.
            if "IMPRECISE" in info:
                continue
            try:
                end_bp = int(info.get("END", 0))
            except ValueError:
                end_bp = 0
            rec = {
                "source": "manta_inv",
                "chrom": parts[0],
                "pos": int(parts[1]),
                "end": end_bp,
                "homlen": info.get("HOMLEN"),
                "homseq": info.get("HOMSEQ"),
                # Manta doesn't emit CONSENSUS but emits the symbolic ALT
                # <INV>; the assembly is implicit in PRECISE + HOMLEN
                "consensus": "<INV>" if parts[4] == "<INV>" else None,
                "cipos": info.get("CIPOS"),
                "ciend": info.get("CIEND"),
                "precise": True,
            }
            out.append(rec)
    finally:
        f.close()
    return out


def match_record_to_candidate(record, cand_row, bp_window):
    """
    Return True if this SV record's boundaries fall within bp_window of the
    candidate's boundary zones.

    For INV records: pos ~ left_bp AND end ~ right_bp.
    For BND records: pos ~ left_bp OR pos ~ right_bp (depending on CT).
    """
    try:
        left = int(cand_row["left_bp"])
        right = int(cand_row["right_bp"])
    except (KeyError, ValueError):
        return False

    if record["source"] in ("delly_inv", "manta_inv"):
        # INV record spans the candidate
        near_left = abs(record["pos"] - left) <= bp_window
        near_right = abs(record["end"] - right) <= bp_window
        return near_left and near_right

    if record["source"] == "delly_bnd":
        # BND is a single junction
        if record["ct"] == "3to3":
            # left junction orientation
            return abs(record["pos"] - left) <= bp_window
        if record["ct"] == "5to5":
            return abs(record["pos"] - right) <= bp_window

    return False


def select_best_record(matching):
    """
    Preference order for the best junction-forensics source:
    1. delly_inv (PRECISE) — has CONSENSUS + HOMLEN
    2. manta_inv (PRECISE) — has HOMLEN, ALT is <INV> symbolic
    3. delly_bnd pair (both 3to3 + 5to5 present for the same candidate)
       — two BND records can each carry CONSENSUS
    4. delly_bnd single-sided — still has CONSENSUS for one boundary
    """
    if not matching:
        return None, "none"

    # 1. Any delly_inv?
    delly_inv = [r for r in matching if r["source"] == "delly_inv"]
    if delly_inv:
        return delly_inv[0], "delly_inv"

    # 2. Any manta_inv?
    manta_inv = [r for r in matching if r["source"] == "manta_inv"]
    if manta_inv:
        return manta_inv[0], "manta_inv"

    # 3/4. DELLY BNDs
    bnds = [r for r in matching if r["source"] == "delly_bnd"]
    left_bnds = [r for r in bnds if r["ct"] == "3to3"]
    right_bnds = [r for r in bnds if r["ct"] == "5to5"]
    if left_bnds and right_bnds:
        return left_bnds[0], "delly_bnd_pair"  # prefer left as primary
    if left_bnds:
        return left_bnds[0], "delly_bnd_single"
    if right_bnds:
        return right_bnds[0], "delly_bnd_single"
    return None, "none"


def build_mechanism_asm_block(cand_row, delly_inv_recs, delly_bnd_recs,
                               manta_inv_recs, bp_window, ref_fallback=None):
    """
    For one candidate, produce the data dict for the mechanism_assembled block.
    """
    matching = []
    for rec in delly_inv_recs + delly_bnd_recs + manta_inv_recs:
        if rec["chrom"] != cand_row["chrom"]:
            continue
        if match_record_to_candidate(rec, cand_row, bp_window):
            matching.append(rec)

    best, src = select_best_record(matching)

    if best is None:
        data = {
            "q4b_asm_precise_record_available": False,
            "q4b_asm_source": "none",
            "q4b_asm_consensus_available": False,
            "q4b_asm_junction_class": "unresolved_no_precise_sv",
            "q4b_asm_vs_ref_concordance": "asm_unavailable",
            "n_matching_records": 0,
        }
    else:
        homlen = best.get("homlen")
        has_consensus = best.get("consensus") not in (None, "", "<INV>")
        jc = classify_junction(homlen, has_consensus)
        data = {
            "q4b_asm_precise_record_available": True,
            "q4b_asm_source": src,
            "q4b_asm_homlen": int(homlen) if homlen is not None else None,
            "q4b_asm_homseq": best.get("homseq"),
            "q4b_asm_consensus_available": has_consensus,
            "q4b_asm_consensus_length": (
                len(best["consensus"]) if has_consensus else None
            ),
            "q4b_asm_junction_class": jc,
            "n_matching_records": len(matching),
        }

        # Concordance with cheat29 reference-based output
        if ref_fallback and cand_row["candidate_id"] in ref_fallback:
            ref = ref_fallback[cand_row["candidate_id"]]
            ref_class = ref.get("junction_class") or ref.get("ref_junction_class")
            if ref_class:
                data["q4b_ref_junction_class"] = ref_class
                if ref_class == jc:
                    data["q4b_asm_vs_ref_concordance"] = "agree"
                else:
                    data["q4b_asm_vs_ref_concordance"] = "disagree"
            else:
                data["q4b_asm_vs_ref_concordance"] = "ref_unavailable"
        else:
            data["q4b_asm_vs_ref_concordance"] = "ref_unavailable"

    return data


def load_candidates(path):
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def load_ref_fallback(path):
    if not path or not os.path.exists(path):
        return {}
    with open(path) as f:
        return json.load(f)


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
        print(f"[cheat29b] registry unavailable: {e}", file=sys.stderr)
        return None


def main():
    p = argparse.ArgumentParser(description=__doc__.split("=" * 77)[0])
    p.add_argument("--candidates", required=True,
                   help="TSV with candidate_id, chrom, left_bp, right_bp")
    p.add_argument("--delly_inv_vcf", default="")
    p.add_argument("--delly_bnd_vcf", default="")
    p.add_argument("--manta_inv_vcf", default="")
    p.add_argument("--manta_raw_vcf", default="",
                   help="(reserved for future use; currently unused)")
    p.add_argument("--ref_fallback_json", default="",
                   help="cheat29 output JSON for concordance check")
    p.add_argument("--outdir", required=True)
    p.add_argument("--registries_root", default="")
    p.add_argument("--bp_window", type=int, default=1000)
    p.add_argument("--dry_run", action="store_true")
    args = p.parse_args()

    cand_rows = load_candidates(args.candidates)
    print(f"[cheat29b] {len(cand_rows)} candidates loaded")

    print("[cheat29b] parsing SV catalogs …")
    delly_inv = parse_delly_inv(args.delly_inv_vcf) if args.delly_inv_vcf else []
    delly_bnd = parse_delly_bnd(args.delly_bnd_vcf) if args.delly_bnd_vcf else []
    manta_inv = parse_manta_inv(args.manta_inv_vcf) if args.manta_inv_vcf else []

    print(f"[cheat29b]   delly INV PRECISE:  {len(delly_inv)}")
    print(f"[cheat29b]   delly BND PRECISE:  {len(delly_bnd)}")
    print(f"[cheat29b]   manta INV PRECISE:  {len(manta_inv)}")

    ref_fallback = load_ref_fallback(args.ref_fallback_json)

    reg = try_registry(args.registries_root)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    n_precise = 0
    for row in cand_rows:
        data = build_mechanism_asm_block(
            row, delly_inv, delly_bnd, manta_inv,
            args.bp_window, ref_fallback
        )
        if data["q4b_asm_precise_record_available"]:
            n_precise += 1

        cid = row["candidate_id"]
        cand_out = outdir / cid / "structured"
        cand_out.mkdir(parents=True, exist_ok=True)
        block = {
            "block_type": "mechanism_assembled",
            "candidate_id": cid,
            "source_script": "cheat29b_assembled_junction.py",
            "data": data,
        }
        if not args.dry_run:
            (cand_out / "mechanism_assembled.json").write_text(
                json.dumps(block, indent=2, default=str)
            )
            if reg is not None:
                try:
                    reg.evidence.write_block(
                        candidate_id=cid,
                        block_type="mechanism_assembled",
                        data=data,
                        source_script="cheat29b_assembled_junction.py",
                    )
                except Exception as e:
                    print(f"[cheat29b] {cid}: registry write failed: {e}",
                          file=sys.stderr)

    print(f"[cheat29b] {n_precise}/{len(cand_rows)} candidates have a "
          f"PRECISE assembled junction record")
    print(f"[cheat29b] {len(cand_rows) - n_precise} remain on reference fallback")


if __name__ == "__main__":
    main()
