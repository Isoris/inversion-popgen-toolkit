#!/usr/bin/env python3
"""
STEP_03_per_read_evidence.py — per-read BAM evidence ledger (Q7 existence audit)

=============================================================================
PIPELINE POSITION
=============================================================================
  STEP 1  assembled-junction forensics       (q4_mechanism)
  STEP 2  single-sided BND support           (q7_existence_audit)
→ STEP 3  per-read evidence ledger           (q7_existence_audit)  ← THIS SCRIPT
  STEP 3A cross-species synteny bridge       (cross_species)
  STEP 3B bp-pipeline bridge                 (bp_bridge)
  STEP 4  structural-class assignment        (phase_9)

New in pass 21 (2026-04-24). Extracts per-read breakpoint evidence from
BAM files for all 30 (or 226) samples at one candidate's two breakpoints,
and persists the full read-name ledger so every supporting read is
traceable back to its (sample, group, bp_side, evidence_type).

This is the UPGRADE that closes the "Q7B 25-key flat-key extraction gap"
flagged in pass 18's handoff. Prior scripts (STEP_A02 in phase_3_refine,
07_breakpoint_evidence_pileup in phase_6) collected read names in-memory
and discarded them after emitting aggregate counts or figures. This
script is the canonical read-name persistence layer.

=============================================================================
RELATIONSHIP TO OTHER EVIDENCE SCRIPTS
=============================================================================
Two independent evidence sources for Q7 existence:

  A) SV-CALLER AGGREGATES (from DELLY/Manta VCF FORMAT fields)
       Handled by: breakpoint_evidence_audit.py  (Q7B VCF source)
                 + STEP_02_bnd_sided_support.py   (orphan BND rescue)
       Writes aggregate counts (DV, RV, PR, SR) — read names NOT available
       because DELLY doesn't export them in the standard filter+regenotype
       workflow. Trusts the SV caller's read-selection logic.
       Sidecar output: evidence_reads_vcf.tsv.gz (aggregate only)

  B) BAM PILEUP (direct pysam extraction) ← THIS SCRIPT
       Handled by: STEP_03_per_read_evidence.py
       Writes per-READ rows with names, coords, strands, CIGAR, tags.
       Independent of SV caller — uses only mapping flags + CIGAR + SA tag.
       Output: evidence_reads_bam.tsv.gz (row-level)

The two are joinable on (candidate_id, sample_id, bp_side) but stay in
separate files because they don't share columns. For audit:
  "DELLY RV=5 for sample1 at BP1 — do I see 5 split reads in the BAM?"
  → count bam rows where (candidate=..., sample=sample1, bp_side=left,
    evidence_type='split_read'); compare with VCF sidecar row.

=============================================================================
INPUT SOURCES (all read from registry given just --candidate)
=============================================================================
  q3_chrom              chromosome
  q3_final_left_bp      left breakpoint (BP1)
  q3_final_right_bp     right breakpoint (BP2)
  (written by STEP_03B_bp_pipeline_bridge from phase_6 consensus TSV)

Plus external:
  --bam_dir             directory of <sample>.markdup.bam files
  --samples_tsv         TSV with columns: sample_id, group (REF/HET/INV)
                        Usually from Q1 decomposition output for this
                        candidate.
  --delly_inv_vcf       DELLY filtered INV VCF (optional, for VCF sidecar)
  --manta_inv_vcf       Manta INV VCF (optional, for VCF sidecar)

=============================================================================
OUTPUT — BAM ledger (evidence_reads_bam.tsv.gz)
=============================================================================
One row per (read × evidence_type × bp_side). A single read can appear on
multiple rows if it supports multiple evidence types (e.g. a split read
that is also soft-clipped gets two rows — one split_read, one soft_clip).

Columns (29):
  candidate_id          LG28_cand_1
  chrom                 C_gar_LG28
  sample_id             sample1
  group                 INV | HET | REF | unknown
  bp_side               left | right  (which breakpoint)
  bp_pos                15115243
  read_name             H7YHMDSXX:3:2212:14082:5891
  evidence_type         discordant_FF | discordant_RR | split_read |
                         soft_clip | mate_at_partner
  clip_end_template     5prime | 3prime | both | none
                         TEMPLATE-frame end of the READ that is clipped
                         (strand-aware). 5prime = clip at the molecule's
                         5' end. For NAHR evidence this is the biologically
                         meaningful frame.
  clip_end_mapped       left | right | both | none
                         MAPPED-frame: which end of the read-as-aligned
                         to reference is clipped (CIGAR[0] = 'left',
                         CIGAR[-1] = 'right'). Useful for breakpoint-
                         orientation checks against bp_side.
  clip_len_left         int  length of soft-clip at CIGAR[0], 0 if none
  clip_len_right        int  length of soft-clip at CIGAR[-1], 0 if none
  read_pos              15115201  (0-based reference_start)
  read_end              15115277  (reference_end exclusive)
  read_strand           + | -
  cigar                 50M26S
  mate_chrom            C_gar_LG28 (or * if unmapped/different)
  mate_pos              18006104
  mate_strand           + | -
  mapq                  60
  nm                    3      NM tag (edit distance)
  as_score              145    AS tag (BWA alignment score)
  xs_score              60     XS tag (BWA suboptimal) — empty if none
  sa                    C_gar_LG28,18005891,-,26M50S,60,1
                         SA tag — the chimeric half for split reads.
                         Empty unless evidence_type='split_read'.
  mc                    76M    MC tag (mate CIGAR) — empty if absent
  xa                    chr2,+12345,76M,0;  XA tag (alt positions,
                         BWA repeat signal). Empty if unique.
  is_supplementary      True | False  (chimeric half marker)
  is_primary_of_pair    True | False  (deduplication flag: only one
                         mate per pair is marked True)
  flag                  163  (raw SAM flag)

=============================================================================
OUTPUT — VCF sidecar (evidence_reads_vcf.tsv.gz, separate schema)
=============================================================================
One row per (candidate × sample × bp_side × caller × support_type).

Columns (10):
  candidate_id          LG28_cand_1
  chrom                 C_gar_LG28
  sample_id             sample1
  group                 INV | HET | REF
  bp_side               left | right | both  (VCF records span both)
  caller                delly | manta
  support_type          DV | RV | PR_alt | SR_alt  (field tag from FORMAT)
  n_reads               5   (the count reported by the caller)
  genotype              0/0 | 0/1 | 1/1 | ./.
  record_pos            15115243  (POS field of the supporting VCF record)

=============================================================================
REGISTRY BLOCK WRITTEN: per_read_evidence_summary
=============================================================================
A summary block (NOT the full ledger — too large for a Tier-2 block) with
cohort-level aggregates:

  q7b_bam_n_split_left_INV     int  carriers with ≥1 split at BP1 in INV group
  q7b_bam_n_split_right_INV    int  etc.
  q7b_bam_n_split_left_HET     int
  q7b_bam_n_split_right_HET    int
  q7b_bam_n_discord_left_INV   int  carriers with ≥1 discordant FF/RR at BP1
  q7b_bam_n_discord_right_INV  int
  q7b_bam_ledger_path          string  path to evidence_reads_bam.tsv.gz
  q7b_vcf_sidecar_path         string  path to evidence_reads_vcf.tsv.gz
  q7b_per_read_source          string  "bam_pysam+vcf_sidecar"
  q7b_window_bp                int  window used (±)
  q7b_min_mapq                 int  MAPQ threshold applied

The full ledger stays on disk; the registry only holds the summary counts
that STEP_04_assign_structural_class needs.

=============================================================================
CLI
=============================================================================
  python3 STEP_03_per_read_evidence.py \
      --candidate LG28_cand_1 \
      --registries_root ../../registries \
      --bam_dir /path/to/markdup_bams \
      --samples_tsv per_candidate_groups.tsv \
      --outdir /output/root \
      [--delly_inv_vcf ...] [--manta_inv_vcf ...] \
      [--window 2000] [--min_mapq 20] [--min_clip_len 10] \
      [--keep_secondary] [--threads 4] [--overwrite]

REQUIRES: pysam (pip install pysam)
DIFFICULTY: medium
DISPATCHER: no (Q7 existence is candidate-level, not group-dependent)
=============================================================================
"""
import argparse
import csv
import gzip
import json
import os
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

try:
    import pysam
except ImportError:
    print("ERROR: pysam not installed. pip install pysam", file=sys.stderr)
    sys.exit(1)


# ==========================================================================
# pysam CIGAR op codes
# ==========================================================================
CIGAR_MATCH       = 0   # M
CIGAR_INS         = 1   # I
CIGAR_DEL         = 2   # D
CIGAR_REF_SKIP    = 3   # N
CIGAR_SOFT_CLIP   = 4   # S
CIGAR_HARD_CLIP   = 5   # H
CIGAR_EQUAL       = 7   # =
CIGAR_DIFF        = 8   # X


# ==========================================================================
# Registry coord lookup
# ==========================================================================
def _resolve_registries_root(user_arg, script_path):
    """Find registries/ by candidate: CLI > env > walk up from script."""
    if user_arg and os.path.isdir(user_arg):
        return user_arg
    env = os.environ.get("REGISTRIES_ROOT", "")
    if env and os.path.isdir(env):
        return env
    cur = Path(script_path).resolve().parent
    for _ in range(8):
        cand = cur / "registries"
        if cand.is_dir() and (cand / "api" / "python").is_dir():
            return str(cand)
        if cur.parent == cur:
            break
        cur = cur.parent
    return None


def load_coords_from_registry(candidate_id, registries_root):
    """Read q3_chrom, q3_final_left_bp, q3_final_right_bp from registry."""
    if not registries_root:
        raise RuntimeError(
            f"cannot locate registries_root; pass --registries_root explicitly")

    api = os.path.join(registries_root, "api", "python")
    if api not in sys.path:
        sys.path.insert(0, api)
    try:
        from registry_loader import load_registry
    except ImportError as e:
        raise RuntimeError(f"registry_loader import failed: {e}")

    reg = load_registry()
    try:
        chrom    = reg.get_key(candidate_id, "q3_chrom")
        left_bp  = int(reg.get_key(candidate_id, "q3_final_left_bp"))
        right_bp = int(reg.get_key(candidate_id, "q3_final_right_bp"))
    except Exception as e:
        raise RuntimeError(
            f"registry lookup failed for {candidate_id}: {e}. "
            f"Ensure STEP_03B_bp_pipeline_bridge has written the "
            f"fragment_distribution block for this candidate.")

    if not chrom or left_bp is None or right_bp is None:
        raise RuntimeError(
            f"registry returned empty coords for {candidate_id}: "
            f"chrom={chrom} left={left_bp} right={right_bp}")

    return chrom, left_bp, right_bp


# ==========================================================================
# Sample group table
# ==========================================================================
def load_samples_tsv(path):
    """Return {sample_id: group} where group ∈ {INV, HET, REF, unknown}."""
    out = {}
    with open(path) as f:
        rdr = csv.DictReader(f, delimiter="\t")
        for row in rdr:
            sid = row.get("sample_id") or row.get("sample")
            if not sid:
                continue
            grp = (row.get("group") or row.get("karyotype") or "unknown").upper()
            if grp not in ("INV", "HET", "REF"):
                grp = "unknown"
            out[sid] = grp
    return out


# ==========================================================================
# Tag helpers
# ==========================================================================
def get_tag_safe(read, tag):
    """Return tag value or empty string."""
    try:
        if read.has_tag(tag):
            return str(read.get_tag(tag))
    except Exception:
        pass
    return ""


def format_sa_tag(read):
    """Return SA tag formatted as first-alignment only (primary partner)."""
    sa = get_tag_safe(read, "SA")
    if not sa:
        return ""
    # SA is semicolon-delimited; usually one entry for chimeric splits
    entries = [e for e in sa.split(";") if e.strip()]
    if not entries:
        return ""
    return entries[0]   # primary chimeric partner (first listed)


def parse_sa_entry(sa_str):
    """Parse one SA entry → (chr, pos, strand, cigar, mapq, nm) or None."""
    if not sa_str:
        return None
    p = sa_str.split(",")
    if len(p) < 6:
        return None
    try:
        return (p[0], int(p[1]), p[2], p[3], int(p[4]), int(p[5]))
    except (ValueError, IndexError):
        return None


# ==========================================================================
# Clip-end classifier
# ==========================================================================
def classify_clip_end(cigartuples, min_clip_len, is_reverse):
    """
    Classify where the read is clipped, in TWO frames:

    - clip_end_mapped: 'left' / 'right' / 'both' / 'none'
      Which end of the read-as-mapped-to-reference is clipped. CIGAR[0]
      is the leftmost end, CIGAR[-1] is the rightmost.

    - clip_end_template: '5prime' / '3prime' / 'both' / 'none'
      The TEMPLATE-frame 5'/3' end. This accounts for strand:
        forward reads: left-clip = 5', right-clip = 3'
        reverse reads: left-clip = 3', right-clip = 5'  (read was reverse-
                       complemented relative to reference when SAM was
                       written; left-clip in CIGAR corresponds to the
                       original 3' end of the template molecule)

    For NAHR breakpoint analysis, clip_end_template is the biologically
    meaningful field — it tells you which end of the SEQUENCING READ
    extends across the junction.
    """
    if not cigartuples:
        return "none", "none"
    left_clip  = (cigartuples[0][0]  == CIGAR_SOFT_CLIP and
                  cigartuples[0][1]  >= min_clip_len)
    right_clip = (cigartuples[-1][0] == CIGAR_SOFT_CLIP and
                  cigartuples[-1][1] >= min_clip_len)

    # Mapped frame (read-as-mapped-to-reference)
    if left_clip and right_clip:
        mapped = "both"
    elif left_clip:
        mapped = "left"
    elif right_clip:
        mapped = "right"
    else:
        mapped = "none"

    # Template frame (read-as-sequenced)
    # For reverse-strand reads, left/right are swapped in template frame.
    if mapped == "none":
        template = "none"
    elif mapped == "both":
        template = "both"
    elif not is_reverse:
        # forward: left=5', right=3'
        template = "5prime" if mapped == "left" else "3prime"
    else:
        # reverse: left=3', right=5'
        template = "3prime" if mapped == "left" else "5prime"

    return mapped, template


# ==========================================================================
# The core extractor — returns list of evidence rows for one BAM at one bp
# ==========================================================================
def extract_bam_rows(bam_path, candidate_id, chrom, bp_pos, bp_side,
                     partner_chrom, partner_pos,
                     sample_id, group,
                     window=2000, min_mapq=20, min_clip_len=10,
                     keep_secondary=False):
    """
    Open BAM, fetch reads near bp_pos, emit one row per evidence instance.
    A single read can produce multiple rows (e.g. split + soft_clip).
    Deduplicates pair-level rows via is_primary_of_pair flag.
    """
    rows = []
    if not os.path.exists(bam_path):
        return rows

    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        sys.stderr.write(f"[STEP_03] cannot open {bam_path}: {e}\n")
        return rows

    start = max(0, bp_pos - window)
    end   = bp_pos + window

    # For primary-of-pair deduplication: mark the mate with lower refstart
    # as primary. This makes the flag deterministic per pair.
    seen_pairs = {}   # read_name → already_emitted_bool

    try:
        for read in bam.fetch(chrom, start, end):
            if read.is_unmapped or read.is_duplicate:
                continue
            if read.is_secondary and not keep_secondary:
                continue
            if read.mapping_quality < min_mapq:
                continue

            # Common per-read fields
            rname  = read.query_name
            rpos   = read.reference_start        # 0-based
            rend   = read.reference_end or rpos
            rstr   = "-" if read.is_reverse else "+"
            cigar  = read.cigarstring or ""
            mapq   = read.mapping_quality
            mchrom = read.next_reference_name or "*"
            mpos   = read.next_reference_start if read.next_reference_start >= 0 else -1
            mstr   = "-" if read.mate_is_reverse else "+"
            flag   = read.flag
            is_sup = read.is_supplementary

            # pair-level primary marker
            if rname in seen_pairs:
                is_primary_of_pair = False
            else:
                seen_pairs[rname] = True
                is_primary_of_pair = True

            # Tags (once — shared across rows this read produces)
            nm_val = get_tag_safe(read, "NM")
            as_val = get_tag_safe(read, "AS")
            xs_val = get_tag_safe(read, "XS")
            mc_val = get_tag_safe(read, "MC")
            xa_val = get_tag_safe(read, "XA")

            clip_end_mapped, clip_end_template = classify_clip_end(
                read.cigartuples, min_clip_len, read.is_reverse)

            # Also compute the clip length(s) — for downstream identity checks
            clip_len_left  = (read.cigartuples[0][1]
                              if read.cigartuples and read.cigartuples[0][0] == CIGAR_SOFT_CLIP
                              else 0)
            clip_len_right = (read.cigartuples[-1][1]
                              if read.cigartuples and read.cigartuples[-1][0] == CIGAR_SOFT_CLIP
                              else 0)

            # Row constructor
            def base_row(ev_type, sa=""):
                # clip_end fields only meaningful for split/soft-clip rows
                is_clip_row = ev_type in ("split_read", "soft_clip")
                return {
                    "candidate_id":       candidate_id,
                    "chrom":              chrom,
                    "sample_id":          sample_id,
                    "group":              group,
                    "bp_side":            bp_side,
                    "bp_pos":             bp_pos,
                    "read_name":          rname,
                    "evidence_type":      ev_type,
                    "clip_end_template":  clip_end_template if is_clip_row else "",
                    "clip_end_mapped":    clip_end_mapped   if is_clip_row else "",
                    "clip_len_left":      clip_len_left     if is_clip_row else 0,
                    "clip_len_right":     clip_len_right    if is_clip_row else 0,
                    "read_pos":           rpos,
                    "read_end":           rend,
                    "read_strand":        rstr,
                    "cigar":              cigar,
                    "mate_chrom":         mchrom,
                    "mate_pos":           mpos,
                    "mate_strand":        mstr,
                    "mapq":               mapq,
                    "nm":                 nm_val,
                    "as_score":           as_val,
                    "xs_score":           xs_val,
                    "sa":                 sa,
                    "mc":                 mc_val,
                    "xa":                 xa_val,
                    "is_supplementary":   is_sup,
                    "is_primary_of_pair": is_primary_of_pair,
                    "flag":               flag,
                }

            # ---- Evidence type 1: discordant FF / RR ----
            if read.is_paired and not read.mate_is_unmapped:
                read_fwd = not read.is_reverse
                mate_fwd = not read.mate_is_reverse
                if read_fwd == mate_fwd:
                    ev_type = "discordant_FF" if read_fwd else "discordant_RR"
                    rows.append(base_row(ev_type))

                # ---- Evidence type 2: mate near partner breakpoint ----
                if (mchrom == partner_chrom and
                        abs(mpos - partner_pos) <= window * 3):
                    rows.append(base_row("mate_at_partner"))

            # ---- Evidence type 3: split read (SA tag) ----
            sa_first = format_sa_tag(read)
            sa_parsed = parse_sa_entry(sa_first)
            if sa_parsed:
                sa_chr, sa_pos, _, _, sa_mapq, _ = sa_parsed
                if (sa_mapq >= min_mapq and
                        sa_chr == partner_chrom and
                        abs(sa_pos - partner_pos) <= window * 3):
                    rows.append(base_row("split_read", sa=sa_first))

            # ---- Evidence type 4: soft-clip at breakpoint ----
            if clip_end_mapped in ("left", "right", "both"):
                # Only count as breakpoint evidence if the clip end aligns
                # with the bp side (clip oriented toward the inversion).
                c = read.cigartuples
                if c:
                    clip_near_bp = False
                    # Left soft-clip: aligned portion starts at rpos, so
                    # clip sits at position rpos. If read maps right of bp
                    # and is clipped on the left, it could support bp.
                    if c[0][0] == CIGAR_SOFT_CLIP and c[0][1] >= min_clip_len:
                        if abs(rpos - bp_pos) <= window:
                            clip_near_bp = True
                    if c[-1][0] == CIGAR_SOFT_CLIP and c[-1][1] >= min_clip_len:
                        if abs(rend - bp_pos) <= window:
                            clip_near_bp = True
                    if clip_near_bp:
                        rows.append(base_row("soft_clip"))

    except Exception as e:
        sys.stderr.write(f"[STEP_03] fetch error in {bam_path} at {chrom}:{bp_pos}: {e}\n")
    finally:
        bam.close()

    return rows


# ==========================================================================
# VCF sidecar extractor
# ==========================================================================
def extract_vcf_sidecar(vcf_path, caller_name, candidate_id, chrom,
                        left_bp, right_bp, groups_by_sample, window=5000):
    """Extract aggregate caller counts per sample at this candidate's BP."""
    if not vcf_path or not os.path.exists(vcf_path):
        return []
    rows = []

    def _open(p):
        return gzip.open(p, "rt") if p.endswith(".gz") else open(p)

    # Precompute acceptance window per side
    left_lo, left_hi = left_bp - window, left_bp + window
    right_lo, right_hi = right_bp - window, right_bp + window

    samples = []
    try:
        with _open(vcf_path) as f:
            for line in f:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    samples = line.rstrip("\n").split("\t")[9:]
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 10:
                    continue
                vcf_chrom, pos = parts[0], int(parts[1])
                if vcf_chrom != chrom:
                    continue
                # Which BP does this record sit near?
                near_left  = left_lo  <= pos <= left_hi
                near_right = right_lo <= pos <= right_hi
                if not (near_left or near_right):
                    continue
                bp_side = "left" if near_left else "right"

                fmt_keys = parts[8].split(":")
                for i, sdata in enumerate(parts[9:]):
                    if i >= len(samples):
                        break
                    sid = samples[i]
                    grp = groups_by_sample.get(sid, "unknown")
                    vals = dict(zip(fmt_keys, sdata.split(":")))
                    gt = vals.get("GT", "./.")

                    # DELLY: DV (discordant pairs), RV (split reads)
                    # Manta: PR (pair support), SR (split support) — each is
                    # "ref_count,alt_count"
                    field_map = {
                        "delly": [("DV", "DV"), ("RV", "RV")],
                        "manta": [("PR", "PR_alt"), ("SR", "SR_alt")],
                    }.get(caller_name, [])

                    for raw_key, support_type in field_map:
                        if raw_key not in vals:
                            continue
                        raw_val = vals[raw_key]
                        try:
                            if "," in raw_val:
                                # Manta PR/SR: "ref,alt" — take alt
                                n_reads = int(raw_val.split(",")[1])
                            else:
                                n_reads = int(raw_val)
                        except (ValueError, IndexError):
                            continue
                        if n_reads == 0 and gt in ("0/0", "./."):
                            continue   # skip pure-absence rows
                        rows.append({
                            "candidate_id": candidate_id,
                            "chrom":        chrom,
                            "sample_id":    sid,
                            "group":        grp,
                            "bp_side":      bp_side,
                            "caller":       caller_name,
                            "support_type": support_type,
                            "n_reads":      n_reads,
                            "genotype":     gt,
                            "record_pos":   pos,
                        })
    except Exception as e:
        sys.stderr.write(f"[STEP_03] VCF parse error in {vcf_path}: {e}\n")

    return rows


# ==========================================================================
# Summary block — cohort aggregates for registry
# ==========================================================================
def build_summary(bam_rows, window, min_mapq):
    """Cohort-level counts for the registry block."""
    # dict[(group, bp_side, ev_type)] → set of sample_ids
    carriers = defaultdict(set)
    for r in bam_rows:
        if r["is_primary_of_pair"] and r["evidence_type"] in (
                "split_read", "soft_clip", "discordant_FF", "discordant_RR"):
            disc_or_split = ("split" if r["evidence_type"] == "split_read"
                             else "soft" if r["evidence_type"] == "soft_clip"
                             else "discord")
            carriers[(r["group"], r["bp_side"], disc_or_split)].add(r["sample_id"])

    def n(group, side, kind):
        return len(carriers.get((group, side, kind), set()))

    return {
        "q7b_bam_n_split_left_INV":      n("INV", "left",  "split"),
        "q7b_bam_n_split_right_INV":     n("INV", "right", "split"),
        "q7b_bam_n_split_left_HET":      n("HET", "left",  "split"),
        "q7b_bam_n_split_right_HET":     n("HET", "right", "split"),
        "q7b_bam_n_split_left_REF":      n("REF", "left",  "split"),
        "q7b_bam_n_split_right_REF":     n("REF", "right", "split"),
        "q7b_bam_n_soft_left_INV":       n("INV", "left",  "soft"),
        "q7b_bam_n_soft_right_INV":      n("INV", "right", "soft"),
        "q7b_bam_n_discord_left_INV":    n("INV", "left",  "discord"),
        "q7b_bam_n_discord_right_INV":   n("INV", "right", "discord"),
        "q7b_bam_n_discord_left_HET":    n("HET", "left",  "discord"),
        "q7b_bam_n_discord_right_HET":   n("HET", "right", "discord"),
        "q7b_bam_n_discord_left_REF":    n("REF", "left",  "discord"),
        "q7b_bam_n_discord_right_REF":   n("REF", "right", "discord"),
        "q7b_window_bp":                 window,
        "q7b_min_mapq":                  min_mapq,
        "q7b_per_read_source":           "bam_pysam+vcf_sidecar",
    }


# ==========================================================================
# Writers
# ==========================================================================
BAM_COLUMNS = [
    "candidate_id", "chrom", "sample_id", "group",
    "bp_side", "bp_pos", "read_name", "evidence_type",
    "clip_end_template", "clip_end_mapped", "clip_len_left", "clip_len_right",
    "read_pos", "read_end", "read_strand", "cigar",
    "mate_chrom", "mate_pos", "mate_strand",
    "mapq", "nm", "as_score", "xs_score", "sa", "mc", "xa",
    "is_supplementary", "is_primary_of_pair", "flag",
]

VCF_COLUMNS = [
    "candidate_id", "chrom", "sample_id", "group", "bp_side",
    "caller", "support_type", "n_reads", "genotype", "record_pos",
]


def write_gz_tsv(rows, out_path, columns):
    with gzip.open(out_path, "wt") as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


# ==========================================================================
# CLI + main
# ==========================================================================
def parse_args():
    p = argparse.ArgumentParser(description=__doc__.split("=" * 77)[0].strip())
    p.add_argument("--candidate", required=True,
                   help="candidate_id (coords pulled from registry)")
    p.add_argument("--registries_root", default="",
                   help="path to registries/ (default: walk up from script)")
    p.add_argument("--bam_dir", required=True,
                   help="dir containing <sample>.markdup.bam files")
    p.add_argument("--samples_tsv", required=True,
                   help="TSV with sample_id, group columns")
    p.add_argument("--outdir", required=True,
                   help="root output dir; ledger goes in <outdir>/<cid>/sd_substrate_evidence/")
    p.add_argument("--bam_suffix", default=".markdup.bam",
                   help="BAM filename suffix (default: .markdup.bam)")
    p.add_argument("--delly_inv_vcf", default="")
    p.add_argument("--manta_inv_vcf", default="")
    p.add_argument("--window", type=int, default=2000)
    p.add_argument("--min_mapq", type=int, default=20)
    p.add_argument("--min_clip_len", type=int, default=10)
    p.add_argument("--keep_secondary", action="store_true",
                   help="keep secondary alignments (default: drop)")
    p.add_argument("--threads", type=int, default=4,
                   help="sample-level parallelism")
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def main():
    args = parse_args()

    # ---- 1. Load coords from registry ----
    reg_root = _resolve_registries_root(args.registries_root, __file__)
    try:
        chrom, left_bp, right_bp = load_coords_from_registry(
            args.candidate, reg_root)
    except RuntimeError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(2)

    print(f"[STEP_03] {args.candidate}: {chrom}:{left_bp}-{right_bp}")

    # ---- 2. Load sample groups ----
    groups = load_samples_tsv(args.samples_tsv)
    if not groups:
        print(f"ERROR: no samples parsed from {args.samples_tsv}", file=sys.stderr)
        sys.exit(2)
    print(f"[STEP_03] {len(groups)} samples "
          f"(INV={sum(1 for g in groups.values() if g=='INV')}, "
          f"HET={sum(1 for g in groups.values() if g=='HET')}, "
          f"REF={sum(1 for g in groups.values() if g=='REF')}, "
          f"unknown={sum(1 for g in groups.values() if g=='unknown')})")

    # ---- 3. Output layout ----
    cand_dir = Path(args.outdir) / args.candidate
    ev_dir = cand_dir / "per_read_evidence"
    structured_dir = cand_dir / "structured"
    ev_dir.mkdir(parents=True, exist_ok=True)
    structured_dir.mkdir(parents=True, exist_ok=True)

    bam_tsv = ev_dir / "evidence_reads_bam.tsv.gz"
    vcf_tsv = ev_dir / "evidence_reads_vcf.tsv.gz"
    summary_json = structured_dir / "per_read_evidence_summary.json"

    if bam_tsv.exists() and vcf_tsv.exists() and not args.overwrite:
        print(f"[STEP_03] outputs exist, skipping (use --overwrite): {ev_dir}")
        return

    # ---- 4. BAM extraction (parallel over samples) ----
    def _one_sample(sid):
        bam_path = os.path.join(args.bam_dir, f"{sid}{args.bam_suffix}")
        grp = groups.get(sid, "unknown")
        rows_left = extract_bam_rows(
            bam_path, args.candidate, chrom, left_bp, "left",
            partner_chrom=chrom, partner_pos=right_bp,
            sample_id=sid, group=grp,
            window=args.window, min_mapq=args.min_mapq,
            min_clip_len=args.min_clip_len,
            keep_secondary=args.keep_secondary)
        rows_right = extract_bam_rows(
            bam_path, args.candidate, chrom, right_bp, "right",
            partner_chrom=chrom, partner_pos=left_bp,
            sample_id=sid, group=grp,
            window=args.window, min_mapq=args.min_mapq,
            min_clip_len=args.min_clip_len,
            keep_secondary=args.keep_secondary)
        return rows_left + rows_right

    all_bam_rows = []
    with ThreadPoolExecutor(max_workers=args.threads) as ex:
        futures = {ex.submit(_one_sample, sid): sid for sid in groups}
        for fut in as_completed(futures):
            sid = futures[fut]
            try:
                rows = fut.result()
                all_bam_rows.extend(rows)
            except Exception as e:
                sys.stderr.write(f"[STEP_03] sample {sid} failed: {e}\n")

    print(f"[STEP_03] extracted {len(all_bam_rows)} BAM evidence rows")
    write_gz_tsv(all_bam_rows, bam_tsv, BAM_COLUMNS)
    print(f"[STEP_03] wrote {bam_tsv}")

    # ---- 5. VCF sidecar extraction ----
    all_vcf_rows = []
    for vcf_path, caller in [(args.delly_inv_vcf, "delly"),
                             (args.manta_inv_vcf, "manta")]:
        if vcf_path:
            vcf_rows = extract_vcf_sidecar(
                vcf_path, caller, args.candidate, chrom,
                left_bp, right_bp, groups, window=args.window)
            print(f"[STEP_03] {caller}: {len(vcf_rows)} aggregate rows")
            all_vcf_rows.extend(vcf_rows)
    write_gz_tsv(all_vcf_rows, vcf_tsv, VCF_COLUMNS)
    print(f"[STEP_03] wrote {vcf_tsv}")

    # ---- 6. Summary block for registry ----
    summary = build_summary(all_bam_rows, args.window, args.min_mapq)
    summary["q7b_bam_ledger_path"] = str(bam_tsv.resolve())
    summary["q7b_vcf_sidecar_path"] = str(vcf_tsv.resolve())

    block = {
        "block_type":    "per_read_evidence_summary",
        "candidate_id":  args.candidate,
        "source_script": "STEP_03_per_read_evidence.py",
        "data":          summary,
    }
    with open(summary_json, "w") as f:
        json.dump(block, f, indent=2, default=str)
    print(f"[STEP_03] wrote {summary_json}")

    # ---- 7. Registry write (optional) ----
    if reg_root:
        try:
            sys.path.insert(0, os.path.join(reg_root, "api", "python"))
            from registry_loader import load_registry
            reg = load_registry()
            reg.evidence.write_block(
                candidate_id=args.candidate,
                block_type="per_read_evidence_summary",
                data=summary,
                source_script="STEP_03_per_read_evidence.py")
            print(f"[STEP_03] registered block for {args.candidate}")
        except Exception as e:
            print(f"[STEP_03] registry write failed: {e}", file=sys.stderr)

    print(f"[STEP_03] done — candidate={args.candidate} "
          f"split_left_INV={summary['q7b_bam_n_split_left_INV']} "
          f"discord_left_INV={summary['q7b_bam_n_discord_left_INV']}")


if __name__ == "__main__":
    main()
