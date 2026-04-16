#!/usr/bin/env python3
"""
STEP02B – Phase-aware local block detection  (v3 – optimized)

Additive upgrade to the Clair3 post-processing workflow.
Does NOT replace STEP02 (local event blocks).  Adds a parallel phase layer.

Strategy:
  1. Extract phase information from Clair3 VCF:
     - Clair3 with --use_whatshap_for_final_output_phasing writes phased GTs
       (pipe-separated, e.g. 0|1) and PS (phase set) tags when WhatsHap succeeds.
     - PS = start position of the phase set (WhatsHap convention)
  2. Build local phase blocks from PS tags where available.
  3. For unphased heterozygous variants, attempt approximate local phase grouping
     using a REGION-CENTRIC read-pair linkage approach (requires pysam).
  4. Merge phase blocks with STEP02 event blocks for enriched annotation.

Phase confidence tiers:
  TIER_1_WHATSHAP  = WhatsHap PS tag present, phased GT (pipe-separated)
  TIER_2_READPAIR  = Approximate phase from read-pair linkage (low-cov)
  UNPHASED         = No phase information

v3 changes vs v2:
  - REMOVED TIER_3_PROXIMITY entirely (not real phasing).
  - Rewrote read-pair rescue: region-centric single-pass BAM scan per chunk
    instead of O(n^2) pairwise BAM fetches.  Uses union-find on shared reads.
  - Only considers HETEROZYGOUS UNPHASED variants for read-pair rescue.
  - Added comprehensive timing, progress, and debug logging.
  - Added checkpoint files per substep.

IMPORTANT CAVEAT for low-coverage Illumina:
  WhatsHap phase blocks in ~5x data will be SHORT (often 2-5 variants).
  This is expected and honest.  We do not overclaim long-range phase.
  The value is in local consistency: nearby variants on the same reads
  are very likely on the same haplotype.

IMPORTANT RULE:
  If the VCF contains NO true phased GTs at all (no 0|1 / 1|0 etc.),
  this script will NOT attempt read-pair fallback phasing.

Outputs:
  phase_blocks.tsv              – one row per phase block with summary
  all_variants_with_phase.tsv   – input table + phase columns added
  phase_block_event_overlap.tsv – overlap between phase blocks and event blocks
  phase_status_summary.tsv      – one-row summary of phase availability

Usage:
  python STEP02B_phase_aware_blocks.py \\
      --vcf         input.vcf.gz \\
      --annotated   all_variants_with_blocks.tsv \\
      --bam         sample.bam \\
      --outdir      output_dir \\
      --sample      SAMPLE_NAME \\
      [--read_phase_dist 500]  # max distance for read-pair phase inference
      [--chunk_size 50000]     # BAM scan chunk size
"""

import os, sys, argparse, time
from collections import defaultdict

import numpy as np
import pandas as pd

try:
    import pysam
    HAS_PYSAM = True
except ImportError:
    HAS_PYSAM = False


# ── Timing helper ────────────────────────────────────────────────────────────

def _now():
    return time.time()

def _elapsed(t0):
    dt = time.time() - t0
    if dt < 60:
        return f"{dt:.1f}s"
    return f"{dt/60:.1f}min"


def _checkpoint(outdir, name):
    """Write a small checkpoint marker so partial progress is visible."""
    path = os.path.join(outdir, f".checkpoint_{name}")
    with open(path, "w") as fh:
        fh.write(f"{name}\t{time.strftime('%Y-%m-%d %H:%M:%S')}\n")


# ── Extract phase info from VCF ──────────────────────────────────────────────

def extract_phase_from_vcf(vcf_path, sample_name):
    """Parse VCF for phased GTs and PS tags.  Returns DataFrame."""
    t0 = _now()
    print(f"[STEP02B]   Parsing VCF: {vcf_path}")

    vf = pysam.VariantFile(vcf_path)
    samples = list(vf.header.samples)
    if sample_name not in samples:
        if len(samples) == 1:
            print(f"[STEP02B]   WARN: sample '{sample_name}' not in VCF header; using '{samples[0]}'")
            sample_name = samples[0]
        else:
            print(f"[STEP02B]   WARN: sample '{sample_name}' not in VCF header; available: {samples}")

    records = []
    n_parsed = 0
    for rec in vf.fetch():
        n_parsed += 1
        s = rec.samples[sample_name]

        gt_tuple = s.get("GT", None)
        if gt_tuple is None:
            continue

        is_phased = s.phased if hasattr(s, 'phased') else False

        ps = None
        try:
            ps = s.get("PS", None)
        except Exception:
            pass

        gt_str = ("|" if is_phased else "/").join(
            str(x) if x is not None else "." for x in gt_tuple
        )

        # Determine zygosity
        alleles = [x for x in gt_tuple if x is not None]
        is_het = len(set(alleles)) > 1 if len(alleles) >= 2 else False

        records.append({
            "CHROM": rec.chrom,
            "POS": int(rec.pos),
            "END": int(rec.stop),
            "REF": rec.ref,
            "ALT1": rec.alts[0] if rec.alts else ".",
            "IS_PHASED": 1 if is_phased else 0,
            "PHASE_GT": gt_str,
            "PS": int(ps) if ps is not None else -1,
            "IS_HET": 1 if is_het else 0,
        })

        if n_parsed % 50000 == 0:
            print(f"[STEP02B]     … parsed {n_parsed} VCF records …")

    vf.close()
    df = pd.DataFrame(records)
    print(f"[STEP02B]   VCF parsing done: {len(df)} records in {_elapsed(t0)}")
    return df


# ── Build phase blocks from PS tags ──────────────────────────────────────────

def build_ps_phase_blocks(phase_df):
    """Group variants by (CHROM, PS) into WhatsHap phase blocks."""
    t0 = _now()
    blocks = []
    block_id = 0

    phased = phase_df[phase_df["PS"] >= 0].copy()
    if phased.empty:
        print(f"[STEP02B]   No PS-tagged variants found → 0 WhatsHap blocks")
        return pd.DataFrame(), block_id

    for (chrom, ps), grp in phased.groupby(["CHROM", "PS"]):
        g = grp.sort_values("POS")
        blocks.append({
            "PHASE_BLOCK_ID": block_id,
            "CHROM": chrom,
            "PS": int(ps),
            "BLOCK_START": int(g["POS"].min()),
            "BLOCK_END": int(g["END"].max()),
            "BLOCK_SPAN": int(g["END"].max()) - int(g["POS"].min()),
            "N_VARIANTS": len(g),
            "N_PHASED": int(g["IS_PHASED"].sum()),
            "PHASE_TIER": "TIER_1_WHATSHAP",
            "POSITIONS": ",".join(g["POS"].astype(str)),
        })
        block_id += 1

    bdf = pd.DataFrame(blocks)
    print(f"[STEP02B]   WhatsHap blocks built: {len(bdf)} in {_elapsed(t0)}")
    return bdf, block_id


# ── Union-Find data structure ────────────────────────────────────────────────

class UnionFind:
    """Efficient union-find with path compression and union-by-rank."""
    def __init__(self):
        self.parent = {}
        self.rank = {}

    def add(self, x):
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x] = 0

    def find(self, x):
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        # path compression
        while self.parent[x] != root:
            nxt = self.parent[x]
            self.parent[x] = root
            x = nxt
        return root

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1

    def groups(self):
        """Return dict: root → list of members (only groups of size ≥ 2)."""
        g = defaultdict(list)
        for x in self.parent:
            g[self.find(x)].append(x)
        return {k: sorted(v) for k, v in g.items() if len(v) >= 2}


# ── Region-centric read-pair phase rescue ────────────────────────────────────

def _region_centric_readpair(bam_path, positions_sorted, chrom, max_dist, chunk_size):
    """
    Efficient region-centric read-pair linkage.

    Instead of fetching reads per variant-pair (O(n^2) BAM seeks), we:
      1. Divide candidate positions into chunks based on gaps and max span.
      2. For each chunk, do ONE bam.fetch() over the chunk region.
      3. For each read, record which candidate positions it spans.
      4. If a read spans ≥2 candidates within max_dist, union them.

    Returns dict: root_pos → [member_positions]  (groups of size ≥ 2)
    """
    if not positions_sorted:
        return {}

    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        print(f"[STEP02B]   WARN: cannot open BAM: {e}")
        return {}

    uf = UnionFind()
    for p in positions_sorted:
        uf.add(p)

    n_reads_scanned = 0
    n_linking_reads = 0
    n_chunks = 0

    # Build chunks: contiguous groups where gap ≤ max_dist, split at chunk_size
    chunks = []
    chunk_start = 0
    for i in range(1, len(positions_sorted)):
        gap = positions_sorted[i] - positions_sorted[i - 1]
        chunk_len = positions_sorted[i] - positions_sorted[chunk_start]
        if gap > max_dist or chunk_len > chunk_size:
            chunks.append((chunk_start, i - 1))
            chunk_start = i
    chunks.append((chunk_start, len(positions_sorted) - 1))

    t0 = _now()
    for ci, (si, ei) in enumerate(chunks):
        region_start = positions_sorted[si]
        region_end = positions_sorted[ei]
        local_positions = positions_sorted[si : ei + 1]

        if len(local_positions) < 2:
            n_chunks += 1
            continue

        # Convert to sorted numpy array for fast lookup
        local_arr = np.array(local_positions, dtype=np.int64)

        # Fetch reads covering the region (with padding for read length)
        fetch_start = max(0, region_start - 1 - 300)  # 0-based, pad for reads starting before
        fetch_end = region_end + 300

        try:
            for read in bam.fetch(chrom, fetch_start, fetch_end):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if read.reference_start is None:
                    continue

                n_reads_scanned += 1
                rstart = read.reference_start      # 0-based inclusive
                rend = read.reference_end           # 0-based exclusive
                if rend is None:
                    continue

                # Which candidate positions does this read span?
                # VCF POS is 1-based; pysam coords are 0-based
                # Read covers 0-based [rstart, rend), so 1-based [rstart+1, rend]
                # A variant at POS p (1-based) is covered if rstart < p <= rend
                # equivalently: p-1 >= rstart and p-1 < rend  →  p > rstart and p <= rend
                lo = np.searchsorted(local_arr, rstart + 1, side='left')   # first p > rstart
                hi = np.searchsorted(local_arr, rend, side='right')        # first p > rend
                hit_positions = local_arr[lo:hi]

                if len(hit_positions) >= 2:
                    n_linking_reads += 1
                    # Union consecutive pairs within max_dist
                    anchor = int(hit_positions[0])
                    for k in range(1, len(hit_positions)):
                        hp = int(hit_positions[k])
                        if hp - anchor <= max_dist:
                            uf.union(anchor, hp)
                        anchor = hp

        except Exception as e:
            print(f"[STEP02B]   WARN: BAM fetch error chunk {ci}: {e}")
            continue

        n_chunks += 1
        if n_chunks % 500 == 0:
            print(f"[STEP02B]     … processed {n_chunks}/{len(chunks)} chunks, "
                  f"{n_reads_scanned} reads scanned, {n_linking_reads} linking reads … "
                  f"({_elapsed(t0)} elapsed)")

    bam.close()

    groups = uf.groups()

    print(f"[STEP02B]   Read-pair scan complete: "
          f"{n_chunks} chunks, {n_reads_scanned} reads scanned, "
          f"{n_linking_reads} linking reads, {len(groups)} groups built "
          f"({_elapsed(t0)})")

    return groups


def build_readpair_phase_blocks(bam_path, unphased_het_df, next_block_id,
                                max_dist=500, chunk_size=50000):
    """Build approximate phase blocks from read-pair linkage for unphased het variants.

    v3: Uses region-centric single-pass BAM scan instead of pairwise fetches.
    Only considers heterozygous unphased variants.
    """
    blocks = []
    block_id = next_block_id
    pos_to_block = {}
    t0 = _now()

    for chrom, grp in unphased_het_df.groupby("CHROM"):
        positions = sorted(grp["POS"].tolist())
        print(f"[STEP02B]   Read-pair rescue on {chrom}: {len(positions)} candidate het positions")

        if len(positions) < 2:
            print(f"[STEP02B]     Fewer than 2 candidates on {chrom}, skipping")
            continue

        if len(positions) > 100000:
            print(f"[STEP02B]   WARNING: {chrom} has {len(positions)} unphased het variants — "
                  f"this is very dense.  Scan may take a while.")

        groups = _region_centric_readpair(bam_path, positions, chrom,
                                           max_dist, chunk_size)

        for root_pos, members in groups.items():
            for p in members:
                pos_to_block[(chrom, p)] = block_id

            blocks.append({
                "PHASE_BLOCK_ID": block_id,
                "CHROM": chrom,
                "PS": -1,
                "BLOCK_START": min(members),
                "BLOCK_END": max(members),
                "BLOCK_SPAN": max(members) - min(members),
                "N_VARIANTS": len(members),
                "N_PHASED": 0,
                "PHASE_TIER": "TIER_2_READPAIR",
                "POSITIONS": ",".join(str(p) for p in members),
            })
            block_id += 1

    print(f"[STEP02B]   Read-pair block construction: {len(blocks)} blocks in {_elapsed(t0)}")
    return pd.DataFrame(blocks), pos_to_block, block_id


# ── Merge with event blocks ──────────────────────────────────────────────────

def compute_phase_event_overlap(phase_blocks, event_blocks):
    """Find overlaps between phase blocks and event blocks.
    Uses sorted intervals with binary search instead of O(n*m) nested loop.
    """
    t0 = _now()
    if phase_blocks.empty or event_blocks.empty:
        return pd.DataFrame()

    overlaps = []

    # Group by chrom for efficiency
    pb_by_chrom = {c: g.sort_values("BLOCK_START") for c, g in phase_blocks.groupby("CHROM")}
    eb_by_chrom = {c: g.sort_values("BLOCK_START") for c, g in event_blocks.groupby("CHROM")}

    for chrom in pb_by_chrom:
        if chrom not in eb_by_chrom:
            continue
        pbs = pb_by_chrom[chrom]
        ebs = eb_by_chrom[chrom]

        eb_starts = ebs["BLOCK_START"].values
        eb_ends = ebs["BLOCK_END"].values
        eb_ids = ebs["BLOCK_ID"].values
        eb_nrecs = ebs["N_RECORDS"].values

        for _, pb in pbs.iterrows():
            ps, pe = pb["BLOCK_START"], pb["BLOCK_END"]
            # Find event blocks that could overlap: eb_end >= ps AND eb_start <= pe
            lo = np.searchsorted(eb_ends, ps, side='left')
            hi = np.searchsorted(eb_starts, pe, side='right')
            for idx in range(lo, min(hi, len(eb_starts))):
                if eb_starts[idx] <= pe and eb_ends[idx] >= ps:
                    overlap_start = max(ps, eb_starts[idx])
                    overlap_end = min(pe, eb_ends[idx])
                    overlaps.append({
                        "PHASE_BLOCK_ID": pb["PHASE_BLOCK_ID"],
                        "EVENT_BLOCK_ID": int(eb_ids[idx]),
                        "CHROM": chrom,
                        "OVERLAP_START": int(overlap_start),
                        "OVERLAP_END": int(overlap_end),
                        "OVERLAP_SPAN": int(overlap_end - overlap_start),
                        "PHASE_TIER": pb["PHASE_TIER"],
                        "PHASE_N_VARIANTS": pb["N_VARIANTS"],
                        "EVENT_N_RECORDS": int(eb_nrecs[idx]),
                    })

    result = pd.DataFrame(overlaps)
    print(f"[STEP02B]   Overlap computation: {len(result)} overlaps in {_elapsed(t0)}")
    return result


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    t_global = _now()

    ap = argparse.ArgumentParser(description="STEP02B: Phase-aware local blocks (v3)")
    ap.add_argument("--vcf", required=True, help="Clair3 VCF (phased output)")
    ap.add_argument("--annotated", required=True,
                    help="all_variants_with_blocks.tsv from STEP02")
    ap.add_argument("--bam", default=None, help="BAM for read-pair phase inference")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--sample", default="SAMPLE")
    ap.add_argument("--read_phase_dist", type=int, default=500,
                    help="Max distance for read-pair phase inference [500]")
    ap.add_argument("--chunk_size", type=int, default=50000,
                    help="Max genomic span per BAM scan chunk [50000]")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # ── 1. Extract phase from VCF ────────────────────────────────────────────
    print("[STEP02B] Extracting phase info from VCF …")
    t1 = _now()
    phase_df = extract_phase_from_vcf(args.vcf, args.sample)
    n_total = len(phase_df)
    n_phased = int(phase_df["IS_PHASED"].sum()) if n_total > 0 else 0
    n_with_ps = int((phase_df["PS"] >= 0).sum()) if n_total > 0 else 0
    n_het = int(phase_df["IS_HET"].sum()) if n_total > 0 else 0
    n_hom = n_total - n_het
    print(f"[STEP02B] VCF records: {n_total}, phased GTs: {n_phased}, "
          f"with PS: {n_with_ps}")
    print(f"[STEP02B]   heterozygous: {n_het}, homozygous: {n_hom}")
    print(f"[STEP02B]   VCF extraction time: {_elapsed(t1)}")
    _checkpoint(args.outdir, "01_vcf_parsed")

    # ── 2. Build WhatsHap phase blocks ───────────────────────────────────────
    print("[STEP02B] Building WhatsHap phase blocks …")
    t2 = _now()
    ws_blocks, next_id = build_ps_phase_blocks(phase_df)
    print(f"[STEP02B] WhatsHap phase blocks: {len(ws_blocks)}")
    print(f"[STEP02B]   WhatsHap block time: {_elapsed(t2)}")
    _checkpoint(args.outdir, "02_whatshap_blocks")

    # ── 3. Read-pair approximate phase for unphased HET variants ─────────────
    unphased = phase_df[phase_df["PS"] < 0].copy()
    n_unphased = len(unphased)
    unphased_het = unphased[unphased["IS_HET"] == 1].copy()
    n_unphased_het = len(unphased_het)
    n_unphased_hom = n_unphased - n_unphased_het
    n_already_phased = n_total - n_unphased

    rp_blocks = pd.DataFrame()
    rp_pos_map = {}

    print(f"[STEP02B] Unphased variants: {n_unphased}")
    print(f"[STEP02B]   already phased (WhatsHap):       {n_already_phased}")
    print(f"[STEP02B]   unphased heterozygous:            {n_unphased_het}  <- candidates for read-pair rescue")
    print(f"[STEP02B]   unphased homozygous/missing:      {n_unphased_hom}  <- skipped (not informative)")

    if n_phased == 0:
        print("[STEP02B] No phased GTs found in VCF; skipping read-pair phase inference.")
        print("[STEP02B]   (This means WhatsHap did not phase any variants for this sample.)")
    elif not args.bam:
        print("[STEP02B] No BAM provided, skipping read-pair phase inference")
    elif not HAS_PYSAM:
        print("[STEP02B] pysam not available, skipping read-pair phase inference")
    elif n_unphased_het == 0:
        print("[STEP02B] No unphased het variants to rescue.")
    else:
        print(f"[STEP02B] Read-pair phase inference for {n_unphased_het} unphased het variants …")
        t3 = _now()
        rp_blocks, rp_pos_map, next_id = build_readpair_phase_blocks(
            args.bam, unphased_het, next_id,
            max_dist=args.read_phase_dist,
            chunk_size=args.chunk_size,
        )
        print(f"[STEP02B] Read-pair phase blocks: {len(rp_blocks)}")
        print(f"[STEP02B]   Read-pair total time: {_elapsed(t3)}")

    _checkpoint(args.outdir, "03_readpair_done")

    # ── 4. Combine all phase blocks ──────────────────────────────────────────
    all_phase_blocks = pd.concat([ws_blocks, rp_blocks], ignore_index=True)

    out_blocks = os.path.join(args.outdir, "phase_blocks.tsv")
    all_phase_blocks.to_csv(out_blocks, sep="\t", index=False)
    print(f"[STEP02B] Phase blocks -> {out_blocks}")

    # ── 5. Annotate variants with phase info ─────────────────────────────────
    print("[STEP02B] Annotating variants with phase columns …")
    t5 = _now()
    df = pd.read_csv(args.annotated, sep="\t")

    # Build fast lookup: (CHROM, POS) → phase info
    phase_lookup = {}
    for _, row in phase_df.iterrows():
        key = (row["CHROM"], int(row["POS"]))
        phase_lookup[key] = {
            "IS_PHASED": row["IS_PHASED"],
            "PHASE_GT": row["PHASE_GT"],
            "PS": row["PS"],
            "IS_HET": row["IS_HET"],
        }

    # Build variant → phase block and tier maps
    var_to_phase_block = {}
    var_to_phase_tier = {}

    if not ws_blocks.empty:
        for _, blk in ws_blocks.iterrows():
            positions = [int(p) for p in str(blk["POSITIONS"]).split(",") if p]
            for p in positions:
                var_to_phase_block[(blk["CHROM"], p)] = blk["PHASE_BLOCK_ID"]
                var_to_phase_tier[(blk["CHROM"], p)] = "TIER_1_WHATSHAP"

    for (chrom, pos), bid in rp_pos_map.items():
        if (chrom, pos) not in var_to_phase_block:
            var_to_phase_block[(chrom, pos)] = bid
            var_to_phase_tier[(chrom, pos)] = "TIER_2_READPAIR"

    # Annotate each variant row
    is_phased_col = []
    phase_gt_col = []
    ps_col = []
    phase_block_id_col = []
    phase_tier_col = []

    for _, row in df.iterrows():
        key = (row["CHROM"], int(row["POS"]))
        pl = phase_lookup.get(key, {})
        is_phased_col.append(pl.get("IS_PHASED", 0))
        phase_gt_col.append(pl.get("PHASE_GT", ""))
        ps_col.append(pl.get("PS", -1))
        phase_block_id_col.append(var_to_phase_block.get(key, -1))

        if key in var_to_phase_tier:
            phase_tier_col.append(var_to_phase_tier[key])
        elif pl.get("IS_PHASED", 0) == 1:
            phase_tier_col.append("TIER_1_WHATSHAP")
        else:
            phase_tier_col.append("UNPHASED")

    df["IS_PHASED"] = is_phased_col
    df["PHASE_GT"] = phase_gt_col
    df["PS"] = ps_col
    df["PHASE_BLOCK_ID"] = phase_block_id_col
    df["PHASE_TIER"] = phase_tier_col

    out_annot = os.path.join(args.outdir, "all_variants_with_phase.tsv")
    df.to_csv(out_annot, sep="\t", index=False)
    print(f"[STEP02B] Annotated variants -> {out_annot}")
    print(f"[STEP02B]   Annotation time: {_elapsed(t5)}")
    _checkpoint(args.outdir, "04_annotation_done")

    # ── 6. Phase–event block overlap ─────────────────────────────────────────
    t6 = _now()
    event_blocks_path = os.path.join(args.outdir, "local_overlap_blocks.tsv")
    if not os.path.isfile(event_blocks_path):
        parent = os.path.dirname(args.annotated)
        event_blocks_path = os.path.join(parent, "local_overlap_blocks.tsv")

    if os.path.isfile(event_blocks_path) and not all_phase_blocks.empty:
        event_blocks = pd.read_csv(event_blocks_path, sep="\t")
        overlap = compute_phase_event_overlap(all_phase_blocks, event_blocks)
        out_overlap = os.path.join(args.outdir, "phase_block_event_overlap.tsv")
        overlap.to_csv(out_overlap, sep="\t", index=False)
        print(f"[STEP02B] Phase-event overlap -> {out_overlap}  ({len(overlap)} overlaps)")
    else:
        print("[STEP02B] Skipping overlap (no event blocks file or no phase blocks)")

    print(f"[STEP02B]   Overlap time: {_elapsed(t6)}")
    _checkpoint(args.outdir, "05_overlap_done")

    # ── 7. Phase status summary ──────────────────────────────────────────────
    if n_phased > 0:
        phase_mode = "WHATSHAP_PRESENT"
    elif len(rp_blocks) > 0:
        phase_mode = "READPAIR_ONLY"
    else:
        phase_mode = "NO_PHASE_IN_VCF"

    status_df = pd.DataFrame([{
        "SAMPLE": args.sample,
        "VCF": args.vcf,
        "N_VCF_RECORDS": n_total,
        "N_PHASED_GT": n_phased,
        "N_WITH_PS": n_with_ps,
        "N_HET": n_het,
        "N_UNPHASED_HET": n_unphased_het,
        "N_WHATSHAP_BLOCKS": len(ws_blocks),
        "N_READPAIR_BLOCKS": len(rp_blocks),
        "PHASE_MODE": phase_mode,
    }])

    out_status = os.path.join(args.outdir, "phase_status_summary.tsv")
    status_df.to_csv(out_status, sep="\t", index=False)
    print(f"[STEP02B] Phase status summary -> {out_status}")

    # ── Summary ──────────────────────────────────────────────────────────────
    print()
    print("[STEP02B] ======= PHASE SUMMARY =======")
    for tier in ["TIER_1_WHATSHAP", "TIER_2_READPAIR", "UNPHASED"]:
        n = sum(1 for t in phase_tier_col if t == tier)
        pct = 100 * n / len(df) if len(df) > 0 else 0
        print(f"  {tier:25s}  {n:>8d}  ({pct:.1f}%)")
    print(f"  {'TOTAL':25s}  {len(df):>8d}")

    if not all_phase_blocks.empty:
        for tier_name in ["TIER_1_WHATSHAP", "TIER_2_READPAIR"]:
            sub = all_phase_blocks[all_phase_blocks["PHASE_TIER"] == tier_name]
            if len(sub) > 0:
                print(f"\n  {tier_name} blocks: {len(sub)}")
                print(f"    median variants/block: {sub['N_VARIANTS'].median():.0f}")
                print(f"    median span:           {sub['BLOCK_SPAN'].median():.0f} bp")
                print(f"    max span:              {sub['BLOCK_SPAN'].max():.0f} bp")
                # Block size distribution
                size_counts = sub["N_VARIANTS"].value_counts().sort_index()
                dist_str = ", ".join(f"{k}var:{v}" for k, v in size_counts.head(8).items())
                print(f"    size distribution:     {dist_str}")

    print()
    print("[STEP02B] IMPORTANT CAVEAT:")
    print("  Phase blocks from ~5x Illumina are SHORT (typically 2-5 variants).")
    print("  This is expected.  Do not overclaim long-range haplotype structure.")
    print("  Value is in LOCAL consistency: nearby variants on same reads = same haplotype.")

    print()
    print(f"[STEP02B] Total STEP02B time: {_elapsed(t_global)}")
    _checkpoint(args.outdir, "06_step02b_complete")


if __name__ == "__main__":
    main()
