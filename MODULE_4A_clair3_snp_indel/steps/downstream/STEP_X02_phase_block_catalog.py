#!/usr/bin/env python3
"""
02_phase_block_catalog.py — Collect phase blocks + prepare inversion marker data

Scans per-sample phase outputs and builds:
  - Phase block catalog across all samples
  - Variant → block membership map
  - Block-level haplotype signature matrix (for sample-distance MDS)
  - Inversion marker preparation table (bi-allelic het SNPs with block context)

This is the bridge between Clair3 phasing and the inversion/local-PCA pipeline.

Outputs:
  phase_prep/phase_block_catalog.tsv
  phase_prep/phase_block_variant_map.tsv
  phase_prep/block_signature_matrix.tsv        (samples × blocks, haplotype code)
  phase_prep/inversion_marker_candidates.tsv   (het SNPs with phase + position)
  phase_prep/phase_block_summary_per_sample.tsv
"""

import os, sys, argparse, time
from collections import defaultdict

def _now(): return time.time()
def _elapsed(t0):
    dt = time.time() - t0
    return f"{dt:.1f}s" if dt < 60 else f"{dt/60:.1f}min"

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--pp_results", required=True, help="postprocess_results/<CHROM>/ dir")
    p.add_argument("--chrom", required=True)
    p.add_argument("--outdir", required=True)
    return p.parse_args()

def main():
    args = parse_args()
    t0 = _now()
    outdir = os.path.join(args.outdir, "phase_prep")
    os.makedirs(outdir, exist_ok=True)

    chrom_dir = args.pp_results
    sample_dirs = []
    for d in sorted(os.listdir(chrom_dir)):
        full = os.path.join(chrom_dir, d)
        if os.path.isdir(full) and not d.startswith("_"):
            sample_dirs.append((d, full))

    print(f"[02] Found {len(sample_dirs)} sample directories")

    # ── Collect phase blocks and variant memberships ──
    all_blocks = []             # list of dicts
    all_var_map = []            # (sample, var_key, block_global_id, tier, phase_gt)
    per_sample_summary = []     # per-sample phase stats
    inversion_candidates = []   # het SNPs with phase info

    global_block_id = 0

    for si, (sample_id, sdir) in enumerate(sample_dirs):
        # Load phase blocks
        pb_path = os.path.join(sdir, "phase_blocks.tsv")
        phase_path = os.path.join(sdir, "all_variants_with_phase.tsv")

        n_t1 = n_t2 = n_unphased = 0
        n_blocks_t1 = n_blocks_t2 = 0
        sample_blocks = []

        if os.path.isfile(pb_path) and os.path.getsize(pb_path) > 0:
            with open(pb_path) as f:
                header = f.readline().strip().split('\t')
                col = {h: i for i, h in enumerate(header)}

                for line in f:
                    p = line.strip().split('\t')
                    if len(p) < len(header):
                        continue

                    tier = p[col["PHASE_TIER"]]
                    positions = p[col.get("POSITIONS", -1)] if "POSITIONS" in col else ""
                    n_var = int(p[col["N_VARIANTS"]])
                    span = int(p[col["BLOCK_SPAN"]])

                    block_entry = {
                        "GLOBAL_BLOCK_ID": global_block_id,
                        "SAMPLE": sample_id,
                        "CHROM": p[col["CHROM"]],
                        "BLOCK_START": int(p[col["BLOCK_START"]]),
                        "BLOCK_END": int(p[col["BLOCK_END"]]),
                        "BLOCK_SPAN": span,
                        "N_VARIANTS": n_var,
                        "PHASE_TIER": tier,
                        "LOCAL_BLOCK_ID": p[col["PHASE_BLOCK_ID"]],
                    }
                    all_blocks.append(block_entry)
                    sample_blocks.append((global_block_id, positions))

                    if tier == "TIER_1_WHATSHAP":
                        n_blocks_t1 += 1
                    elif tier == "TIER_2_READPAIR":
                        n_blocks_t2 += 1

                    global_block_id += 1

        # Load phased variants
        if os.path.isfile(phase_path) and os.path.getsize(phase_path) > 0:
            with open(phase_path) as f:
                header = f.readline().strip().split('\t')
                col = {h: i for i, h in enumerate(header)}

                for line in f:
                    p = line.strip().split('\t')
                    if len(p) < len(header):
                        continue

                    chrom = p[col["CHROM"]]
                    pos = p[col["POS"]]
                    ref = p[col.get("REF", -1)] if "REF" in col else "."
                    alt = p[col.get("ALT1", col.get("ALT", -1))] if ("ALT1" in col or "ALT" in col) else "."
                    gt = p[col.get("GT", col.get("PHASE_GT", -1))] if ("GT" in col or "PHASE_GT" in col) else "./."
                    phase_tier = p[col["PHASE_TIER"]] if "PHASE_TIER" in col else "UNPHASED"
                    phase_gt = p[col["PHASE_GT"]] if "PHASE_GT" in col else ""
                    block_id = p[col["PHASE_BLOCK_ID"]] if "PHASE_BLOCK_ID" in col else "-1"
                    var_type = p[col.get("VAR_TYPE", -1)] if "VAR_TYPE" in col else ""
                    gt_class = p[col.get("GT_CLASS", -1)] if "GT_CLASS" in col else ""

                    if phase_tier == "TIER_1_WHATSHAP":
                        n_t1 += 1
                    elif phase_tier == "TIER_2_READPAIR":
                        n_t2 += 1
                    else:
                        n_unphased += 1

                    var_key = f"{chrom}:{pos}:{ref}:{alt}"

                    # Map to global block ID
                    gbl_id = -1
                    if block_id != "-1":
                        local_bid = int(block_id)
                        # Find matching global block for this sample
                        for gbid, positions_str in sample_blocks:
                            if positions_str and pos in positions_str.split(","):
                                gbl_id = gbid
                                break

                    all_var_map.append({
                        "SAMPLE": sample_id,
                        "VAR_KEY": var_key,
                        "CHROM": chrom,
                        "POS": int(pos),
                        "REF": ref,
                        "ALT": alt,
                        "VAR_TYPE": var_type,
                        "PHASE_TIER": phase_tier,
                        "PHASE_GT": phase_gt,
                        "LOCAL_BLOCK_ID": block_id,
                        "GLOBAL_BLOCK_ID": gbl_id,
                    })

                    # Inversion marker candidate: phased het SNP
                    is_het = "|" in phase_gt and phase_gt in ("0|1", "1|0")
                    is_snp = len(ref) == 1 and len(alt) == 1
                    if is_het and is_snp and phase_tier == "TIER_1_WHATSHAP":
                        inversion_candidates.append({
                            "SAMPLE": sample_id,
                            "CHROM": chrom,
                            "POS": int(pos),
                            "REF": ref,
                            "ALT": alt,
                            "PHASE_GT": phase_gt,
                            "PHASE_TIER": phase_tier,
                            "LOCAL_BLOCK_ID": block_id,
                            "GLOBAL_BLOCK_ID": gbl_id,
                        })

        per_sample_summary.append({
            "SAMPLE": sample_id,
            "N_TIER1_VARIANTS": n_t1,
            "N_TIER2_VARIANTS": n_t2,
            "N_UNPHASED": n_unphased,
            "N_TIER1_BLOCKS": n_blocks_t1,
            "N_TIER2_BLOCKS": n_blocks_t2,
        })

        if (si + 1) % 50 == 0:
            print(f"[02]   Processed {si+1}/{len(sample_dirs)} samples …")

    # ── Write outputs ──

    # Phase block catalog
    out_blocks = os.path.join(outdir, "phase_block_catalog.tsv")
    with open(out_blocks, 'w') as o:
        cols = ["GLOBAL_BLOCK_ID", "SAMPLE", "CHROM", "BLOCK_START", "BLOCK_END",
                "BLOCK_SPAN", "N_VARIANTS", "PHASE_TIER", "LOCAL_BLOCK_ID"]
        o.write("\t".join(cols) + "\n")
        for b in all_blocks:
            o.write("\t".join(str(b[c]) for c in cols) + "\n")
    print(f"[02] Phase block catalog: {len(all_blocks)} blocks → {out_blocks}")

    # Variant → block map
    out_map = os.path.join(outdir, "phase_block_variant_map.tsv")
    with open(out_map, 'w') as o:
        cols = ["SAMPLE", "VAR_KEY", "CHROM", "POS", "REF", "ALT", "VAR_TYPE",
                "PHASE_TIER", "PHASE_GT", "LOCAL_BLOCK_ID", "GLOBAL_BLOCK_ID"]
        o.write("\t".join(cols) + "\n")
        for v in all_var_map:
            o.write("\t".join(str(v[c]) for c in cols) + "\n")
    print(f"[02] Variant-block map: {len(all_var_map)} entries → {out_map}")

    # Inversion marker candidates
    out_inv = os.path.join(outdir, "inversion_marker_candidates.tsv")
    with open(out_inv, 'w') as o:
        cols = ["SAMPLE", "CHROM", "POS", "REF", "ALT", "PHASE_GT",
                "PHASE_TIER", "LOCAL_BLOCK_ID", "GLOBAL_BLOCK_ID"]
        o.write("\t".join(cols) + "\n")
        for v in inversion_candidates:
            o.write("\t".join(str(v[c]) for c in cols) + "\n")
    print(f"[02] Inversion marker candidates: {len(inversion_candidates)} → {out_inv}")

    # Per-sample phase summary
    out_ps = os.path.join(outdir, "phase_block_summary_per_sample.tsv")
    with open(out_ps, 'w') as o:
        cols = ["SAMPLE", "N_TIER1_VARIANTS", "N_TIER2_VARIANTS", "N_UNPHASED",
                "N_TIER1_BLOCKS", "N_TIER2_BLOCKS"]
        o.write("\t".join(cols) + "\n")
        for s in per_sample_summary:
            o.write("\t".join(str(s[c]) for c in cols) + "\n")
    print(f"[02] Per-sample phase summary: {len(per_sample_summary)} → {out_ps}")

    print(f"\n[02] Total time: {_elapsed(t0)}")

if __name__ == "__main__":
    main()
