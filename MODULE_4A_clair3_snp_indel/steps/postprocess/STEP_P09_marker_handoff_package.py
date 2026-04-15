#!/usr/bin/env python3
"""
STEP09 – Generate marker handoff package ("Zoo folder").

Creates a clean, self-contained downstream marker package from the finished
Clair3 SNP/indel classification results.  Designed so another lab can open
one folder and understand what markers exist, what each is for, and which
are ready for assay development.

Inputs:
  --classified   final_variant_classification.tsv  (from STEP08)
  --phase        all_variants_with_phase.tsv       (from STEP02B, optional)
  --ref          reference.fa
  --outdir       01_Use_for_markers_multiplex/
  --sample_id    CGA009
  [--flank_bp    150]   flanking sequence for primer design
  [--gene_bed    gene_annotation.bed]  (optional, for gene overlap)

Outputs (in outdir):
  README_markers.txt
  marker_master_table.tsv
  marker_flanking_sequences.fa
  marker_coordinates.bed
  primer_design_notes.tsv
  marker_groups/
    high_priority_markers.tsv
    direct_snp_markers.tsv
    direct_indel_markers.tsv
    exploratory_markers.tsv

Usage:
  python STEP09_marker_handoff_package.py \
      --classified  final_variant_classification.tsv \
      --ref         reference.fa \
      --outdir      01_Use_for_markers_multiplex \
      --sample_id   CGA009 \
      [--phase      all_variants_with_phase.tsv] \
      [--flank_bp   150] \
      [--gene_bed   genes.bed]
"""

import os, sys, argparse
import numpy as np
import pandas as pd

try:
    import pysam
    HAS_PYSAM = True
except ImportError:
    HAS_PYSAM = False


# ── Marker classification logic ──────────────────────────────────────────────

def assign_marker_class(row):
    """Assign marker class A-E and assay logic.

    Class A – direct deployable: PASS or strong rescued, clean context
    Class B – deployable proxy: not directly assayable but proxy possible
    Class C – candidate: promising, needs validation
    Class D – exploratory only: kept for completeness
    Class E – not recommended: too weak/messy
    """
    fc = row.get("FINAL_CLASS", "")
    vt = row.get("VAR_TYPE", "")
    qual = row.get("QUAL", 0)
    dp = row.get("DP", 0)
    af1 = row.get("AF1", 0)
    ad_alt = row.get("AD_ALT", 0)
    n_local = row.get("N_LOCAL_VAR_10BP", 0)
    abs_ilen = row.get("ABS_INDEL_LEN", 0)

    # defaults
    marker_class = "D"
    marker_subclass = "exploratory_only"
    assay_logic = "exploratory_only"
    priority = "low"
    confidence = "low"
    deployment = "exploratory_only"
    direct_assay = "unknown"

    # ── Class A: direct deployable ──
    if fc == "STRICT_PASS" and vt == "SNP":
        if qual >= 30 and dp >= 5 and af1 >= 0.2 and n_local <= 2:
            marker_class = "A"
            marker_subclass = "direct_snp_high_conf"
            assay_logic = "direct_SNP"
            priority = "high"
            confidence = "high"
            deployment = "deployable"
            direct_assay = "yes"

    elif fc == "STRICT_PASS" and vt in ("INS", "DEL"):
        if qual >= 30 and dp >= 5 and af1 >= 0.2 and n_local <= 2:
            if pd.notna(abs_ilen) and abs_ilen <= 30:
                marker_class = "A"
                marker_subclass = "direct_indel_high_conf"
                assay_logic = "direct_indel"
                priority = "high"
                confidence = "high"
                deployment = "deployable"
                direct_assay = "yes"
            elif pd.notna(abs_ilen) and abs_ilen <= 100:
                marker_class = "B"
                marker_subclass = "direct_indel_medium"
                assay_logic = "direct_indel"
                priority = "medium"
                confidence = "medium"
                deployment = "candidate"
                direct_assay = "conditional"

    # ── Class A/B: rescued strong ──
    elif fc == "RESCUED_STRONG_SINGLE_SAMPLE":
        rescue_reason = row.get("RESCUE_REASON", "")

        if vt == "SNP" and "snp_q12" in rescue_reason:
            if qual >= 15 and dp >= 4 and n_local <= 3:
                marker_class = "B"
                marker_subclass = "rescued_snp_strong"
                assay_logic = "direct_SNP"
                priority = "medium"
                confidence = "medium"
                deployment = "candidate"
                direct_assay = "yes"
            else:
                marker_class = "C"
                marker_subclass = "rescued_snp_conditional"
                assay_logic = "direct_SNP"
                priority = "low"
                confidence = "low"
                deployment = "candidate"
                direct_assay = "conditional"

        elif vt in ("INS", "DEL"):
            if qual >= 10 and dp >= 4 and af1 >= 0.25 and n_local <= 3:
                if pd.notna(abs_ilen) and abs_ilen <= 30:
                    marker_class = "B"
                    marker_subclass = "rescued_indel_strong"
                    assay_logic = "direct_indel"
                    priority = "medium"
                    confidence = "medium"
                    deployment = "candidate"
                    direct_assay = "conditional"
                else:
                    marker_class = "C"
                    marker_subclass = "rescued_indel_large"
                    assay_logic = "direct_indel"
                    priority = "low"
                    confidence = "low"
                    deployment = "candidate"
                    direct_assay = "conditional"
            else:
                marker_class = "C"
                marker_subclass = "rescued_indel_weak"
                assay_logic = "direct_indel"
                priority = "low"
                confidence = "low"
                deployment = "candidate"
                direct_assay = "conditional"

    # ── Class C: population-rescued ──
    elif fc == "RESCUED_POPULATION_SIGNATURE":
        marker_class = "C"
        marker_subclass = "population_rescued"
        assay_logic = "direct_indel" if vt in ("INS", "DEL") else "direct_SNP"
        priority = "low"
        confidence = "low"
        deployment = "candidate"
        direct_assay = "conditional"

    # ── Class D/E: weak or discarded ──
    elif fc == "WEAK_SIGNATURE_CANDIDATE":
        marker_class = "D"
        marker_subclass = "weak_candidate"
        assay_logic = "exploratory_only"
        priority = "none"
        confidence = "very_low"
        deployment = "exploratory_only"
        direct_assay = "no"

    elif fc in ("FILTERED_DISCARD", "COMPLEX_LOCAL_BLOCK_ONLY"):
        marker_class = "E"
        marker_subclass = "not_recommended"
        assay_logic = "not_recommended"
        priority = "none"
        confidence = "none"
        deployment = "not_recommended"
        direct_assay = "no"

    # ── PASS variants that didn't hit Class A (e.g. messy context) ──
    elif fc == "STRICT_PASS":
        marker_class = "B"
        marker_subclass = "pass_messy_context"
        assay_logic = "direct_SNP" if vt == "SNP" else "direct_indel"
        priority = "medium"
        confidence = "medium"
        deployment = "candidate"
        direct_assay = "conditional"

    return pd.Series({
        "MARKER_CLASS": marker_class,
        "MARKER_SUBCLASS": marker_subclass,
        "ASSAY_LOGIC": assay_logic,
        "PRIORITY": priority,
        "CONFIDENCE": confidence,
        "DEPLOYMENT_STATUS": deployment,
        "DIRECT_ASSAY_POSSIBLE": direct_assay,
    })


def generate_marker_id(row, idx):
    """Generate human-readable marker ID."""
    vt = row.get("VAR_TYPE", "UNK")
    chrom = str(row.get("CHROM", "unk"))
    pos = int(row.get("POS", 0))
    # shorten chrom name
    chrom_short = chrom.replace("C_gar_LG", "LG")
    return f"MRK_{chrom_short}_{pos}_{vt}_{idx:06d}"


# ── Flanking sequence extraction ─────────────────────────────────────────────

def extract_flanking(ref_fa, chrom, pos, ref_allele, alt_allele, flank_bp=150):
    """Extract flanking sequence for primer design."""
    if ref_fa is None:
        return "", "", ""
    try:
        pos0 = pos - 1
        chrom_len = ref_fa.get_reference_length(chrom)
        left_start = max(0, pos0 - flank_bp)
        left = ref_fa.fetch(chrom, left_start, pos0).upper()
        right_start = pos0 + len(ref_allele)
        right_end = min(chrom_len, right_start + flank_bp)
        right = ref_fa.fetch(chrom, right_start, right_end).upper()
        return left, ref_allele, right
    except Exception:
        return "", "", ""


# ── Primer design notes ──────────────────────────────────────────────────────

def assess_primer_feasibility(row):
    """Quick assessment of PCR feasibility."""
    n_local_10 = row.get("N_LOCAL_VAR_10BP", 0)
    n_local_20 = row.get("N_LOCAL_VAR_20BP", 0)
    abs_ilen = row.get("ABS_INDEL_LEN", 0)
    vt = row.get("VAR_TYPE", "")

    if vt == "SNP":
        preferred_assay = "KASP_or_TaqMan"
    elif pd.notna(abs_ilen) and abs_ilen <= 5:
        preferred_assay = "KASP_or_fragment_analysis"
    elif pd.notna(abs_ilen) and abs_ilen <= 30:
        preferred_assay = "fragment_analysis"
    else:
        preferred_assay = "gel_or_capillary"

    easy = "yes" if n_local_10 <= 1 and n_local_20 <= 3 else "no"
    hard = "yes" if n_local_10 > 3 or n_local_20 > 6 else "no"

    caution = []
    if n_local_10 > 2:
        caution.append("high_local_variant_density")
    if pd.notna(abs_ilen) and abs_ilen > 50:
        caution.append("large_indel")
    if row.get("IS_OVERLAPPING", 0) == 1:
        caution.append("overlapping_events")

    return pd.Series({
        "PREFERRED_ASSAY": preferred_assay,
        "LIKELY_EASY_PCR": easy,
        "LIKELY_HARD_PCR": hard,
        "CAUTION_NOTE": "; ".join(caution) if caution else "none",
    })


# ── README ───────────────────────────────────────────────────────────────────

README_TEXT = """
================================================================================
 MARKER HANDOFF PACKAGE – SNPs and Indels
 Clair3 post-processing pipeline
================================================================================

PURPOSE:
  This folder contains a curated marker package from the Clair3 SNP/indel
  discovery and rescue pipeline.  It is designed for downstream use in
  multiplex marker development, PCR assay design, routine screening, or
  applied genotyping.

WHAT IS INCLUDED:
  - SNP markers (strict PASS + rescued QUAL >= 12)
  - Small indel markers (strict PASS + rescued)
  - Marker classification (Class A–E) with confidence and priority
  - Flanking sequences for primer design
  - Primer feasibility notes
  - Coordinate files for genome browsers

MARKER CLASSES:
  Class A = Direct deployable, high confidence
  Class B = Deployable (direct or proxy), medium confidence
  Class C = Candidate, needs validation
  Class D = Exploratory only
  Class E = Not recommended

IMPORTANT WARNINGS:
  1. Class D and E markers are NOT recommended for routine use.
  2. Low-coverage (~5x) data means some markers may have uncertain genotypes.
  3. Indels in repetitive regions may be unstable for PCR.
  4. Always validate candidate markers (Class B/C) before deployment.
  5. Phase information is LOCAL only (short blocks, not chromosome-scale).

FILE DESCRIPTIONS:
  marker_master_table.tsv       One row per marker with all metadata
  marker_flanking_sequences.fa  FASTA with flanking + alleles for primer design
  marker_coordinates.bed        BED file for IGV / genome browsers
  primer_design_notes.tsv       Feasibility notes for primer development
  marker_groups/                Subgroup tables by class and type

FOR WET-LAB USE:
  Start with marker_master_table.tsv, filter to Class A or B,
  then use marker_flanking_sequences.fa for primer design.

FOR BIOINFORMATICS:
  Use marker_coordinates.bed or marker_master_table.tsv for intersections.

TRACEABILITY:
  Each marker retains FINAL_CLASS and RESCUE_REASON from the discovery pipeline.
  marker_id format: MRK_<chrom>_<pos>_<type>_<index>

GENERATED BY: STEP09_marker_handoff_package.py
================================================================================
"""


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description="STEP09: Marker handoff package")
    ap.add_argument("--classified", required=True,
                    help="final_variant_classification.tsv")
    ap.add_argument("--phase", default=None,
                    help="all_variants_with_phase.tsv (optional)")
    ap.add_argument("--ref", required=True, help="Reference FASTA")
    ap.add_argument("--outdir", required=True,
                    help="Output directory (e.g. 01_Use_for_markers_multiplex)")
    ap.add_argument("--sample_id", default="SAMPLE")
    ap.add_argument("--flank_bp", type=int, default=150)
    ap.add_argument("--gene_bed", default=None, help="Gene annotation BED (optional)")
    # only include deployable markers (Class A-C) in main outputs
    ap.add_argument("--include_exploratory", action="store_true",
                    help="Also include Class D/E in outputs")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    groups_dir = os.path.join(args.outdir, "marker_groups")
    os.makedirs(groups_dir, exist_ok=True)

    # load data
    df = pd.read_csv(args.classified, sep="\t")
    print(f"[STEP09] Input variants: {len(df)}")

    # merge phase info if available
    if args.phase and os.path.isfile(args.phase):
        phase_df = pd.read_csv(args.phase, sep="\t")
        phase_cols = ["CHROM", "POS", "IS_PHASED", "PHASE_GT", "PS",
                      "PHASE_BLOCK_ID", "PHASE_TIER"]
        phase_cols = [c for c in phase_cols if c in phase_df.columns]
        if len(phase_cols) > 2:
            df = df.merge(
                phase_df[phase_cols],
                on=["CHROM", "POS"], how="left", suffixes=("", "_phase")
            )
            print(f"[STEP09] Phase info merged")

    # exclude HOM_REF and OTHER genotypes
    df = df[df["GT_CLASS"].isin(["HET", "HOM_VAR"])].copy()
    print(f"[STEP09] After genotype filter (HET/HOM_VAR): {len(df)}")

    # exclude Class E unless requested
    print("[STEP09] Assigning marker classes …")
    marker_info = df.apply(assign_marker_class, axis=1)
    df = pd.concat([df, marker_info], axis=1)

    if not args.include_exploratory:
        df = df[df["MARKER_CLASS"].isin(["A", "B", "C"])].copy()
        print(f"[STEP09] After class filter (A/B/C): {len(df)}")

    # assign marker IDs
    df = df.reset_index(drop=True)
    df["MARKER_ID"] = [generate_marker_id(row, i) for i, (_, row) in enumerate(df.iterrows())]

    # ── marker_master_table.tsv ──
    master_cols = [
        "MARKER_ID", "MARKER_CLASS", "MARKER_SUBCLASS", "ASSAY_LOGIC",
        "PRIORITY", "CONFIDENCE", "DEPLOYMENT_STATUS", "DIRECT_ASSAY_POSSIBLE",
        "FINAL_CLASS", "FINAL_LAYER", "RESCUE_REASON",
        "CHROM", "POS", "END", "REF", "ALT1", "VAR_TYPE",
        "INDEL_LEN", "ABS_INDEL_LEN",
        "QUAL", "GQ", "DP", "AD_REF", "AD_ALT", "AF1",
        "GT", "GT_CLASS",
        "N_LOCAL_VAR_10BP", "N_LOCAL_VAR_20BP", "IS_OVERLAPPING",
        "BLOCK_ID",
    ]
    # add phase cols if present
    for pc in ["IS_PHASED", "PHASE_TIER", "PHASE_BLOCK_ID"]:
        if pc in df.columns:
            master_cols.append(pc)

    available_cols = [c for c in master_cols if c in df.columns]
    master = df[available_cols].copy()

    out_master = os.path.join(args.outdir, "marker_master_table.tsv")
    master.to_csv(out_master, sep="\t", index=False)
    print(f"[STEP09] Master table → {out_master}  ({len(master)} markers)")

    # ── marker_flanking_sequences.fa ──
    ref_fa = None
    if HAS_PYSAM:
        try:
            ref_fa = pysam.FastaFile(args.ref)
        except Exception as e:
            print(f"[WARN] Cannot open ref: {e}", file=sys.stderr)

    out_fa = os.path.join(args.outdir, "marker_flanking_sequences.fa")
    with open(out_fa, "w") as fh:
        for _, row in df.iterrows():
            left, ref_seq, right = extract_flanking(
                ref_fa, row["CHROM"], int(row["POS"]),
                str(row["REF"]), str(row["ALT1"]), args.flank_bp
            )
            if not left and not right:
                continue
            header = (f">{row['MARKER_ID']} class={row['MARKER_CLASS']} "
                      f"type={row['VAR_TYPE']} {row['CHROM']}:{row['POS']} "
                      f"ref={row['REF']} alt={row['ALT1']}")
            # format: LEFT[REF/ALT]RIGHT
            seq = f"{left}[{row['REF']}/{row['ALT1']}]{right}"
            fh.write(f"{header}\n{seq}\n")

    if ref_fa:
        ref_fa.close()
    print(f"[STEP09] Flanking sequences → {out_fa}")

    # ── marker_coordinates.bed ──
    out_bed = os.path.join(args.outdir, "marker_coordinates.bed")
    with open(out_bed, "w") as fh:
        fh.write("#chrom\tstart\tend\tmarker_id\tclass\tstrand\n")
        for _, row in df.iterrows():
            start = int(row["POS"]) - 1
            end = int(row["END"])
            fh.write(f"{row['CHROM']}\t{start}\t{end}\t{row['MARKER_ID']}\t"
                      f"{row['MARKER_CLASS']}\t.\n")
    print(f"[STEP09] Coordinates → {out_bed}")

    # ── primer_design_notes.tsv ──
    primer_info = df.apply(assess_primer_feasibility, axis=1)
    primer_df = pd.concat([
        df[["MARKER_ID", "CHROM", "POS", "VAR_TYPE", "ABS_INDEL_LEN",
            "MARKER_CLASS", "N_LOCAL_VAR_10BP", "N_LOCAL_VAR_20BP",
            "IS_OVERLAPPING"]],
        primer_info
    ], axis=1)

    out_primer = os.path.join(args.outdir, "primer_design_notes.tsv")
    primer_df.to_csv(out_primer, sep="\t", index=False)
    print(f"[STEP09] Primer notes → {out_primer}")

    # ── marker_groups/ ──
    group_filters = {
        "high_priority_markers.tsv": df["PRIORITY"] == "high",
        "direct_snp_markers.tsv": (df["ASSAY_LOGIC"] == "direct_SNP") & (df["MARKER_CLASS"].isin(["A", "B"])),
        "direct_indel_markers.tsv": (df["ASSAY_LOGIC"] == "direct_indel") & (df["MARKER_CLASS"].isin(["A", "B"])),
        "exploratory_markers.tsv": df["MARKER_CLASS"].isin(["D"]),
    }

    for fname, mask in group_filters.items():
        sub = df[mask]
        if len(sub) > 0:
            sub[available_cols + ["MARKER_ID"]].to_csv(
                os.path.join(groups_dir, fname), sep="\t", index=False
            )
            print(f"[STEP09] Group: {fname} ({len(sub)} markers)")

    # ── README ──
    out_readme = os.path.join(args.outdir, "README_markers.txt")
    with open(out_readme, "w") as fh:
        fh.write(README_TEXT)
    print(f"[STEP09] README → {out_readme}")

    # ── Summary ──
    print("\n[STEP09] ═══════ MARKER PACKAGE SUMMARY ═══════")
    class_counts = df["MARKER_CLASS"].value_counts().sort_index()
    for cls, n in class_counts.items():
        print(f"  Class {cls}: {n:>8d}")
    print(f"  TOTAL:   {len(df):>8d}")

    type_counts = df.groupby(["MARKER_CLASS", "VAR_TYPE"]).size().reset_index(name="N")
    print("\n  By class × type:")
    for _, row in type_counts.iterrows():
        print(f"    {row['MARKER_CLASS']} / {row['VAR_TYPE']:5s}: {row['N']:>6d}")


if __name__ == "__main__":
    main()
