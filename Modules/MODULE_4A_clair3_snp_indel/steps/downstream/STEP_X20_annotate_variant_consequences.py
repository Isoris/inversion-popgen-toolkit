#!/usr/bin/env python3
"""
20_annotate_variant_consequences.py — Variant consequence annotation layer

Takes the population-level variant catalog and annotates each variant with:
  - Gene overlap context (exon/CDS/intron/splice/intergenic)
  - bcftools csq consequence string (if available)
  - Splice proximity flag
  - Variant class assignment (SNP / small_indel / splice / DEL / DUP / INV / INS / BND)
  - Disruption severity tier (very_high / high / moderate / low / minimal)
  - Zygosity per sample (from GT matrix)
  - ROH overlap flag per sample (if ROH BED available)
  - Cohort recurrence (carrier count)

This script is the INPUT FEEDER for 21_score_breeding_concern.py.
It does NOT assign BC scores — only annotation layers.

Inputs:
  --catalog          {TYPE}_catalog.tsv from 01
  --gt_matrix        GT_matrix.{TYPE}.tsv from 01
  --functional_bed   BED with gene/exon/CDS features (from DELLY annotation layer or GFF-derived)
  --csq_tsv          Optional: pre-computed bcftools csq output
  --roh_dir          Optional: per-sample ROH BED files (from ngsF-HMM or bcftools roh)
  --chrom            Chromosome
  --outdir           Output directory

Outputs:
  vef/annotated_variants.{TYPE}.tsv   — one row per variant with all annotation columns
  vef/per_sample_genotypes.{TYPE}.tsv — one row per (variant × sample) with zygosity + ROH
"""

import os, sys, argparse, time, glob
from collections import defaultdict

def _now(): return time.time()
def _elapsed(t0):
    dt = time.time() - t0
    return f"{dt:.1f}s" if dt < 60 else f"{dt/60:.1f}min"

# ── Variant class assignment ─────────────────────────────────────────────────

def assign_variant_class(var_type, ref, alt, indel_len, csq_consequence=""):
    """Assign variant class from the class table (Section 1 of framework)."""
    vt = str(var_type).upper()
    csq = str(csq_consequence).lower()

    # Structural variants first
    if vt in ("DEL", "DUP", "INV", "INS", "BND", "TRA"):
        return vt

    # Splice events (detected from csq)
    if any(s in csq for s in ["splice_donor", "splice_acceptor", "splice_region"]):
        return "splice"

    # Small indel
    if vt in ("INDEL", "MNP") or (len(ref) != len(alt) and len(ref) <= 50 and len(alt) <= 50):
        return "small_indel"

    # SNP
    if len(ref) == 1 and len(alt) == 1:
        return "SNP"

    # Fallback
    try:
        if abs(int(indel_len)) > 0:
            return "small_indel"
    except (ValueError, TypeError):
        pass

    return "SNP" if vt == "SNP" else "unknown"


# ── Disruption severity ──────────────────────────────────────────────────────

# CSQ consequence → severity mapping
# Based on Ensembl VEP severity + framework Section 2A
CSQ_SEVERITY = {
    # Very high
    "stop_gained": "very_high",
    "frameshift_variant": "very_high",
    "start_lost": "very_high",
    "splice_donor_variant": "very_high",
    "splice_acceptor_variant": "very_high",
    "transcript_ablation": "very_high",
    "feature_truncation": "very_high",

    # High
    "stop_lost": "high",
    "initiator_codon_variant": "high",
    "inframe_insertion": "high",
    "inframe_deletion": "high",
    "disruptive_inframe_insertion": "high",
    "disruptive_inframe_deletion": "high",
    "protein_altering_variant": "high",

    # Moderate
    "missense_variant": "moderate",
    "splice_region_variant": "moderate",
    "incomplete_terminal_codon_variant": "moderate",
    "coding_sequence_variant": "moderate",

    # Low
    "synonymous_variant": "low",
    "5_prime_utr_variant": "low",
    "3_prime_utr_variant": "low",
    "non_coding_transcript_exon_variant": "low",
    "intron_variant": "low",
    "nmd_transcript_variant": "low",

    # Minimal
    "upstream_gene_variant": "minimal",
    "downstream_gene_variant": "minimal",
    "intergenic_variant": "minimal",
    "non_coding_transcript_variant": "minimal",
}

def assign_disruption_severity_small(csq_consequence, feature_context):
    """Assign disruption severity for small variants (Section 2A)."""
    csq = str(csq_consequence).lower()

    # Check each csq term (can be multiple, e.g. "missense_variant&splice_region_variant")
    terms = csq.replace("&", ",").replace("|", ",").split(",")
    worst = "minimal"
    severity_rank = {"very_high": 5, "high": 4, "moderate": 3, "low": 2, "minimal": 1}

    for term in terms:
        term = term.strip()
        sev = CSQ_SEVERITY.get(term, None)
        if sev and severity_rank.get(sev, 0) > severity_rank.get(worst, 0):
            worst = sev

    # If no csq available, use feature context
    if worst == "minimal" and csq in ("", ".", "none", "nan"):
        ctx = str(feature_context).lower()
        if "cds" in ctx:
            worst = "moderate"
        elif "exon" in ctx:
            worst = "moderate"
        elif "splice" in ctx:
            worst = "moderate"
        elif "intron" in ctx:
            worst = "low"

    return worst

def assign_disruption_severity_sv(feature_context, sv_type, sv_len=0):
    """Assign disruption severity for structural variants (Section 2B)."""
    ctx = str(feature_context).lower()
    svt = str(sv_type).upper()

    if "cds" in ctx or "exon" in ctx:
        if svt == "DEL":
            return "very_high"
        elif svt in ("DUP", "INV"):
            return "high"
        elif svt == "INS":
            return "high"
        else:
            return "high"
    elif "splice" in ctx:
        return "high"
    elif "intron" in ctx or "gene" in ctx:
        return "moderate"
    elif "promoter" in ctx:
        return "moderate"
    else:
        return "low"


# ── Splice proximity ────────────────────────────────────────────────────────

def check_splice_proximity(csq_consequence):
    """Return splice flag from csq consequence string."""
    csq = str(csq_consequence).lower()
    if "splice_donor" in csq or "splice_acceptor" in csq:
        return "essential_splice"
    if "splice_region" in csq:
        return "splice_region"
    return "none"


# ── Evidence codes ───────────────────────────────────────────────────────────

def assign_evidence_codes(variant_class, disruption_severity, csq_consequence,
                          splice_flag, feature_context):
    """Assign initial evidence codes (Section 3) based on annotation alone.
    Zygosity/ROH/ESM codes are added in the scoring step."""
    codes = []
    csq = str(csq_consequence).lower()

    # Strong disruption evidence
    if "frameshift" in csq:
        codes.append("SD1")
    if "stop_gained" in csq:
        codes.append("SD2")
    if splice_flag == "essential_splice":
        codes.append("SD3")
    if variant_class in ("DEL", "DUP", "INV", "INS") and feature_context in ("CDS_overlap", "exon_overlap"):
        codes.append("SD4")

    # Moderate disruption evidence
    if "missense" in csq:
        codes.append("MD1")
    if "inframe" in csq:
        codes.append("MD2")
    if splice_flag == "splice_region":
        codes.append("MD3")
    if variant_class in ("DEL", "DUP", "INV", "INS") and "intron" in str(feature_context).lower():
        codes.append("MD4")

    # Supportive
    # SU1 (ESM), SU2 (rare), SU6 (pathway), SU7 (ROH) added later

    return codes


# ── Main ─────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--catalog", required=True)
    p.add_argument("--gt_matrix", required=True)
    p.add_argument("--chrom", required=True)
    p.add_argument("--outdir", required=True)
    # Optional annotation inputs
    p.add_argument("--csq_tsv", default="", help="Pre-computed csq consequences TSV")
    p.add_argument("--functional_class_tsv", default="",
                   help="Functional class from DELLY annotation (del_key, func_class, gene_names)")
    p.add_argument("--roh_dir", default="", help="Directory with per-sample ROH BEDs")
    p.add_argument("--esm_tsv", default="", help="ESM missense scores TSV (var_key, esm_score)")
    return p.parse_args()

def main():
    args = parse_args()
    t0 = _now()
    outdir = os.path.join(args.outdir, "vef")
    os.makedirs(outdir, exist_ok=True)

    # ── Load catalog ──
    print("[VEF-20] Loading variant catalog …")
    variants = []
    with open(args.catalog) as f:
        header = f.readline().strip().split('\t')
        col = {h: i for i, h in enumerate(header)}
        for line in f:
            p = line.strip().split('\t')
            variants.append(dict(zip(header, p)))
    print(f"[VEF-20]   {len(variants)} variants loaded")

    # ── Load optional CSQ ──
    csq_lookup = {}
    if args.csq_tsv and os.path.isfile(args.csq_tsv):
        print(f"[VEF-20] Loading CSQ from {args.csq_tsv} …")
        with open(args.csq_tsv) as f:
            ch = f.readline().strip().split('\t')
            cc = {h: i for i, h in enumerate(ch)}
            for line in f:
                p = line.strip().split('\t')
                # Expect columns: CHROM, POS, REF, ALT, CSQ_CONSEQUENCE, GENE, FEATURE
                key = f"{p[cc.get('CHROM',0)]}:{p[cc.get('POS',1)]}:{p[cc.get('REF',2)]}:{p[cc.get('ALT',3)]}"
                csq_lookup[key] = {
                    "csq": p[cc.get("CSQ_CONSEQUENCE", cc.get("consequence", 4))],
                    "gene": p[cc.get("GENE", cc.get("gene", 5))] if len(p) > 5 else ".",
                    "feature": p[cc.get("FEATURE", cc.get("feature_context", 6))] if len(p) > 6 else ".",
                }
        print(f"[VEF-20]   CSQ loaded for {len(csq_lookup)} variants")

    # ── Load optional functional class ──
    func_lookup = {}
    if args.functional_class_tsv and os.path.isfile(args.functional_class_tsv):
        print(f"[VEF-20] Loading functional class …")
        with open(args.functional_class_tsv) as f:
            fh = f.readline().strip().split('\t')
            fc = {h: i for i, h in enumerate(fh)}
            for line in f:
                p = line.strip().split('\t')
                # Try to match by del_id or by chrom:pos
                for key_col in ["del_key", "id", "del_id"]:
                    if key_col in fc:
                        func_lookup[p[fc[key_col]]] = {
                            "func_class": p[fc.get("class", fc.get("func_class", -1))],
                            "gene_names": p[fc.get("gene_names", -1)] if "gene_names" in fc else ".",
                        }

    # ── Load optional ESM scores ──
    esm_lookup = {}
    if args.esm_tsv and os.path.isfile(args.esm_tsv):
        print(f"[VEF-20] Loading ESM scores …")
        with open(args.esm_tsv) as f:
            eh = f.readline().strip().split('\t')
            ec = {h: i for i, h in enumerate(eh)}
            for line in f:
                p = line.strip().split('\t')
                vk = p[ec.get("VAR_KEY", ec.get("var_key", 0))]
                score = p[ec.get("ESM_SCORE", ec.get("esm_score", 1))]
                esm_lookup[vk] = float(score) if score not in (".", "NA", "") else None
        print(f"[VEF-20]   ESM scores for {len(esm_lookup)} variants")

    # ── Load ROH regions ──
    roh_regions = {}  # sample → list of (start, end)
    if args.roh_dir and os.path.isdir(args.roh_dir):
        print(f"[VEF-20] Loading ROH BEDs …")
        for bed_file in glob.glob(os.path.join(args.roh_dir, "*.roh.bed")) + \
                        glob.glob(os.path.join(args.roh_dir, "*", "*.roh.bed")):
            sample_id = os.path.basename(bed_file).replace(".roh.bed", "")
            regions = []
            with open(bed_file) as f:
                for line in f:
                    p = line.strip().split('\t')
                    if p[0] == args.chrom:
                        regions.append((int(p[1]), int(p[2])))
            if regions:
                roh_regions[sample_id] = sorted(regions)
        print(f"[VEF-20]   ROH data for {len(roh_regions)} samples")

    # ── Load GT matrix for zygosity ──
    print("[VEF-20] Loading GT matrix …")
    gt_samples = []
    gt_data = {}  # var_key → {sample: gt_string}
    with open(args.gt_matrix) as f:
        gh = f.readline().strip().split('\t')
        gt_samples = gh[5:]
        for line in f:
            p = line.strip().split('\t')
            vk = p[4]
            gts = p[5:]
            gt_data[vk] = dict(zip(gt_samples, gts))
    print(f"[VEF-20]   GT data for {len(gt_data)} variants × {len(gt_samples)} samples")

    # ── Annotate each variant ──
    print("[VEF-20] Annotating variants …")

    out_var_path = os.path.join(outdir, f"annotated_variants.tsv")
    out_gt_path = os.path.join(outdir, f"per_sample_genotypes.tsv")

    var_cols = [
        "VAR_KEY", "CHROM", "POS", "REF", "ALT",
        "VARIANT_CLASS", "VAR_TYPE_ORIG",
        "GENE", "FEATURE_CONTEXT", "CSQ_CONSEQUENCE", "SPLICE_FLAG",
        "ESM_SCORE",
        "DISRUPTION_SEVERITY",
        "EVIDENCE_CODES",
        "N_CARRIERS_ALL", "CARRIER_FREQ_ALL", "FREQ_CLASS",
        "INDEL_LEN", "ABS_INDEL_LEN",
    ]

    gt_cols = [
        "VAR_KEY", "SAMPLE", "GT", "ZYGOSITY", "IN_ROH",
    ]

    n_annotated = 0
    with open(out_var_path, 'w') as fvar, open(out_gt_path, 'w') as fgt:
        fvar.write("\t".join(var_cols) + "\n")
        fgt.write("\t".join(gt_cols) + "\n")

        for v in variants:
            vk = v.get("VAR_KEY", "")
            chrom = v.get("CHROM", args.chrom)
            pos = int(v.get("POS", 0))
            ref = v.get("REF", ".")
            alt = v.get("ALT", ".")
            vtype_orig = v.get("VAR_TYPE", "")
            indel_len = v.get("INDEL_LEN", 0)
            abs_ilen = v.get("ABS_INDEL_LEN", 0)
            n_carriers = v.get("N_CARRIERS_ALL", 0)
            cf = v.get("CARRIER_FREQ_ALL", 0)
            fc = v.get("FREQ_CLASS", "unknown")

            # CSQ lookup
            csq_info = csq_lookup.get(vk, {"csq": ".", "gene": ".", "feature": "."})
            csq_consequence = csq_info["csq"]
            gene = csq_info["gene"]
            feature_context = csq_info["feature"]

            # Functional class fallback
            if feature_context == "." and vk in func_lookup:
                feature_context = func_lookup[vk].get("func_class", ".")
                if gene == ".":
                    gene = func_lookup[vk].get("gene_names", ".")

            # ESM
            esm_score = esm_lookup.get(vk, None)
            esm_str = f"{esm_score:.4f}" if esm_score is not None else "."

            # Variant class
            variant_class = assign_variant_class(vtype_orig, ref, alt, indel_len, csq_consequence)

            # Splice flag
            splice_flag = check_splice_proximity(csq_consequence)

            # Disruption severity
            if variant_class in ("DEL", "DUP", "INV", "INS", "BND", "TRA"):
                disruption_sev = assign_disruption_severity_sv(feature_context, variant_class,
                                                                int(abs_ilen) if abs_ilen not in ("", ".", "0") else 0)
            else:
                disruption_sev = assign_disruption_severity_small(csq_consequence, feature_context)

            # Evidence codes (annotation-based only)
            ev_codes = assign_evidence_codes(variant_class, disruption_sev,
                                              csq_consequence, splice_flag, feature_context)
            # ESM supportive
            if esm_score is not None and esm_score > 0.5:
                ev_codes.append("SU1")
            # Rarity supportive
            try:
                if int(n_carriers) <= 5 and int(n_carriers) > 0:
                    ev_codes.append("SU2")
            except (ValueError, TypeError):
                pass

            ev_str = ";".join(ev_codes) if ev_codes else "."

            fvar.write("\t".join(str(x) for x in [
                vk, chrom, pos, ref, alt,
                variant_class, vtype_orig,
                gene, feature_context, csq_consequence, splice_flag,
                esm_str,
                disruption_sev,
                ev_str,
                n_carriers, cf, fc,
                indel_len, abs_ilen,
            ]) + "\n")

            # Per-sample genotypes
            gts = gt_data.get(vk, {})
            for sample_id in gt_samples:
                gt = gts.get(sample_id, "./.")
                gt_clean = gt.replace("|", "/")
                alleles = gt_clean.split("/")
                is_missing = all(a in (".", "") for a in alleles)
                is_carrier = any(a not in ("0", ".", "") for a in alleles)

                if is_missing:
                    zyg = "missing"
                elif not is_carrier:
                    continue  # skip hom-ref to keep file manageable
                else:
                    non_ref = [a for a in alleles if a not in ("0", ".", "")]
                    if len(non_ref) == len(alleles):
                        zyg = "hom_alt"
                    else:
                        zyg = "het"

                # ROH check
                in_roh = "no"
                if sample_id in roh_regions:
                    for (rs, re) in roh_regions[sample_id]:
                        if rs <= pos <= re:
                            in_roh = "yes"
                            break

                fgt.write("\t".join(str(x) for x in [
                    vk, sample_id, gt, zyg, in_roh,
                ]) + "\n")

            n_annotated += 1
            if n_annotated % 50000 == 0:
                print(f"[VEF-20]   Annotated {n_annotated}/{len(variants)} …")

    print(f"[VEF-20] Annotated variants → {out_var_path}")
    print(f"[VEF-20] Per-sample genotypes → {out_gt_path}")
    print(f"[VEF-20] Total time: {_elapsed(t0)}")

if __name__ == "__main__":
    main()
