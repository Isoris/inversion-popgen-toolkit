#!/usr/bin/env python3
"""
21_score_breeding_concern.py — Variant Evidence Framework scoring engine

Takes annotated variants from 20_annotate_variant_consequences.py and assigns:
  - Final evidence codes (adds zygosity/ROH/recurrence-dependent codes)
  - Technical confidence tier (high / medium / low)
  - Zygosity-context tier (ZC0–ZC3)
  - Final breeding concern tier (BC0–BC4)

Then builds:
  - Per-variant scored table
  - Per-fish burden summary

Framework sections implemented:
  3. Evidence codes (SD, MD, SU)
  4. Confidence tiers
  5. Zygosity-context (ZC0–ZC3)
  6. Breeding concern (BC0–BC4)
  7. Combination rules
  8. Fish-level burden

Inputs:
  --annotated_variants   vef/annotated_variants.tsv (from 20)
  --per_sample_genotypes vef/per_sample_genotypes.tsv (from 20)
  --outdir               Output directory

Outputs:
  vef/scored_variants.tsv          — per variant × sample with BC score
  vef/variant_master_table.tsv     — per variant with worst-case BC
  vef/fish_burden_summary.tsv      — per fish aggregate burden
  vef/high_concern_report.tsv      — BC3 + BC4 variants only (for review)
"""

import os, sys, argparse, time
from collections import defaultdict, Counter

def _now(): return time.time()
def _elapsed(t0):
    dt = time.time() - t0
    return f"{dt:.1f}s" if dt < 60 else f"{dt/60:.1f}min"

# ── Evidence code classification ─────────────────────────────────────────────

STRONG_CODES = {"SD1", "SD2", "SD3", "SD4", "SD5", "SD6", "SD7"}
MODERATE_CODES = {"MD1", "MD2", "MD3", "MD4", "MD5", "MD6", "MD7"}
SUPPORTIVE_CODES = {"SU1", "SU2", "SU3", "SU4", "SU5", "SU6", "SU7"}

def count_evidence(codes):
    """Count strong, moderate, supportive evidence codes."""
    code_set = set(codes)
    n_strong = len(code_set & STRONG_CODES)
    n_moderate = len(code_set & MODERATE_CODES)
    n_supportive = len(code_set & SUPPORTIVE_CODES)
    return n_strong, n_moderate, n_supportive


# ── Zygosity-context tier ────────────────────────────────────────────────────

def assign_zygosity_context(zygosity, in_roh):
    """Assign ZC0–ZC3 (Section 5)."""
    is_hom = zygosity == "hom_alt"
    is_roh = in_roh == "yes"

    if is_hom and is_roh:
        return "ZC3"
    elif is_hom:
        return "ZC2"
    elif is_roh:  # het in ROH-adjacent
        return "ZC1"
    else:
        return "ZC0"


# ── Technical confidence ─────────────────────────────────────────────────────

def assign_confidence(freq_class, n_carriers, disruption_severity):
    """Assign confidence tier (Section 4).

    In the absence of per-variant quality metrics in this aggregated view,
    we use population-level proxies. When integrated with per-sample caller
    quality (Clair3 PASS/rescue, DELLY PRECISE, read support), this can
    be refined.
    """
    try:
        nc = int(n_carriers)
    except (ValueError, TypeError):
        nc = 0

    fc = str(freq_class).lower()

    # Common variants with many carriers = high confidence genotyping
    if nc >= 10:
        return "high"
    elif nc >= 3:
        return "medium"
    elif nc >= 1:
        # Singletons/doubletons: could be real rare or artifact
        if disruption_severity in ("very_high", "high"):
            return "medium"  # benefit of doubt for severe events
        return "low"
    else:
        return "low"


# ── Breeding concern scoring (Section 7 combination rules) ───────────────────

def score_breeding_concern(n_strong, n_moderate, n_supportive,
                            disruption_severity, zc_tier, confidence,
                            variant_class, feature_context):
    """
    Apply combination rules from Section 7.

    BC4 — Strong candidate for breeding caution
    BC3 — High breeding concern
    BC2 — Moderate breeding concern
    BC1 — Low / uncertain breeding concern
    BC0 — Minimal current concern
    """

    sev = str(disruption_severity).lower()
    conf = str(confidence).lower()
    zc = str(zc_tier)
    ctx = str(feature_context).lower()

    # ── BC4: Strong candidate for breeding caution ──
    # Rule: ≥1 strong + high confidence + ZC3
    if n_strong >= 1 and conf == "high" and zc == "ZC3":
        return "BC4"
    # Rule: ≥2 strong evidence lines
    if n_strong >= 2:
        return "BC4"
    # Rule: very high disruption SV affecting CDS/exon with strong support
    if sev == "very_high" and variant_class in ("DEL", "DUP", "INV", "INS") and \
       ("cds" in ctx or "exon" in ctx) and n_strong >= 1:
        return "BC4"
    # Rule: severe csq + hom + ROH
    if sev == "very_high" and zc == "ZC3":
        return "BC4"

    # ── BC3: High breeding concern ──
    # Rule: 1 strong + medium confidence
    if n_strong >= 1 and conf in ("high", "medium"):
        return "BC3"
    # Rule: 2 moderate
    if n_moderate >= 2:
        return "BC3"
    # Rule: 1 moderate + 2 supportive
    if n_moderate >= 1 and n_supportive >= 2:
        return "BC3"
    # Rule: moderate/high sev + homozygous
    if sev in ("very_high", "high") and zc in ("ZC2", "ZC3"):
        return "BC3"
    # Rule: moderate sev + hom + ROH
    if sev == "moderate" and zc == "ZC3":
        return "BC3"

    # ── BC2: Moderate breeding concern ──
    # Rule: 1 moderate
    if n_moderate >= 1:
        return "BC2"
    # Rule: ≥3 supportive
    if n_supportive >= 3:
        return "BC2"
    # Rule: moderate event with uncertain confidence
    if sev == "moderate" and conf == "low":
        return "BC2"
    # Rule: het but recurrent suspicious coding event
    if sev in ("moderate", "high") and n_supportive >= 1:
        return "BC2"

    # ── BC1: Low / uncertain breeding concern ──
    # Rule: only supportive evidence
    if n_supportive >= 1:
        return "BC1"
    # Rule: moderate annotation but weak confidence
    if sev == "moderate":
        return "BC1"
    # Rule: missense with no extra support
    if sev == "moderate" and n_moderate == 0 and n_strong == 0:
        return "BC1"

    # ── BC0: Minimal current concern ──
    return "BC0"


# ── Main ─────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--annotated_variants", required=True)
    p.add_argument("--per_sample_genotypes", required=True)
    p.add_argument("--outdir", required=True)
    return p.parse_args()

def main():
    args = parse_args()
    t0 = _now()
    outdir = os.path.join(args.outdir, "vef")
    os.makedirs(outdir, exist_ok=True)

    # ── Load annotated variants ──
    print("[VEF-21] Loading annotated variants …")
    variants = {}
    with open(args.annotated_variants) as f:
        header = f.readline().strip().split('\t')
        col = {h: i for i, h in enumerate(header)}
        for line in f:
            p = line.strip().split('\t')
            d = dict(zip(header, p))
            variants[d["VAR_KEY"]] = d
    print(f"[VEF-21]   {len(variants)} variants")

    # ── Load per-sample genotypes ──
    print("[VEF-21] Loading per-sample genotypes …")
    sample_gts = []  # list of dicts
    with open(args.per_sample_genotypes) as f:
        gh = f.readline().strip().split('\t')
        gc = {h: i for i, h in enumerate(gh)}
        for line in f:
            p = line.strip().split('\t')
            sample_gts.append(dict(zip(gh, p)))
    print(f"[VEF-21]   {len(sample_gts)} genotype entries")

    # ── Score each variant × sample ──
    print("[VEF-21] Scoring …")

    scored_path = os.path.join(outdir, "scored_variants.tsv")
    scored_cols = [
        "VAR_KEY", "SAMPLE", "CHROM", "POS", "REF", "ALT",
        "VARIANT_CLASS", "GENE", "FEATURE_CONTEXT",
        "CSQ_CONSEQUENCE", "SPLICE_FLAG", "ESM_SCORE",
        "DISRUPTION_SEVERITY",
        "GT", "ZYGOSITY", "IN_ROH", "ZC_TIER",
        "EVIDENCE_CODES", "N_STRONG", "N_MODERATE", "N_SUPPORTIVE",
        "TECH_CONFIDENCE",
        "FINAL_BREEDING_CONCERN",
        "N_CARRIERS_ALL", "FREQ_CLASS",
    ]

    # Track per-fish burden
    fish_burden = defaultdict(lambda: {
        "N_BC4": 0, "N_BC3": 0, "N_BC2": 0, "N_BC1": 0, "N_BC0": 0,
        "N_strong_disruption": 0, "N_hom_high": 0, "N_high_in_ROH": 0,
        "N_gene_disruptive_SV": 0, "N_splice_flags": 0,
        "N_high_ESM_missense": 0,
        "burden_score": 0.0,
        "burden_in_ROH": 0.0,
        "burden_small_var": 0.0,
        "burden_SV": 0.0,
        "burden_splice": 0.0,
    })

    # Track worst BC per variant
    var_worst_bc = {}

    n_scored = 0
    with open(scored_path, 'w') as fout:
        fout.write("\t".join(scored_cols) + "\n")

        for sg in sample_gts:
            vk = sg["VAR_KEY"]
            sample = sg["SAMPLE"]
            gt = sg["GT"]
            zyg = sg["ZYGOSITY"]
            in_roh = sg["IN_ROH"]

            if zyg == "missing":
                continue

            v = variants.get(vk)
            if v is None:
                continue

            # Get base evidence codes from annotation
            base_codes = v.get("EVIDENCE_CODES", "").split(";")
            base_codes = [c for c in base_codes if c and c != "."]

            # Add zygosity/ROH-dependent codes
            ev_codes = list(base_codes)

            sev = v.get("DISRUPTION_SEVERITY", "minimal")
            vclass = v.get("VARIANT_CLASS", "unknown")
            splice = v.get("SPLICE_FLAG", "none")

            # SD6: Homozygous very-high-disruption inside ROH
            if zyg == "hom_alt" and in_roh == "yes" and sev == "very_high":
                ev_codes.append("SD6")

            # MD5: Homozygous moderate-disruption
            if zyg == "hom_alt" and sev in ("moderate", "high"):
                ev_codes.append("MD5")

            # MD6: Variant lies inside ROH
            if in_roh == "yes":
                ev_codes.append("MD6")

            # SU7: ROH context
            if in_roh == "yes" and "MD6" not in ev_codes:
                ev_codes.append("SU7")

            # Deduplicate
            ev_codes = sorted(set(ev_codes))

            # Count evidence
            n_s, n_m, n_su = count_evidence(ev_codes)

            # Zygosity context
            zc = assign_zygosity_context(zyg, in_roh)

            # Confidence
            conf = assign_confidence(v.get("FREQ_CLASS", ""), v.get("N_CARRIERS_ALL", 0), sev)

            # Breeding concern
            bc = score_breeding_concern(
                n_s, n_m, n_su, sev, zc, conf,
                vclass, v.get("FEATURE_CONTEXT", "")
            )

            fout.write("\t".join(str(x) for x in [
                vk, sample, v["CHROM"], v["POS"], v["REF"], v["ALT"],
                vclass, v.get("GENE", "."), v.get("FEATURE_CONTEXT", "."),
                v.get("CSQ_CONSEQUENCE", "."), splice, v.get("ESM_SCORE", "."),
                sev,
                gt, zyg, in_roh, zc,
                ";".join(ev_codes) if ev_codes else ".",
                n_s, n_m, n_su,
                conf,
                bc,
                v.get("N_CARRIERS_ALL", 0), v.get("FREQ_CLASS", ""),
            ]) + "\n")

            # ── Update fish burden ──
            fb = fish_burden[sample]
            fb[f"N_{bc}"] += 1

            bc_weight = {"BC4": 8, "BC3": 5, "BC2": 2, "BC1": 0.5, "BC0": 0}
            w = bc_weight.get(bc, 0)
            fb["burden_score"] += w

            if n_s > 0:
                fb["N_strong_disruption"] += 1
            if zyg == "hom_alt" and sev in ("high", "very_high"):
                fb["N_hom_high"] += 1
            if in_roh == "yes" and bc in ("BC3", "BC4"):
                fb["N_high_in_ROH"] += 1
                fb["burden_in_ROH"] += w
            if vclass in ("DEL", "DUP", "INV", "INS") and "exon" in str(v.get("FEATURE_CONTEXT", "")).lower():
                fb["N_gene_disruptive_SV"] += 1
                fb["burden_SV"] += w
            elif vclass in ("DEL", "DUP", "INV", "INS"):
                fb["burden_SV"] += w
            else:
                fb["burden_small_var"] += w
            if splice != "none":
                fb["N_splice_flags"] += 1
                fb["burden_splice"] += w
            esm_val = v.get("ESM_SCORE", ".")
            if esm_val not in (".", "NA", "") and float(esm_val) > 0.7:
                fb["N_high_ESM_missense"] += 1

            # Track worst BC per variant
            bc_rank = {"BC4": 4, "BC3": 3, "BC2": 2, "BC1": 1, "BC0": 0}
            if vk not in var_worst_bc or bc_rank.get(bc, 0) > bc_rank.get(var_worst_bc[vk], 0):
                var_worst_bc[vk] = bc

            n_scored += 1
            if n_scored % 100000 == 0:
                print(f"[VEF-21]   Scored {n_scored} entries …")

    print(f"[VEF-21] Scored variants → {scored_path}  ({n_scored} entries)")

    # ── Variant master table ──
    master_path = os.path.join(outdir, "variant_master_table.tsv")
    with open(master_path, 'w') as fout:
        mcols = [
            "VAR_KEY", "CHROM", "POS", "REF", "ALT",
            "VARIANT_CLASS", "GENE", "FEATURE_CONTEXT",
            "CSQ_CONSEQUENCE", "SPLICE_FLAG", "ESM_SCORE",
            "DISRUPTION_SEVERITY", "EVIDENCE_CODES",
            "N_CARRIERS_ALL", "FREQ_CLASS",
            "WORST_BC",
        ]
        fout.write("\t".join(mcols) + "\n")
        for vk, v in variants.items():
            worst = var_worst_bc.get(vk, "BC0")
            fout.write("\t".join(str(x) for x in [
                vk, v["CHROM"], v["POS"], v["REF"], v["ALT"],
                v.get("VARIANT_CLASS", "."), v.get("GENE", "."), v.get("FEATURE_CONTEXT", "."),
                v.get("CSQ_CONSEQUENCE", "."), v.get("SPLICE_FLAG", "."), v.get("ESM_SCORE", "."),
                v.get("DISRUPTION_SEVERITY", "."), v.get("EVIDENCE_CODES", "."),
                v.get("N_CARRIERS_ALL", 0), v.get("FREQ_CLASS", ""),
                worst,
            ]) + "\n")
    print(f"[VEF-21] Variant master → {master_path}")

    # ── Fish burden summary ──
    burden_path = os.path.join(outdir, "fish_burden_summary.tsv")
    burden_cols = [
        "SAMPLE",
        "N_BC4", "N_BC3", "N_BC2", "N_BC1", "N_BC0",
        "N_strong_disruption", "N_hom_high", "N_high_in_ROH",
        "N_gene_disruptive_SV", "N_splice_flags", "N_high_ESM_missense",
        "burden_score", "burden_in_ROH", "burden_small_var", "burden_SV", "burden_splice",
    ]
    with open(burden_path, 'w') as fout:
        fout.write("\t".join(burden_cols) + "\n")
        for sample in sorted(fish_burden.keys()):
            fb = fish_burden[sample]
            fout.write("\t".join(str(x) for x in [
                sample,
                fb["N_BC4"], fb["N_BC3"], fb["N_BC2"], fb["N_BC1"], fb["N_BC0"],
                fb["N_strong_disruption"], fb["N_hom_high"], fb["N_high_in_ROH"],
                fb["N_gene_disruptive_SV"], fb["N_splice_flags"], fb["N_high_ESM_missense"],
                round(fb["burden_score"], 1), round(fb["burden_in_ROH"], 1),
                round(fb["burden_small_var"], 1), round(fb["burden_SV"], 1),
                round(fb["burden_splice"], 1),
            ]) + "\n")
    print(f"[VEF-21] Fish burden → {burden_path}")

    # ── High concern report (BC3 + BC4 only) ──
    high_path = os.path.join(outdir, "high_concern_report.tsv")
    n_high = 0
    with open(scored_path) as fin, open(high_path, 'w') as fout:
        header_line = fin.readline()
        fout.write(header_line)
        for line in fin:
            if "\tBC4\t" in line or "\tBC3\t" in line:
                fout.write(line)
                n_high += 1
    print(f"[VEF-21] High concern report: {n_high} entries → {high_path}")

    # ── Summary stats ──
    print(f"\n[VEF-21] ======= SCORING SUMMARY =======")
    bc_total = Counter()
    for vk, bc in var_worst_bc.items():
        bc_total[bc] += 1
    for bc in ["BC4", "BC3", "BC2", "BC1", "BC0"]:
        print(f"  {bc}: {bc_total.get(bc, 0)} variants")

    burden_scores = [fb["burden_score"] for fb in fish_burden.values()]
    if burden_scores:
        import statistics
        print(f"\n  Fish burden scores:")
        print(f"    median: {statistics.median(burden_scores):.1f}")
        print(f"    mean:   {statistics.mean(burden_scores):.1f}")
        print(f"    max:    {max(burden_scores):.1f}")
        print(f"    min:    {min(burden_scores):.1f}")

    print(f"\n[VEF-21] Total time: {_elapsed(t0)}")

if __name__ == "__main__":
    main()
