#!/usr/bin/env bash
#SBATCH --job-name=CONS_14_score
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=slurm_logs/14_score_%j.out
#SBATCH --error=slurm_logs/14_score_%j.err
set -euo pipefail

# =============================================================================
# STEP 14 — Variant Master Table + Scoring (CORE-ONLY)
# =============================================================================
# Merges annotation layers into one master variant table, then scores:
#   A. Intrinsic deleterious evidence: consequence + SIFT4G + VESM + splice
#   B. Hatchery risk (computed in STEP 15 with genotype data)
#   C. Final priority → Class A/B/C/D
#
# NO GERP. NO phastCons. NO orthology tiers. NO CAFE. NO DupGen.
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config/pipeline.config.sh"

mkdir -p "${DIR_MERGED}" "${DIR_LOGS}"

echo "=== STEP 14: Building variant master table + scoring (CORE-ONLY) ==="

python3 - <<'PYEOF'
"""
Build variant_master_scored.tsv by merging:
- SnpEff consequences (canonical per variant)
- bcftools csq consequences (haplotype-aware, complementary)
- SIFT4G scores (missense only)
- VESM scores (protein language model LLR, from 05B)
- Splice module classifications (from Clair3 splice annotation)

Then score each variant:
- deleterious_evidence_score = csq + SIFT + VESM + splice_bonus
- hatchery_risk_score = 0 (placeholder — STEP 15 fills this)
- final_priority_score → priority_class (A/B/C/D)
"""
import os, sys, gzip, glob
from collections import defaultdict

base = os.environ['MODCONS']
merged_dir = os.path.join(base, '16_merged_variant_tables')
snpeff_dir = os.path.join(base, '04_snpeff')
sift_dir = os.path.join(base, '05_sift4g')
vesm_dir = os.path.join(base, '05B_vesm')
var_dir = os.path.join(base, '03_variants', 'normalized')

# Splice module cohort results
clair3_base = os.environ.get('CLAIR3_DIR', '')
splice_base = os.path.join(clair3_base, 'splice_annotation', 'results') if clair3_base else ''

# Priority thresholds (from config, with defaults)
THRESH_A = int(os.environ.get('PRIORITY_THRESHOLD_A', '16'))
THRESH_B = int(os.environ.get('PRIORITY_THRESHOLD_B', '10'))
THRESH_C = int(os.environ.get('PRIORITY_THRESHOLD_C', '5'))

# ─── Load scoring weights ───
score_config = os.path.join(base, '00_config', 'scoring_weights.tsv')
weights = {}
if os.path.exists(score_config):
    with open(score_config) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                cat, comp, cond, pts = parts[0], parts[1], parts[2], parts[3]
                try:
                    weights[(cat, comp)] = float(pts)
                except ValueError:
                    pass
    print(f"  Loaded {len(weights)} scoring weights")

# ─── Helpers ───

def severity_rank(impact):
    return {'HIGH': 4, 'MODERATE': 3, 'LOW': 2, 'MODIFIER': 1}.get(impact, 0)

def safe_float(val, default=None):
    try:
        return float(val) if val and val != '' else default
    except (ValueError, TypeError):
        return default

def score_consequence(annotation):
    """Score based on SnpEff annotation string."""
    annotation = annotation.lower() if annotation else ''
    if any(x in annotation for x in ['stop_gained', 'frameshift']):
        return weights.get(('consequence', 'stop_gained'), 8)
    elif any(x in annotation for x in ['splice_acceptor', 'splice_donor']):
        return weights.get(('consequence', 'splice_acceptor'), 7)
    elif 'start_lost' in annotation:
        return weights.get(('consequence', 'start_lost'), 6)
    elif 'stop_lost' in annotation:
        return weights.get(('consequence', 'stop_lost'), 5)
    elif 'missense' in annotation:
        return weights.get(('consequence', 'missense'), 4)
    elif 'inframe' in annotation:
        return weights.get(('consequence', 'inframe_insertion'), 2)
    elif 'splice_region' in annotation:
        return weights.get(('consequence', 'splice_region'), 1)
    return 0

def score_sift(sift_score, sift_class):
    s = safe_float(sift_score)
    if s is not None and s < 0.05:
        return weights.get(('sift', 'deleterious'), 4)
    elif s is not None and s < 0.1:
        return weights.get(('sift', 'borderline'), 2)
    return 0

def score_vesm(vesm_llr, vesm_class):
    llr = safe_float(vesm_llr)
    if llr is not None and llr < -7:
        return weights.get(('vesm', 'likely_damaging'), 5)
    elif llr is not None and llr < -3:
        return weights.get(('vesm', 'possibly_damaging'), 2)
    return 0

def score_splice_module(splice_info):
    """Score fine-grained splice effects that generic SnpEff underweights.
    Returns ADDITIONAL points beyond what score_consequence() gives.
    Capped at 5 to avoid splice_region + bonus outscooring frameshift."""
    if not splice_info:
        return 0
    effect = splice_info.get('EFFECT', '').lower()
    subclass = splice_info.get('SPLICE_SUBCLASS', '')
    priority = splice_info.get('PRIORITY_CLASS', splice_info.get('BEST_PRIORITY', 'D'))
    bonus = 0
    if 'splice_donor_5th_base' in effect:
        bonus = max(bonus, weights.get(('splice', 'donor_5th_base'), 3))
    if 'splice_branch' in effect:
        bonus = max(bonus, weights.get(('splice', 'branch_point'), 3))
    if 'splice_polypyrimidine' in effect:
        bonus = max(bonus, weights.get(('splice', 'polypyrimidine'), 2))
    if subclass in ('splice_donor', 'splice_acceptor') and priority == 'A':
        bonus = max(bonus, weights.get(('splice', 'class_A'), 4))
    elif subclass == 'splice_region' and priority in ('A', 'B'):
        bonus = max(bonus, weights.get(('splice', 'class_B'), 2))
    return min(bonus, 5)

def classify_priority(total_score):
    if total_score >= THRESH_A: return 'A'
    elif total_score >= THRESH_B: return 'B'
    elif total_score >= THRESH_C: return 'C'
    return 'D'

def build_reason(csq_pts, sift_pts, vesm_pts, splice_pts):
    parts = []
    if csq_pts >= 7: parts.append('severe_consequence')
    elif csq_pts >= 4: parts.append('protein_changing')
    if sift_pts >= 4: parts.append('SIFT_damaging')
    elif sift_pts >= 2: parts.append('SIFT_borderline')
    if vesm_pts >= 5: parts.append('VESM_damaging')
    elif vesm_pts >= 2: parts.append('VESM_suspicious')
    if splice_pts >= 3: parts.append('splice_module_strong')
    elif splice_pts >= 1: parts.append('splice_module_hit')
    return '; '.join(parts) if parts else 'low_evidence'

# ─── Load SnpEff canonical consequences ───
snpeff_canon = os.path.join(snpeff_dir, 'snpeff_canonical_consequences.tsv')
snpeff_data = {}

if os.path.exists(snpeff_canon):
    with open(snpeff_canon) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            d = dict(zip(header, parts))
            vk = d.get('var_key', '')
            if vk:
                if vk not in snpeff_data or \
                   severity_rank(d.get('impact','')) > severity_rank(snpeff_data[vk].get('impact','')):
                    snpeff_data[vk] = d
    print(f"  SnpEff annotations: {len(snpeff_data)} variants")
else:
    print(f"  WARNING: SnpEff canonical consequences not found at {snpeff_canon}")

# ─── Load bcftools csq consequences (complementary) ───
csq_path = os.path.join(snpeff_dir, 'bcftools_csq', 'bcftools_csq_consequences.tsv')
csq_data = {}
if os.path.exists(csq_path):
    with open(csq_path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            d = dict(zip(header, parts))
            vk = d.get('var_key', '')
            if vk and vk not in csq_data:
                csq_data[vk] = d
    print(f"  bcftools csq annotations: {len(csq_data)} variants")

# ─── Load SIFT4G scores ───
sift_path = os.path.join(sift_dir, 'sift4g_scores.tsv')
sift_data = {}
if os.path.exists(sift_path):
    with open(sift_path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            d = dict(zip(header, parts))
            vk = d.get('var_key', '')
            if vk:
                sift_data[vk] = d
    print(f"  SIFT4G scores: {len(sift_data)} variants")
else:
    print(f"  WARNING: SIFT4G scores not found — run STEP_05 first")

# ─── Load VESM scores ───
vesm_path = os.path.join(vesm_dir, 'vesm_variant_scores.tsv')
vesm_data = {}
if os.path.exists(vesm_path):
    with open(vesm_path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            d = dict(zip(header, parts))
            vk = d.get('var_key', '')
            if vk:
                vesm_data[vk] = d
    print(f"  VESM scores: {len(vesm_data)} variants")
else:
    print(f"  WARNING: VESM scores not found — run STEP_05B first (GPU)")

# ─── Load splice module cohort results ───
splice_data = {}
if splice_base and os.path.isdir(splice_base):
    splice_files = glob.glob(os.path.join(splice_base, '*', '_cohort', 'cohort_splice_by_variant.tsv'))
    for sf in splice_files:
        try:
            with open(sf) as f:
                sh = f.readline().strip().split('\t')
                for line in f:
                    parts = line.strip().split('\t')
                    d = dict(zip(sh, parts))
                    vk = d.get('VARIANT_ID', '')
                    if vk:
                        splice_data[vk] = d
        except Exception:
            pass
    print(f"  Splice module: {len(splice_data)} variants from {len(splice_files)} chrom files")

    # Fallback: per-sample splice_summary
    if not splice_data:
        for sf in glob.glob(os.path.join(splice_base, '*', '*', 'clair3.splice_summary.tsv')):
            try:
                with open(sf) as f:
                    sh = f.readline().strip().split('\t')
                    for line in f:
                        parts = line.strip().split('\t')
                        d = dict(zip(sh, parts))
                        vk = d.get('VARIANT_ID', '')
                        if vk and vk not in splice_data:
                            splice_data[vk] = d
            except Exception:
                pass
        if splice_data:
            print(f"  Splice module (per-sample fallback): {len(splice_data)} variants")
else:
    print(f"  Splice module: not found — continuing without splice bonus")

# ─── Build variant master + score ───
print("\n  Building variant master table...")

# Collect all variant keys from SnpEff + bcftools csq
all_var_keys = set(snpeff_data.keys()) | set(csq_data.keys())
print(f"  Total variants to score: {len(all_var_keys)}")

out_path = os.path.join(merged_dir, 'variant_master_scored.tsv')
out_header = [
    'var_key', 'chr', 'pos', 'ref', 'alt',
    # SnpEff
    'snpeff_annotation', 'snpeff_impact', 'gene_id', 'gene_name',
    'transcript_id', 'hgvs_c', 'hgvs_p',
    # bcftools csq
    'bcsq_consequence', 'bcsq_gene', 'bcsq_transcript',
    # SIFT
    'sift_score', 'sift_class',
    # VESM
    'vesm_llr', 'vesm_class',
    # Splice module
    'splice_subclass', 'splice_priority', 'splice_direct_hit', 'splice_effect_detail',
    # Score components (visible for auditability)
    'csq_score', 'sift_score_pts', 'vesm_score_pts', 'splice_score_pts',
    'deleterious_evidence_score',
    'hatchery_risk_score', 'final_priority_score',
    'priority_class', 'priority_reason'
]

priority_counts = defaultdict(int)

with open(out_path, 'w') as out:
    out.write('\t'.join(out_header) + '\n')

    for vk in sorted(all_var_keys):
        se = snpeff_data.get(vk, {})
        cq = csq_data.get(vk, {})
        si = sift_data.get(vk, {})

        # Parse var_key
        vk_parts = vk.split(':')
        chrom = vk_parts[0] if len(vk_parts) > 0 else ''
        pos = vk_parts[1] if len(vk_parts) > 1 else ''
        ref = vk_parts[2] if len(vk_parts) > 2 else ''
        alt = vk_parts[3] if len(vk_parts) > 3 else ''

        # SnpEff fields
        annotation = se.get('annotation', '')
        impact = se.get('impact', '')
        gene_id = se.get('gene_id', '')
        gene_name = se.get('gene_name', '')
        tx_id = se.get('feature_id', '')
        hgvs_c = se.get('hgvs_c', '')
        hgvs_p = se.get('hgvs_p', '')

        # bcftools csq fields
        bcsq_csq = cq.get('consequence', '')
        bcsq_gene = cq.get('gene', '')
        bcsq_tx = cq.get('transcript', '')

        # If no SnpEff gene, try bcftools csq gene
        if not gene_id and bcsq_gene:
            gene_id = bcsq_gene

        # SIFT
        sift_sc = si.get('sift_score', '')
        sift_cl = si.get('sift_class', '')

        # VESM
        vm = vesm_data.get(vk, {})
        vesm_llr = vm.get('vesm_llr', '')
        vesm_cl = vm.get('vesm_class', '')

        # Splice module
        sp = splice_data.get(vk, {})
        splice_sub = sp.get('SPLICE_SUBCLASS', '')
        splice_pri = sp.get('BEST_PRIORITY', sp.get('PRIORITY_CLASS', ''))
        splice_direct = sp.get('DIRECT_HIT', sp.get('DIRECT_SPLICE_HIT', ''))
        splice_eff = sp.get('EFFECT', '')

        # --- SCORING ---
        csq_pts = score_consequence(annotation)

        # If bcftools csq gives a MORE severe consequence, use that
        bcsq_pts = score_consequence(bcsq_csq)
        if bcsq_pts > csq_pts:
            csq_pts = bcsq_pts

        sift_pts = score_sift(sift_sc, sift_cl)
        vesm_pts = score_vesm(vesm_llr, vesm_cl)
        splice_pts = score_splice_module(sp)

        del_evidence = csq_pts + sift_pts + vesm_pts + splice_pts

        # Hatchery risk: placeholder — STEP 15 fills with actual genotype data
        hatch_risk = 0

        final_score = del_evidence + hatch_risk
        pclass = classify_priority(final_score)
        reason = build_reason(csq_pts, sift_pts, vesm_pts, splice_pts)

        priority_counts[pclass] += 1

        row = [vk, chrom, pos, ref, alt,
               annotation, impact, gene_id, gene_name,
               tx_id, hgvs_c, hgvs_p,
               bcsq_csq, bcsq_gene, bcsq_tx,
               sift_sc, sift_cl,
               vesm_llr, vesm_cl,
               splice_sub, splice_pri, splice_direct, splice_eff,
               str(csq_pts), str(sift_pts), str(vesm_pts), str(splice_pts),
               f"{del_evidence:.1f}",
               f"{hatch_risk:.1f}", f"{final_score:.1f}",
               pclass, reason]
        out.write('\t'.join(str(x) for x in row) + '\n')

print(f"\n  Priority class distribution (deleterious evidence only, before hatchery risk):")
for c in ['A', 'B', 'C', 'D']:
    print(f"    Class {c}: {priority_counts.get(c, 0)}")
print(f"\n  NOTE: These classes will shift upward after STEP 15 adds hatchery risk scores.")
print(f"\n  Output: {out_path}")
PYEOF

echo "=== STEP 14 complete ==="
