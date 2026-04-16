#!/usr/bin/env bash
#SBATCH --job-name=CONS_15_burden
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=slurm_logs/15_burden_%j.out
#SBATCH --error=slurm_logs/15_burden_%j.err
set -euo pipefail

# =============================================================================
# STEP 15 — Burden Tables (CORE-ONLY)
# =============================================================================
# Reads variant_master_scored.tsv + genotype data from VCFs + ROH BEDs
# Produces:
#   A. Hatchery risk scores → updates variant_master_final.tsv
#   B. deleterious_variant_genotype_matrix.tsv
#   C. deleterious_burden_per_individual.tsv
#   D. deleterious_burden_per_ROH.tsv
#   E. deleterious_burden_per_gene.tsv (basic, no orthology/CAFE)
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config/pipeline.config.sh"

mkdir -p "${DIR_BURDEN}" "${DIR_ROH}" "${DIR_LOGS}"

echo "=== STEP 15: Building burden tables (CORE-ONLY) ==="

python3 - <<'PYEOF'
"""
Build all burden tables from scored variants + genotypes + ROH overlaps.
CORE-ONLY: no gene_context_master, no orthology tiers, no CAFE.
"""
import os, sys, gzip, glob
from collections import defaultdict

base = os.environ['MODCONS']
merged_dir = os.path.join(base, '16_merged_variant_tables')
burden_dir = os.path.join(base, '17_burden_tables')
roh_overlap_dir = os.path.join(base, '18_roh_overlap')
var_dir = os.path.join(base, '03_variants', 'normalized')
roh_dir = os.environ.get('ROH_DIR', '')
samples_all = os.environ.get('SAMPLES_ALL', '')
samples_unrel = os.environ.get('SAMPLES_UNRELATED', '')

# Priority thresholds
THRESH_A = int(os.environ.get('PRIORITY_THRESHOLD_A', '16'))
THRESH_B = int(os.environ.get('PRIORITY_THRESHOLD_B', '10'))
THRESH_C = int(os.environ.get('PRIORITY_THRESHOLD_C', '5'))
MIN_DP = int(os.environ.get('MIN_GT_DP', '3'))

os.makedirs(burden_dir, exist_ok=True)
os.makedirs(roh_overlap_dir, exist_ok=True)

# ─── Load sample lists ───
all_samples = []
if os.path.exists(samples_all):
    with open(samples_all) as f:
        all_samples = [l.strip() for l in f if l.strip()]
print(f"  Samples: {len(all_samples)}")

unrel_samples = set()
if os.path.exists(samples_unrel):
    with open(samples_unrel) as f:
        unrel_samples = {l.strip() for l in f if l.strip()}

# ─── Load scored variants ───
var_master_path = os.path.join(merged_dir, 'variant_master_scored.tsv')
variants = {}
if os.path.exists(var_master_path):
    with open(var_master_path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            d = dict(zip(header, parts))
            variants[d['var_key']] = d
    print(f"  Scored variants: {len(variants)}")
else:
    print("  ERROR: variant_master_scored.tsv not found — run STEP 14")
    sys.exit(1)

# ─── Load ROH segments per sample ───
roh_segments = defaultdict(list)
if roh_dir and os.path.exists(roh_dir):
    for roh_file in glob.glob(os.path.join(roh_dir, '**', '*.roh*.bed'), recursive=True) + \
                    glob.glob(os.path.join(roh_dir, '**', '*ROH*'), recursive=True):
        try:
            with open(roh_file) as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 4:
                        sample = parts[0]
                        chrom, start, end = parts[1], int(parts[2]), int(parts[3])
                        roh_segments[sample].append((chrom, start, end))
                    elif len(parts) == 3:
                        sample = os.path.basename(roh_file).split('.')[0]
                        chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                        roh_segments[sample].append((chrom, start, end))
        except (ValueError, IndexError):
            continue
    print(f"  ROH segments loaded for {len(roh_segments)} samples")
else:
    print(f"  WARNING: ROH directory not found — ROH overlap columns will be zero")

# ─── Extract genotypes from VCFs ───
print("  Extracting genotypes from VCFs...")

target_vars = set(variants.keys())
chr_vars = defaultdict(set)
for vk in target_vars:
    parts = vk.split(':')
    if len(parts) >= 2:
        chr_vars[parts[0]].add(vk)

genotypes = defaultdict(dict)  # var_key → {sample: GT}
gt_depths = defaultdict(dict)  # var_key → {sample: DP}
chroms = [f"C_gar_LG{i:02d}" for i in range(1, 29)]

for chrom in chroms:
    vcf_path = os.path.join(var_dir, f'{chrom}.clair3.norm.vcf.gz')
    if not os.path.exists(vcf_path) or chrom not in chr_vars:
        continue

    needed = chr_vars[chrom]
    with gzip.open(vcf_path, 'rt') as f:
        vcf_samples = []
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                parts = line.strip().split('\t')
                vcf_samples = parts[9:]
                continue
            parts = line.strip().split('\t')
            if len(parts) < 10:
                continue
            c, pos, _, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            vk = f"{c}:{pos}:{ref}:{alt}"
            if vk not in needed:
                continue

            fmt = parts[8].split(':')
            gt_idx = fmt.index('GT') if 'GT' in fmt else 0
            dp_idx = fmt.index('DP') if 'DP' in fmt else -1

            for i, sample_data in enumerate(parts[9:]):
                if i >= len(vcf_samples):
                    break
                sample = vcf_samples[i]
                fields = sample_data.split(':')
                gt = fields[gt_idx] if gt_idx < len(fields) else './.'

                # Extract depth
                dp = 0
                if dp_idx >= 0 and dp_idx < len(fields):
                    try:
                        dp = int(fields[dp_idx])
                    except (ValueError, TypeError):
                        dp = 0
                gt_depths[vk][sample] = dp

                # Apply minimum depth filter
                if dp < MIN_DP:
                    genotypes[vk][sample] = -1
                    continue

                if gt in ('./.', '.'):
                    genotypes[vk][sample] = -1
                elif gt in ('0/0', '0|0'):
                    genotypes[vk][sample] = 0
                elif gt in ('0/1', '1/0', '0|1', '1|0'):
                    genotypes[vk][sample] = 1
                elif gt in ('1/1', '1|1'):
                    genotypes[vk][sample] = 2
                else:
                    genotypes[vk][sample] = -1

print(f"  Genotypes extracted for {len(genotypes)} variants (DP >= {MIN_DP})")

# ─── ROH helpers ───
def in_roh(sample, chrom, pos):
    for c, s, e in roh_segments.get(sample, []):
        if c == chrom and s <= pos <= e:
            return True
    return False

def in_long_roh(sample, chrom, pos, min_len=1_000_000):
    for c, s, e in roh_segments.get(sample, []):
        if c == chrom and s <= pos <= e and (e - s) >= min_len:
            return True
    return False

# ═══ A. Hatchery risk scores ═══
print("\n  Computing hatchery risk scores...")

for vk, vdata in variants.items():
    gt = genotypes.get(vk, {})
    n_hom = sum(1 for v in gt.values() if v == 2)
    n_het = sum(1 for v in gt.values() if v == 1)
    n_called = sum(1 for v in gt.values() if v >= 0)
    maf = (n_het + 2*n_hom) / (2*n_called) if n_called > 0 else 0

    chrom = vdata.get('chr', '')
    pos = int(vdata.get('pos', 0)) if vdata.get('pos', '').isdigit() else 0

    n_hom_in_roh = sum(1 for s in all_samples if gt.get(s, -1) == 2 and in_roh(s, chrom, pos))
    n_hom_in_long_roh = sum(1 for s in all_samples if gt.get(s, -1) == 2 and in_long_roh(s, chrom, pos))

    hr = 0
    if maf < 0.01: hr += 3
    elif maf < 0.05: hr += 2
    elif maf < 0.20: hr += 1

    if n_hom >= 10: hr += 3
    elif n_hom >= 3: hr += 2
    elif n_hom >= 1: hr += 1

    if n_hom_in_long_roh > 0: hr += 3
    elif n_hom_in_roh > 0: hr += 2

    vdata['hatchery_risk_score'] = str(hr)
    vdata['n_hom_alt'] = str(n_hom)
    vdata['n_het'] = str(n_het)
    vdata['n_called'] = str(n_called)
    vdata['maf'] = f"{maf:.4f}"
    vdata['n_hom_in_roh'] = str(n_hom_in_roh)
    vdata['n_hom_in_long_roh'] = str(n_hom_in_long_roh)

    try:
        del_ev = float(vdata.get('deleterious_evidence_score', '0'))
    except ValueError:
        del_ev = 0
    final = del_ev + hr
    if final >= THRESH_A: pc = 'A'
    elif final >= THRESH_B: pc = 'B'
    elif final >= THRESH_C: pc = 'C'
    else: pc = 'D'
    vdata['final_priority_score'] = f"{final:.1f}"
    vdata['priority_class'] = pc

# Write updated variant master
updated_path = os.path.join(merged_dir, 'variant_master_final.tsv')
extra_cols = ['n_hom_alt', 'n_het', 'n_called', 'maf', 'n_hom_in_roh', 'n_hom_in_long_roh']
with open(var_master_path) as f:
    orig_header = f.readline().strip().split('\t')
full_header = orig_header + [c for c in extra_cols if c not in orig_header]

with open(updated_path, 'w') as out:
    out.write('\t'.join(full_header) + '\n')
    for vk in sorted(variants.keys()):
        row = [variants[vk].get(h, '') for h in full_header]
        out.write('\t'.join(row) + '\n')

# Count final classes
final_counts = defaultdict(int)
for v in variants.values():
    final_counts[v.get('priority_class', 'D')] += 1
print(f"\n  Final priority class distribution (with hatchery risk):")
for c in ['A', 'B', 'C', 'D']:
    print(f"    Class {c}: {final_counts.get(c, 0)}")
print(f"  Updated variant master: {updated_path}")

# ═══ B. Genotype matrix ═══
print("\n  Building genotype matrix...")
gt_matrix_path = os.path.join(burden_dir, 'deleterious_variant_genotype_matrix.tsv')

candidate_vks = [vk for vk, vd in variants.items() if vd.get('priority_class') in ('A', 'B')]
print(f"  Class A/B variants: {len(candidate_vks)}")

with open(gt_matrix_path, 'w') as out:
    out.write('var_key\tchr\tpos\tref\talt\tpriority_class\t' + '\t'.join(all_samples) + '\n')
    for vk in sorted(candidate_vks):
        vd = variants[vk]
        gt = genotypes.get(vk, {})
        vals = [str(gt.get(s, -1)) for s in all_samples]
        # -1 → NA for downstream R
        vals_clean = ['NA' if v == '-1' else v for v in vals]
        out.write(f"{vk}\t{vd.get('chr','')}\t{vd.get('pos','')}\t"
                  f"{vd.get('ref','')}\t{vd.get('alt','')}\t"
                  f"{vd.get('priority_class','')}\t"
                  + '\t'.join(vals_clean) + '\n')

# ═══ C. Per-individual burden ═══
print("\n  Computing per-individual burden...")
burden_path = os.path.join(burden_dir, 'deleterious_burden_per_individual.tsv')

genome_size = 0
fai_path = os.environ.get('REF_FAI', '')
if os.path.exists(fai_path):
    with open(fai_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if parts[0].startswith('C_gar_LG'):
                genome_size += int(parts[1])
if genome_size == 0:
    genome_size = 800_000_000

with open(burden_path, 'w') as out:
    out.write('sample_id\tis_unrelated\tn_del_total\tn_del_het\tn_del_hom\t'
              'sum_del_evidence\tsum_hatchery_risk\tsum_final_priority\t'
              'n_classA\tn_classB\tn_classA_hom\tn_classA_in_ROH\t'
              'n_classB_in_ROH\tn_del_in_ROH\tn_del_hom_in_ROH\tn_del_in_long_ROH\t'
              'total_ROH_bp\tn_ROH_segments\tFROH\n')

    for sample in all_samples:
        is_unrel = 'yes' if sample in unrel_samples else 'no'

        rohs = roh_segments.get(sample, [])
        total_roh_bp = sum(e - s for _, s, e in rohs)
        n_roh = len(rohs)
        froh = total_roh_bp / genome_size if genome_size > 0 else 0

        n_total = n_het = n_hom = 0
        sum_del = sum_hr = sum_final = 0.0
        n_A = n_B = n_A_hom = n_A_roh = n_B_roh = 0
        n_del_roh = n_del_hom_roh = n_del_long_roh = 0

        for vk, vd in variants.items():
            gt_val = genotypes.get(vk, {}).get(sample, -1)
            if gt_val <= 0:
                continue

            pc = vd.get('priority_class', 'D')
            chrom = vd.get('chr', '')
            pos = int(vd.get('pos', 0)) if vd.get('pos', '').isdigit() else 0

            n_total += 1
            if gt_val == 1: n_het += 1
            if gt_val == 2: n_hom += 1

            try: sum_del += float(vd.get('deleterious_evidence_score', '0'))
            except: pass
            try: sum_hr += float(vd.get('hatchery_risk_score', '0'))
            except: pass
            try: sum_final += float(vd.get('final_priority_score', '0'))
            except: pass

            if pc == 'A':
                n_A += 1
                if gt_val == 2: n_A_hom += 1
                if in_roh(sample, chrom, pos): n_A_roh += 1
            if pc == 'B':
                n_B += 1
                if in_roh(sample, chrom, pos): n_B_roh += 1

            if in_roh(sample, chrom, pos):
                n_del_roh += 1
                if gt_val == 2: n_del_hom_roh += 1
            if in_long_roh(sample, chrom, pos):
                n_del_long_roh += 1

        out.write(f"{sample}\t{is_unrel}\t{n_total}\t{n_het}\t{n_hom}\t"
                  f"{sum_del:.1f}\t{sum_hr:.1f}\t{sum_final:.1f}\t"
                  f"{n_A}\t{n_B}\t{n_A_hom}\t{n_A_roh}\t{n_B_roh}\t"
                  f"{n_del_roh}\t{n_del_hom_roh}\t{n_del_long_roh}\t"
                  f"{total_roh_bp}\t{n_roh}\t{froh:.6f}\n")

print(f"  Output: {burden_path}")

# ═══ D. Per-ROH burden ═══
print("\n  Computing per-ROH burden...")
roh_burden_path = os.path.join(burden_dir, 'deleterious_burden_per_ROH.tsv')

with open(roh_burden_path, 'w') as out:
    out.write('sample_id\tchr\troh_start\troh_end\troh_length\t'
              'n_del_variants\tn_del_hom\tsum_priority_score\t'
              'n_classA\tn_classB\tgenes_hit\n')

    for sample in all_samples:
        for chrom, rstart, rend in roh_segments.get(sample, []):
            rlen = rend - rstart
            n_del = n_del_hom = 0
            sum_ps = 0.0
            n_ca = n_cb = 0
            genes = set()

            for vk, vd in variants.items():
                vc = vd.get('chr', '')
                vpos = int(vd.get('pos', 0)) if vd.get('pos', '').isdigit() else 0
                if vc != chrom or vpos < rstart or vpos > rend:
                    continue
                gt_val = genotypes.get(vk, {}).get(sample, -1)
                if gt_val <= 0:
                    continue
                n_del += 1
                if gt_val == 2: n_del_hom += 1
                try: sum_ps += float(vd.get('final_priority_score', '0'))
                except: pass
                pc = vd.get('priority_class', 'D')
                if pc == 'A': n_ca += 1
                if pc == 'B': n_cb += 1
                g = vd.get('gene_id', '')
                if g: genes.add(g)

            out.write(f"{sample}\t{chrom}\t{rstart}\t{rend}\t{rlen}\t"
                      f"{n_del}\t{n_del_hom}\t{sum_ps:.1f}\t{n_ca}\t{n_cb}\t"
                      f"{';'.join(sorted(genes)) if genes else 'none'}\n")

print(f"  Output: {roh_burden_path}")

# ═══ E. Per-gene burden (basic — no orthology/CAFE) ═══
print("\n  Computing per-gene burden...")
gene_burden_path = os.path.join(burden_dir, 'deleterious_burden_per_gene.tsv')

gene_stats = defaultdict(lambda: {'n_var': 0, 'n_A': 0, 'n_B': 0,
                                   'n_hom_total': 0, 'max_indiv_hom': 0,
                                   'n_roh_overlaps': 0})

for vk, vd in variants.items():
    gid = vd.get('gene_id', '')
    if not gid:
        continue
    pc = vd.get('priority_class', 'D')
    gs = gene_stats[gid]
    gs['n_var'] += 1
    if pc == 'A': gs['n_A'] += 1
    if pc == 'B': gs['n_B'] += 1

    gt = genotypes.get(vk, {})
    hom_samples = [s for s, v in gt.items() if v == 2]
    gs['n_hom_total'] += len(hom_samples)
    gs['max_indiv_hom'] = max(gs['max_indiv_hom'], len(hom_samples))

    chrom = vd.get('chr', '')
    pos = int(vd.get('pos', 0)) if vd.get('pos', '').isdigit() else 0
    for s in hom_samples:
        if in_roh(s, chrom, pos):
            gs['n_roh_overlaps'] += 1

with open(gene_burden_path, 'w') as out:
    out.write('gene_id\tn_variants\tn_classA\tn_classB\t'
              'n_hom_alt_total\tmax_indiv_with_hom_alt\tn_roh_overlaps\n')
    for gid in sorted(gene_stats.keys()):
        gs = gene_stats[gid]
        out.write(f"{gid}\t{gs['n_var']}\t{gs['n_A']}\t{gs['n_B']}\t"
                  f"{gs['n_hom_total']}\t{gs['max_indiv_hom']}\t{gs['n_roh_overlaps']}\n")

print(f"  Output: {gene_burden_path}")

# ═══ Summary ═══
n_a = sum(1 for v in variants.values() if v.get('priority_class') == 'A')
n_b = sum(1 for v in variants.values() if v.get('priority_class') == 'B')
print(f"\n=== STEP 15 Summary ===")
print(f"  Class A variants: {n_a}")
print(f"  Class B variants: {n_b}")
print(f"  Burden tables in: {burden_dir}/")
print(f"  Variant master:   {updated_path}")
print(f"\n  Three output files for C01f_c:")
print(f"    {updated_path}")
print(f"    {gt_matrix_path}")
print(f"    {burden_path}")
PYEOF

echo "=== STEP 15 complete ==="
