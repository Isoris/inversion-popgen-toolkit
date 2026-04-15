#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=64G
#SBATCH -t 0-04:00:00
#SBATCH -J manta_sum
#SBATCH -o logs/05_summary.%j.out
#SBATCH -e logs/05_summary.%j.err
# =============================================================================
# 05_summary_report.sh — v3: integrates with tiered filtering output
# =============================================================================
#
# Reads from:
#   05_final_catalogs/   — PASS VCFs + GT matrices (from step 03)
#   10_filtered_catalogs/ — summary tables with evidence columns (from filter script)
#
# Produces in 08_summary/:
#   - Per-sample burden tables (raw PASS + per tier)
#   - Per-type size distributions
#   - Sharing spectrum
#   - Per-chromosome counts
#   - Tier composition breakdown
#   - Combined multi-type tables
#   - Pairwise shared DEL matrix
#   - 1-Mb window density tables
#   - Text report
#
# Run AFTER: 03_merge_and_split.sh AND filter_all_sv_types.sh
# =============================================================================
set -euo pipefail
source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${SCRIPT_DIR}/../00_module4g_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
set -a
source "${CONFIG}"
set +a

mv_init_dirs
FILT_DIR="${OUTDIR}/10_filtered_catalogs"
mv_log "=== SUMMARY REPORT v3 ==="

SV_TYPES=(DEL DUP INV INS_small INS_large BND)

# =============================================================================
# SECTION 1: Master counts table (raw PASS → tier1 → tier2 → tier3)
# =============================================================================
mv_log "--- Section 1: Master counts ---"

MASTER_COUNTS="${DIR_SUMMARY}/master_counts.tsv"
{
  echo -e "cohort\tsvtype\traw_total\tPASS\ttier1\ttier2\ttier3\tprecise_pass\timprecise_pass"
  for cohort_n in 226 81; do
    for svtype in "${SV_TYPES[@]}"; do
      raw_vcf="${DIR_SPLIT}/catalog_${cohort_n}.${svtype}.vcf.gz"
      pass_vcf="${DIR_FINAL}/catalog_${cohort_n}.${svtype}.PASS.vcf.gz"
      summary="${FILT_DIR}/${cohort_n}.${svtype}.summary.tsv"

      n_raw=0; [[ -f "${raw_vcf}" ]] && n_raw=$(bcftools view -H "${raw_vcf}" 2>/dev/null | wc -l)
      n_pass=0; [[ -f "${pass_vcf}" ]] && n_pass=$(bcftools view -H "${pass_vcf}" 2>/dev/null | wc -l)

      n_t1=0; n_t2=0; n_t3=0; n_prec=0; n_imp=0
      if [[ -f "${summary}" ]]; then
        n_prec=$(tail -n +2 "${summary}" | awk -F'\t' '$9=="PRECISE"' | wc -l)
        n_imp=$(tail -n +2 "${summary}" | awk -F'\t' '$9=="IMPRECISE"' | wc -l)
      fi
      for tf in "${FILT_DIR}/${cohort_n}.${svtype}.tier1_"*.tsv; do
        [[ -f "$tf" ]] && n_t1=$(tail -n +2 "$tf" | wc -l); done
      for tf in "${FILT_DIR}/${cohort_n}.${svtype}.tier2_"*.tsv; do
        [[ -f "$tf" ]] && n_t2=$(tail -n +2 "$tf" | wc -l); done
      for tf in "${FILT_DIR}/${cohort_n}.${svtype}.tier3_"*.tsv; do
        [[ -f "$tf" ]] && n_t3=$(tail -n +2 "$tf" | wc -l); done

      echo -e "${cohort_n}\t${svtype}\t${n_raw}\t${n_pass}\t${n_t1}\t${n_t2}\t${n_t3}\t${n_prec}\t${n_imp}"
    done
  done
} > "${MASTER_COUNTS}"

mv_log "Master counts:"
column -t "${MASTER_COUNTS}" | while IFS= read -r line; do mv_log "  ${line}"; done

# =============================================================================
# SECTION 2: Per-sample SV burden at each tier
# =============================================================================
mv_log "--- Section 2: Per-sample burden per tier ---"

for cohort_n in 226 81; do
  python3 << PYEOF
import os, glob

filt_dir = "${FILT_DIR}"
summary_dir = "${DIR_SUMMARY}"
cohort = "${cohort_n}"
sv_types = ["DEL", "DUP", "INV", "INS_small", "INS_large", "BND"]

# For each tier level, build per-sample counts from the summary tables
for tier_label in ["PASS", "tier1", "tier2", "tier3"]:
    data = {}  # sample -> {type: count}
    samples_ordered = None

    for svtype in sv_types:
        if tier_label == "PASS":
            # Use the full summary table, count all
            path = os.path.join(filt_dir, f"{cohort}.{svtype}.summary.tsv")
        else:
            # Find the tier file
            pattern = os.path.join(filt_dir, f"{cohort}.{svtype}.{tier_label}_*.tsv")
            files = glob.glob(pattern)
            path = files[0] if files else None

        if not path or not os.path.exists(path):
            continue

        # Read the summary table and count carriers per sample
        # We need the GT matrix for per-sample counts — use the PASS VCF GT matrix
        # For tiered data, we count from the summary table's carrier count
        # But for per-SAMPLE breakdown we need the actual GT matrix
        pass

    # Actually, per-sample counts need the GT matrix. Let's use the PASS GT matrices
    # and also build tiered per-sample by subsetting sites.

# Simpler approach: build per-sample counts from GT matrices for PASS level,
# then from tiered VCFs for tier2/tier3.
import gzip

for tier_label in ["PASS", "tier2", "tier3"]:
    per_sample = {}

    for svtype in sv_types:
        if tier_label == "PASS":
            vcf_path = f"${DIR_FINAL}/catalog_{cohort}.{svtype}.PASS.vcf.gz"
        else:
            pattern = os.path.join(filt_dir, f"{cohort}.{svtype}.{tier_label}.vcf.gz")
            vcf_files = glob.glob(pattern)
            vcf_path = vcf_files[0] if vcf_files else None

        if not vcf_path or not os.path.exists(vcf_path):
            continue

        try:
            with gzip.open(vcf_path, 'rt') as f:
                sample_names = []
                for line in f:
                    if line.startswith('##'):
                        continue
                    if line.startswith('#CHROM'):
                        sample_names = line.strip().split('\t')[9:]
                        for s in sample_names:
                            if s not in per_sample:
                                per_sample[s] = {t: 0 for t in sv_types}
                        continue
                    fields = line.strip().split('\t')
                    gts = fields[9:]
                    for i, gt_field in enumerate(gts):
                        gt = gt_field.split(':')[0]
                        if gt not in ('./.', '0/0', './0', '0/.', '.', '0|0'):
                            if i < len(sample_names):
                                per_sample[sample_names[i]][svtype] += 1
        except Exception as e:
            print(f"  Warning: {vcf_path}: {e}")

    if not per_sample:
        continue

    out = os.path.join(summary_dir, f"per_sample_{tier_label}_{cohort}.tsv")
    with open(out, 'w') as o:
        o.write("sample\t" + "\t".join(sv_types) + "\ttotal\n")
        for sample in sorted(per_sample.keys()):
            counts = [per_sample[sample].get(t, 0) for t in sv_types]
            total = sum(counts)
            o.write(sample + "\t" + "\t".join(str(c) for c in counts) + f"\t{total}\n")
    print(f"  {out}: {len(per_sample)} samples × {tier_label}")

PYEOF
done

# =============================================================================
# SECTION 3: Size distributions + sharing + per-chr (from filter summary tables)
# =============================================================================
mv_log "--- Section 3: Distributions from filter summary tables ---"

for cohort_n in 226 81; do
  for svtype in DEL DUP INV INS_small INS_large; do
    summary="${FILT_DIR}/${cohort_n}.${svtype}.summary.tsv"
    [[ -f "${summary}" ]] || continue

    # SVLEN distribution
    tail -n +2 "${summary}" | awk -F'\t' '$12 > 0 {print $12}' \
      > "${DIR_SUMMARY}/${svtype}_svlen_${cohort_n}.tsv"

    # Per-chromosome
    tail -n +2 "${summary}" | awk -F'\t' '{print $1}' | sort | uniq -c \
      | awk '{print $2"\t"$1}' > "${DIR_SUMMARY}/${svtype}_per_chr_${cohort_n}.tsv"

    # Sharing (carrier count distribution)
    tail -n +2 "${summary}" | awk -F'\t' '{print $14}' | sort -n | uniq -c \
      | awk '{print $2"\t"$1}' > "${DIR_SUMMARY}/${svtype}_carrier_dist_${cohort_n}.tsv"

    # Private vs shared
    python3 << PYEOF
summary = "${summary}"
out = "${DIR_SUMMARY}/${svtype}_sharing_${cohort_n}.tsv"
sharing = []
with open(summary) as f:
    next(f)
    for line in f:
        p = line.strip().split('\t')
        try:
            sharing.append(int(p[13]))  # carriers column
        except:
            pass
ns = int(p[12]) if sharing else 0  # n_samples
priv = sum(1 for s in sharing if s == 1)
s2_5 = sum(1 for s in sharing if 2 <= s <= 5)
s6_20 = sum(1 for s in sharing if 6 <= s <= 20)
s21p = sum(1 for s in sharing if s > 20)
fixed = sum(1 for s in sharing if s == ns)
n = len(sharing)
with open(out, 'w') as o:
    o.write("category\tcount\tpercent\n")
    for cat, cnt in [("private", priv), ("shared_2_5", s2_5), ("shared_6_20", s6_20),
                     ("shared_21plus", s21p), ("fixed", fixed)]:
        pct = 100.0 * cnt / n if n > 0 else 0
        o.write(f"{cat}\t{cnt}\t{pct:.1f}\n")
PYEOF

    # Evidence distribution (sum_total_alt)
    tail -n +2 "${summary}" | awk -F'\t' '{print $21}' | sort -n | uniq -c \
      | awk '{print $2"\t"$1}' > "${DIR_SUMMARY}/${svtype}_evidence_dist_${cohort_n}.tsv"

    # Precision breakdown
    tail -n +2 "${summary}" | awk -F'\t' '{print $9}' | sort | uniq -c \
      | awk '{print $2"\t"$1}' > "${DIR_SUMMARY}/${svtype}_precision_${cohort_n}.tsv"
  done
done
mv_log "  Distributions written for all types × both cohorts"

# =============================================================================
# SECTION 4: Tier composition tables (for stacked bar plots)
# =============================================================================
mv_log "--- Section 4: Tier composition ---"

for cohort_n in 226 81; do
  TIER_COMP="${DIR_SUMMARY}/tier_composition_${cohort_n}.tsv"
  {
    echo -e "svtype\tPASS_total\ttier1\ttier2\ttier3\texcluded_by_tier1\ttier1_only\ttier2_only\tprecise_in_pass\timprecise_in_pass"
    for svtype in "${SV_TYPES[@]}"; do
      summary="${FILT_DIR}/${cohort_n}.${svtype}.summary.tsv"
      [[ -f "${summary}" ]] || { echo -e "${svtype}\t0\t0\t0\t0\t0\t0\t0\t0\t0"; continue; }

      n_pass=$(tail -n +2 "${summary}" | wc -l)
      n_prec=$(tail -n +2 "${summary}" | awk -F'\t' '$9=="PRECISE"' | wc -l)
      n_imp=$((n_pass - n_prec))

      n_t1=0; n_t2=0; n_t3=0
      for tf in "${FILT_DIR}/${cohort_n}.${svtype}.tier1_"*.tsv; do
        [[ -f "$tf" ]] && n_t1=$(tail -n +2 "$tf" | wc -l); done
      for tf in "${FILT_DIR}/${cohort_n}.${svtype}.tier2_"*.tsv; do
        [[ -f "$tf" ]] && n_t2=$(tail -n +2 "$tf" | wc -l); done
      for tf in "${FILT_DIR}/${cohort_n}.${svtype}.tier3_"*.tsv; do
        [[ -f "$tf" ]] && n_t3=$(tail -n +2 "$tf" | wc -l); done

      excluded=$((n_pass - n_t1))
      t1_only=$((n_t1 - n_t2))
      t2_only=$((n_t2 - n_t3))

      echo -e "${svtype}\t${n_pass}\t${n_t1}\t${n_t2}\t${n_t3}\t${excluded}\t${t1_only}\t${t2_only}\t${n_prec}\t${n_imp}"
    done
  } > "${TIER_COMP}"
  mv_log "  ${cohort_n}: ${TIER_COMP}"
  column -t "${TIER_COMP}" | while IFS= read -r line; do mv_log "    ${line}"; done
done

# =============================================================================
# SECTION 5: 1-Mb window density (for genome heatmap)
# =============================================================================
mv_log "--- Section 5: Window density ---"

for cohort_n in 226 81; do
  for svtype in DEL DUP INV INS_small INS_large; do
    summary="${FILT_DIR}/${cohort_n}.${svtype}.summary.tsv"
    [[ -f "${summary}" ]] || continue

    tail -n +2 "${summary}" | awk -F'\t' 'BEGIN{OFS="\t"}{
      w = int($2 / 1000000) * 1000000
      print $1, w
    }' | sort | uniq -c | awk 'BEGIN{OFS="\t"}{print $2,$3,$3+1000000,$1}' \
    > "${DIR_SUMMARY}/${svtype}_windows_1Mb_${cohort_n}.tsv"
  done
done
mv_log "  Window tables written"

# =============================================================================
# SECTION 6: Pairwise shared DEL matrix (226 cohort, for heatmap)
# =============================================================================
mv_log "--- Section 6: Pairwise DEL sharing (226) ---"

DEL_MATRIX="${DIR_FINAL}/catalog_226.DEL.PASS.GT_matrix.tsv"
PAIRWISE_OUT="${DIR_SUMMARY}/pairwise_shared_DEL_226.tsv"

if [[ -f "${DEL_MATRIX}" ]]; then
  python3 << PYEOF
import numpy as np

matrix = "${DEL_MATRIX}"
out = "${PAIRWISE_OUT}"

with open(matrix) as f:
    header = f.readline().strip().split('\t')
    samples = header[6:]
    ns = len(samples)
    gt_rows = []
    for line in f:
        p = line.strip().split('\t')
        gts = p[6:]
        binary = []
        for gt in gts:
            binary.append(1 if gt not in ('./.', '0/0', './0', '0/.', '.', '0|0') else 0)
        gt_rows.append(binary)

gt_arr = np.array(gt_rows, dtype=np.uint8)
print(f"  Matrix: {gt_arr.shape[0]} DELs × {gt_arr.shape[1]} samples")
shared = gt_arr.T @ gt_arr

with open(out, 'w') as o:
    o.write("sample_i\tsample_j\tn_shared\n")
    for i in range(ns):
        for j in range(i, ns):
            o.write(f"{samples[i]}\t{samples[j]}\t{shared[i,j]}\n")

print(f"  Pairwise: {out}")
PYEOF
else
  mv_log "  DEL GT matrix not found — skipping pairwise"
fi

# =============================================================================
# SECTION 7: INS_large assembly summary
# =============================================================================
mv_log "--- Section 7: INS_large assembly ---"

for cohort_n in 226 81; do
  INS_VCF="${DIR_FINAL}/catalog_${cohort_n}.INS_large.PASS.vcf.gz"
  [[ -f "${INS_VCF}" ]] || continue

  python3 << PYEOF
import gzip
vcf = "${INS_VCF}"
out = "${DIR_SUMMARY}/INS_large_assembly_${cohort_n}.tsv"
with gzip.open(vcf, 'rt') as v, open(out, 'w') as o:
    o.write("chrom\tpos\tid\thas_left\thas_right\tleft_len\tright_len\ttotal_assembled_bp\n")
    for l in v:
        if l.startswith('#'): continue
        p = l.strip().split('\t')
        info = dict(f.split('=', 1) for f in p[7].split(';') if '=' in f)
        ls = info.get('LEFT_SVINSSEQ', '')
        rs = info.get('RIGHT_SVINSSEQ', '')
        ll, rl = (len(ls) if ls else 0), (len(rs) if rs else 0)
        o.write(f"{p[0]}\t{p[1]}\t{p[2]}\t{'Y' if ll>0 else 'N'}\t{'Y' if rl>0 else 'N'}\t{ll}\t{rl}\t{ll+rl}\n")
PYEOF
  mv_log "  INS_large assembly ${cohort_n}: done"
done

# =============================================================================
# SECTION 8: Text report
# =============================================================================
mv_log "--- Section 8: Text report ---"

REPORT="${DIR_SUMMARY}/manta_SV_summary_report.txt"
cat > "${REPORT}" << 'REOF_START'
================================================================================
MANTA SV CATALOG — SUMMARY REPORT (v3)
================================================================================
REOF_START

cat >> "${REPORT}" << REOF
Generated: $(date '+%F %T')
Reference: fClaHyb_Gar_LG.fa | Cohort: 226 (81 unrelated) | ~9X lcWGS
Manta config: minCandidateVariantSize = 50
Filtering: evidence-based (PR + SR alt counts), NOT precision-gated

MASTER COUNTS:
$(column -t "${MASTER_COUNTS}")

TIER COMPOSITION (226):
$(column -t "${DIR_SUMMARY}/tier_composition_226.tsv")

TIER COMPOSITION (81):
$(column -t "${DIR_SUMMARY}/tier_composition_81.tsv")

TIER DEFINITIONS:
  tier1: sum_alt≥4 + type-specific size + ≥1 carrier
  tier2: sum_alt≥6 + type-specific size + ≥2 carriers + QUAL≥50
  tier3: sum_alt≥10 + min_carrier_alt≥3 + size + ≥3 carriers + QUAL≥100

  IMPRECISE = breakpoint coords approximate (no split-read assembly).
  IMPRECISE calls are KEPT — they represent real SVs at ~9X coverage.

NOTE ON SV TYPES:
  DEL       = Deletions (≥50 bp)
  DUP       = Tandem duplications
  INV       = Inversions (BND→INV via convertInversion_py3.py)
  BND       = Inter-chromosomal translocations
  INS_small = Fully assembled insertions (SVLEN present)
  INS_large = Incompletely assembled insertions (LEFT/RIGHT_SVINSSEQ)
================================================================================
REOF

mv_log "=== SUMMARY COMPLETE ==="
cat "${REPORT}"
