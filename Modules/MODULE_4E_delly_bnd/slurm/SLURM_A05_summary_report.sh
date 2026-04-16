#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=32G
#SBATCH -t 0-02:00:00
#SBATCH -J bnd_sum
#SBATCH -o logs/05_summary.%j.out
#SBATCH -e logs/05_summary.%j.err
set -euo pipefail
source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${SCRIPT_DIR}/../00_module4e_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
set -a
source "${CONFIG}"
set +a

dv_init_dirs
dv_log "=== FINAL SUMMARY (BND) ==="

REPORT="${DIR_SUMMARY}/delly_BND_summary_report.txt"
COUNTS="${DIR_SUMMARY}/delly_BND_counts.tsv"

# Counts
FINAL_226_VCF="${DIR_FINAL}/catalog_226.BND.vcf.gz"
N_226=0
[[ -f "${FINAL_226_VCF}" ]] && N_226=$(bcftools view -H "${FINAL_226_VCF}" | wc -l)

FINAL_81_PASS=$(ls "${DIR_FINAL}"/catalog_81.BND.*.PASS.vcf.gz 2>/dev/null | head -1 || true)
N_81_PASS=0
[[ -n "${FINAL_81_PASS}" && -f "${FINAL_81_PASS}" ]] && N_81_PASS=$(bcftools view -H "${FINAL_81_PASS}" | wc -l)

FUNC="${DIR_ANNOT}/catalog_226.functional_class.tsv"
N_CDS=0; N_EXON=0; N_INTRONIC=0; N_INTERGENIC=0; N_GENES=0
if [[ -f "${FUNC}" ]]; then
  N_CDS=$(awk -F'\t' 'NR>1&&$10=="CDS_overlap"' "${FUNC}" | wc -l)
  N_EXON=$(awk -F'\t' 'NR>1&&$10=="exon_overlap"' "${FUNC}" | wc -l)
  N_INTRONIC=$(awk -F'\t' 'NR>1&&$10=="intronic"' "${FUNC}" | wc -l)
  N_INTERGENIC=$(awk -F'\t' 'NR>1&&$10=="intergenic"' "${FUNC}" | wc -l)
  N_GENES=$(awk -F'\t' 'NR>1&&$10!="intergenic"{print $11}' "${FUNC}" | tr ',' '\n' | grep -v '^\.$' | sort -u | wc -l)
fi

N_REP=0; N_NOREP=0
[[ -f "${DIR_ANNOT}/catalog_226.BNDs_in_repeats.bed" ]] && N_REP=$(wc -l < "${DIR_ANNOT}/catalog_226.BNDs_in_repeats.bed")
[[ -f "${DIR_ANNOT}/catalog_226.BNDs_not_in_repeats.bed" ]] && N_NOREP=$(wc -l < "${DIR_ANNOT}/catalog_226.BNDs_not_in_repeats.bed")

# Per-sample + sharing
MATRIX_226="${DIR_FINAL}/catalog_226.BND.GT_matrix.tsv"
PER_SAMPLE="${DIR_SUMMARY}/per_sample_BND_counts.tsv"
PRIV_SHARED="${DIR_SUMMARY}/private_vs_shared_BND.tsv"
PER_CHR="${DIR_SUMMARY}/per_chromosome_BND_counts.tsv"
PAIRWISE_SHARED="${DIR_SUMMARY}/pairwise_shared_BND.tsv"
PCA_INPUT="${DIR_SUMMARY}/BND_binary_genotype_matrix.tsv"
WINDOW_COUNTS="${DIR_SUMMARY}/BND_window_counts_1Mb.tsv"

if [[ -f "${MATRIX_226}" ]]; then
  python3 << PYEOF
import numpy as np
from collections import defaultdict

matrix = "${MATRIX_226}"
per_out = "${PER_SAMPLE}"
ps_out = "${PRIV_SHARED}"
chr_out = "${PER_CHR}"
pairwise_out = "${PAIRWISE_SHARED}"
pca_out = "${PCA_INPUT}"
window_out = "${WINDOW_COUNTS}"

with open(matrix) as f:
    header = f.readline().rstrip("\n").split("\t")
    samples = header[5:]
    ns = len(samples)

    sample_counts = [0] * ns
    sharing = []
    chr_counts = defaultdict(int)
    gt_matrix = []
    window_data = defaultdict(int)

    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 6:
            continue

        chrom = p[0]
        pos = int(p[1])
        gts = p[5:]

        binary = []
        n_alt = 0
        for i, gt in enumerate(gts):
            has_alt = gt not in ('./.','0/0','./0','0/.','.','0|0')
            if has_alt:
                sample_counts[i] += 1
                n_alt += 1
            binary.append(1 if has_alt else 0)

        sharing.append(n_alt)
        chr_counts[chrom] += 1
        gt_matrix.append(binary)

        w_start = (pos // 1000000) * 1000000
        window_data[(chrom, w_start)] += 1

with open(per_out, "w") as o:
    o.write("sample\tn_BNDs\n")
    for s, c in zip(samples, sample_counts):
        o.write(f"{s}\t{c}\n")

priv = sum(1 for s in sharing if s == 1)
s2_5 = sum(1 for s in sharing if 2 <= s <= 5)
s6_20 = sum(1 for s in sharing if 6 <= s <= 20)
s21p = sum(1 for s in sharing if s > 20)
fixed = sum(1 for s in sharing if s == ns)

with open(ps_out, "w") as o:
    o.write("category\tcount\n")
    o.write(f"private\t{priv}\n")
    o.write(f"shared_2_5\t{s2_5}\n")
    o.write(f"shared_6_20\t{s6_20}\n")
    o.write(f"shared_21plus\t{s21p}\n")
    o.write(f"fixed\t{fixed}\n")

with open(chr_out, "w") as o:
    o.write("chrom\tn_BNDs\n")
    for c in sorted(chr_counts, key=lambda x: (x.replace('C_gar_LG',''), x)):
        o.write(f"{c}\t{chr_counts[c]}\n")

with open(window_out, "w") as o:
    o.write("chrom\twindow_start\twindow_end\tn_SVs\n")
    for (c, ws), cnt in sorted(window_data.items()):
        o.write(f"{c}\t{ws}\t{ws+1000000}\t{cnt}\n")

print("Computing pairwise shared BNDs...")
if len(gt_matrix) > 0:
    gt_arr = np.array(gt_matrix, dtype=np.uint8)
    shared = gt_arr.T @ gt_arr
else:
    shared = np.zeros((ns, ns), dtype=np.uint32)

with open(pairwise_out, "w") as o:
    o.write("sample_i\tsample_j\tn_shared\n")
    for i in range(ns):
        for j in range(i, ns):
            o.write(f"{samples[i]}\t{samples[j]}\t{int(shared[i,j])}\n")

with open(pca_out, "w") as o:
    o.write("BND_ID\t" + "\t".join(samples) + "\n")
    with open(matrix) as fm:
        fm.readline()
        for line in fm:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            did = f"{p[0]}:{p[1]}"
            gts = p[5:]
            binary = ["1" if gt not in ('./.','0/0','./0','0/.','.','0|0') else "0" for gt in gts]
            o.write(did + "\t" + "\t".join(binary) + "\n")

print("Summary tables complete.")
PYEOF
fi

cat > "${REPORT}" << REOF
================================================================================
DELLY BND CATALOG — SUMMARY REPORT
Generated: $(date '+%F %T')
Reference: fClaHyb_Gar_LG.fa | Cohort: 226 (81 unrelated) | ~9X lcWGS
Type: Breakend variants (inter/intra-chromosomal rearrangements)
================================================================================
226 catalog:     ${N_226} BNDs
81 PASS catalog: ${N_81_PASS} BNDs

Functional: CDS=${N_CDS} exon=${N_EXON} intronic=${N_INTRONIC} intergenic=${N_INTERGENIC}
Unique genes: ${N_GENES}
Repeat >=50%: ${N_REP} | No repeat: ${N_NOREP}

Caveats: Catalog only. ~9X limits het sensitivity.
================================================================================
REOF

cat > "${COUNTS}" << CEOF
metric	value
total_BND_226	${N_226}
PASS_BND_81	${N_81_PASS}
CDS_overlap	${N_CDS}
exon_overlap	${N_EXON}
intronic	${N_INTRONIC}
intergenic	${N_INTERGENIC}
unique_genes	${N_GENES}
repeat_50pct	${N_REP}
no_repeat	${N_NOREP}
CEOF

dv_log "=== SUMMARY COMPLETE ==="
cat "${REPORT}"
