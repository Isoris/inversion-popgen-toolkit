#!/usr/bin/env bash
#SBATCH --job-name=CONS_04B_csq
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=slurm_logs/04B_csq_%j.out
#SBATCH --error=slurm_logs/04B_csq_%j.err
set -euo pipefail

# =============================================================================
# STEP 04B — bcftools csq (Haplotype-Aware Consequence Calling)
# =============================================================================
# Runs bcftools csq on normalized Clair3 VCFs.
# bcftools csq is haplotype-aware: it considers phased variants on the same
# haplotype together when predicting consequences (e.g., two nearby frameshifts
# that restore the frame). Complementary to SnpEff.
#
# NOTE: bcftools csq handles strand-aware CDS mapping internally.
#       Input must be left-aligned (done in STEP03). Do NOT right-align.
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config/pipeline.config.sh"

OUTDIR="${DIR_SNPEFF}/bcftools_csq"
mkdir -p "${OUTDIR}" "${DIR_LOGS}"

# ─────────────────────────────────────────
# Validate GFF3 for bcftools csq
# ─────────────────────────────────────────
echo "=== Validating GFF3 for bcftools csq ==="

# bcftools csq requires Ensembl-style GFF3 with:
# - gene features with ID=
# - mRNA/transcript features with Parent=gene_id
# - CDS features with Parent=transcript_id
# - proper biotype attributes

if [[ ! -f "${REF_GFF3}" ]]; then
    echo "ERROR: GFF3 not found: ${REF_GFF3}"
    echo "  Run 00_setup.sh first"
    exit 1
fi

echo "  GFF3: $(grep -c $'\tCDS\t' "${REF_GFF3}") CDS features"

# ─────────────────────────────────────────
# Run bcftools csq per chromosome
# ─────────────────────────────────────────
echo ""
echo "=== Running bcftools csq per chromosome ==="

for CHR in "${CHROMS[@]}"; do
    INPUT="${DIR_VARIANTS}/normalized/${CHR}.clair3.norm.vcf.gz"
    OUTPUT="${OUTDIR}/${CHR}.csq.vcf.gz"
    TSV="${OUTDIR}/${CHR}.csq.tsv.gz"

    if [[ ! -f "${INPUT}" ]]; then
        echo "  [SKIP] ${CHR} — no normalized VCF"
        continue
    fi
    if [[ -f "${OUTPUT}" ]]; then
        echo "  [SKIP] ${CHR} — already annotated"
        continue
    fi

    echo "  [CSQ] ${CHR}..."

    # Run bcftools csq
    # -l = local haplotype-aware mode (considers nearby phased variants)
    # -p a = phase: use all phase information available
    bcftools csq \
        -f "${REF_FASTA}" \
        -g "${REF_GFF3}" \
        -p a \
        -l \
        --threads 4 \
        -Oz -o "${OUTPUT}" \
        "${INPUT}" \
        2>"${DIR_LOGS}/04B_csq_${CHR}.err"

    bcftools index -t "${OUTPUT}"

    # Extract BCSQ tags to TSV
    bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/BCSQ\n' \
        "${OUTPUT}" \
    | gzip > "${TSV}"

    N=$(zcat "${TSV}" | wc -l)
    echo "    ${CHR}: ${N} annotated variants"
done

# ─────────────────────────────────────────
# Parse BCSQ → flat TSV
# ─────────────────────────────────────────
echo ""
echo "=== Parsing BCSQ tags → TSV ==="

python3 - "${OUTDIR}" <<'PYEOF'
"""Parse bcftools csq BCSQ tags into flat TSV."""
import sys, gzip, os, glob

outdir = sys.argv[1]
out_path = os.path.join(outdir, 'bcftools_csq_consequences.tsv')

# BCSQ format: consequence|gene|transcript|biotype|strand|amino_acid_change|dna_change
BCSQ_FIELDS = ['consequence', 'gene', 'transcript', 'biotype', 'strand', 'aa_change', 'dna_change']

header = ['chr', 'pos', 'ref', 'alt', 'var_key'] + BCSQ_FIELDS

n = 0
with open(out_path, 'w') as out:
    out.write('\t'.join(header) + '\n')

    for tsv_gz in sorted(glob.glob(os.path.join(outdir, '*.csq.tsv.gz'))):
        with gzip.open(tsv_gz, 'rt') as fh:
            for line in fh:
                parts = line.strip().split('\t')
                if len(parts) < 5:
                    continue
                chrom, pos, ref, alt, bcsq = parts[0], parts[1], parts[2], parts[3], parts[4]
                var_key = f"{chrom}:{pos}:{ref}:{alt}"

                if bcsq == '.' or not bcsq:
                    continue

                # BCSQ can have multiple annotations separated by ','
                for ann in bcsq.split(','):
                    fields = ann.split('|')
                    while len(fields) < len(BCSQ_FIELDS):
                        fields.append('')
                    row = [chrom, pos, ref, alt, var_key] + fields[:len(BCSQ_FIELDS)]
                    out.write('\t'.join(row) + '\n')
                    n += 1

print(f"  Total BCSQ annotations: {n}")
print(f"  Output: {out_path}")
PYEOF

echo "=== STEP 04B complete ==="
