#!/usr/bin/env bash
#SBATCH --job-name=CONS_04_snpeff
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=slurm_logs/04_snpeff_%j.out
#SBATCH --error=slurm_logs/04_snpeff_%j.err
set -euo pipefail

# =============================================================================
# STEP 04 — Build SnpEff Custom DB + Annotate Variants
# =============================================================================
# A. Build custom SnpEff database for fClaHyb_Gar_LG
# B. Annotate Clair3 VCFs (per-chrom, SLURM array-compatible)
# C. Annotate SV VCFs
# D. Parse ANN fields → flat TSV for merging
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config/pipeline.config.sh"

GENOME_NAME="fClaHyb_Gar_LG"
SNPEFF_DATA="${SNPEFF_DIR}/data"
DB_DIR="${SNPEFF_DATA}/${GENOME_NAME}"

mkdir -p "${DIR_SNPEFF}" "${DB_DIR}" "${DIR_LOGS}"

# ─────────────────────────────────────────
# A. BUILD CUSTOM SNPEFF DATABASE
# ─────────────────────────────────────────
echo "=== A. Building SnpEff database: ${GENOME_NAME} ==="

if [[ -f "${DB_DIR}/snpEffectPredictor.bin" ]]; then
    echo "  [SKIP] Database already built"
else
    # 1. Add genome entry to snpEff.config
    CONFIG="${SNPEFF_DIR}/snpEff.config"
    if ! grep -q "^${GENOME_NAME}.genome" "${CONFIG}" 2>/dev/null; then
        echo "" >> "${CONFIG}"
        echo "# Catfish hybrid (Clarias gariepinus × C. macrocephalus)" >> "${CONFIG}"
        echo "${GENOME_NAME}.genome : Clarias hybrid (Gar haplotype)" >> "${CONFIG}"
        echo "${GENOME_NAME}.chromosomes : $(IFS=','; echo "${CHROMS[*]}")" >> "${CONFIG}"
        echo "  Added ${GENOME_NAME} to snpEff.config"
    fi

    # 2. Copy files into SnpEff data directory
    # SnpEff expects: data/<genome>/genes.gff (GFF3), sequences.fa (reference)
    cp "${REF_GFF3}" "${DB_DIR}/genes.gff"
    cp "${REF_FASTA}" "${DB_DIR}/sequences.fa"
    [[ -f "${REF_CDS}" ]] && cp "${REF_CDS}" "${DB_DIR}/cds.fa"
    [[ -f "${REF_PEP}" ]] && cp "${REF_PEP}" "${DB_DIR}/protein.fa"

    # 3. Build database
    echo "  Building database..."
    java -Xmx24g -jar "${SNPEFF_JAR}" build \
        -gff3 \
        -v "${GENOME_NAME}" \
        2>&1 | tee "${DIR_LOGS}/04_snpeff_build.log"

    if [[ -f "${DB_DIR}/snpEffectPredictor.bin" ]]; then
        echo "  Database built successfully"
    else
        echo "  ERROR: Database build failed"
        exit 1
    fi
fi

# ─────────────────────────────────────────
# B. ANNOTATE CLAIR3 VARIANTS (per-chrom)
# ─────────────────────────────────────────
echo ""
echo "=== B. Annotating Clair3 variants with SnpEff ==="

for CHR in "${CHROMS[@]}"; do
    INPUT="${DIR_VARIANTS}/normalized/${CHR}.clair3.norm.vcf.gz"
    OUTPUT="${DIR_SNPEFF}/${CHR}.snpeff.vcf.gz"
    STATS="${DIR_SNPEFF}/${CHR}.snpeff.stats.html"

    if [[ ! -f "${INPUT}" ]]; then
        echo "  [SKIP] ${CHR} — no normalized VCF"
        continue
    fi
    if [[ -f "${OUTPUT}" ]]; then
        echo "  [SKIP] ${CHR} — already annotated"
        continue
    fi

    echo "  [ANNOTATE] ${CHR}..."
    java -Xmx8g -jar "${SNPEFF_JAR}" ann \
        -v \
        -stats "${STATS}" \
        -csvStats "${DIR_SNPEFF}/${CHR}.snpeff.stats.csv" \
        -nodownload \
        -canon \
        "${GENOME_NAME}" \
        "${INPUT}" \
    | bcftools view -Oz -o "${OUTPUT}"
    bcftools index -t "${OUTPUT}"

    N=$(bcftools view -H "${OUTPUT}" | wc -l)
    echo "    ${CHR}: ${N} annotated variants"
done

# ─────────────────────────────────────────
# C. ANNOTATE SV VCFs
# ─────────────────────────────────────────
echo ""
echo "=== C. Annotating SV variants with SnpEff ==="

for SV_VCF in "${DIR_VARIANTS}/sv_combined/all_sv.sorted.vcf.gz" \
              "${DIR_VARIANTS}/sv_combined/delly_all_sv.sorted.vcf.gz" \
              "${DIR_VARIANTS}/manta_merged/manta_all_sv.sorted.vcf.gz"; do
    if [[ ! -f "${SV_VCF}" ]]; then
        continue
    fi

    BASENAME=$(basename "${SV_VCF}" .vcf.gz)
    OUTPUT="${DIR_SNPEFF}/${BASENAME}.snpeff.vcf.gz"

    if [[ -f "${OUTPUT}" ]]; then
        echo "  [SKIP] ${BASENAME} — already annotated"
        continue
    fi

    echo "  [ANNOTATE] ${BASENAME}..."
    java -Xmx8g -jar "${SNPEFF_JAR}" ann \
        -v \
        -nodownload \
        -canon \
        "${GENOME_NAME}" \
        "${SV_VCF}" \
    | bcftools view -Oz -o "${OUTPUT}"
    bcftools index -t "${OUTPUT}"
done

# ─────────────────────────────────────────
# D. PARSE ANN FIELDS → FLAT TSV
# ─────────────────────────────────────────
echo ""
echo "=== D. Parsing SnpEff ANN → TSV ==="

python3 - "${DIR_SNPEFF}" "${CANONICAL_TX}" <<'PYEOF'
"""Parse SnpEff ANN fields from annotated VCFs into flat TSV."""
import sys, gzip, os, glob

snpeff_dir = sys.argv[1]
canon_file = sys.argv[2]

# Load canonical transcripts
canonical = set()
if os.path.exists(canon_file):
    with open(canon_file) as f:
        canonical = {line.strip() for line in f if line.strip()}
    print(f"  Loaded {len(canonical)} canonical transcripts")

# ANN field format (SnpEff):
# Allele|Annotation|Annotation_Impact|Gene_Name|Gene_ID|Feature_Type|Feature_ID|
# Transcript_BioType|Rank|HGVS.c|HGVS.p|cDNA.pos/cDNA.length|CDS.pos/CDS.length|
# AA.pos/AA.length|Distance|ERRORS
ANN_FIELDS = ['allele', 'annotation', 'impact', 'gene_name', 'gene_id',
              'feature_type', 'feature_id', 'transcript_biotype', 'rank',
              'hgvs_c', 'hgvs_p', 'cdna_pos', 'cds_pos', 'aa_pos', 'distance', 'errors']

header = ['chr', 'pos', 'ref', 'alt', 'var_key'] + ANN_FIELDS + ['is_canonical']

out_all = os.path.join(snpeff_dir, 'snpeff_all_consequences.tsv')
out_canon = os.path.join(snpeff_dir, 'snpeff_canonical_consequences.tsv')

n_total = 0
n_canon = 0

with open(out_all, 'w') as fa, open(out_canon, 'w') as fc:
    fa.write('\t'.join(header) + '\n')
    fc.write('\t'.join(header) + '\n')

    # Process all annotated VCFs
    vcf_files = sorted(glob.glob(os.path.join(snpeff_dir, '*.snpeff.vcf.gz')))
    for vcf_path in vcf_files:
        basename = os.path.basename(vcf_path)
        print(f"  Parsing {basename}...")

        with gzip.open(vcf_path, 'rt') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 8:
                    continue
                chrom, pos, _, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
                info = parts[7]
                var_key = f"{chrom}:{pos}:{ref}:{alt}"

                # Extract ANN field
                ann_str = ''
                for field in info.split(';'):
                    if field.startswith('ANN='):
                        ann_str = field[4:]
                        break

                if not ann_str:
                    continue

                # Parse each annotation
                for ann in ann_str.split(','):
                    fields = ann.split('|')
                    # Pad to expected length
                    while len(fields) < len(ANN_FIELDS):
                        fields.append('')

                    is_canon = 'yes' if fields[6] in canonical else 'no'
                    row = [chrom, pos, ref, alt, var_key] + fields[:len(ANN_FIELDS)] + [is_canon]
                    fa.write('\t'.join(row) + '\n')
                    n_total += 1

                    if is_canon == 'yes' or not canonical:
                        fc.write('\t'.join(row) + '\n')
                        n_canon += 1

print(f"\n  Total annotations: {n_total}")
print(f"  Canonical annotations: {n_canon}")
print(f"  All → {out_all}")
print(f"  Canonical → {out_canon}")
PYEOF

echo "=== STEP 04 complete ==="
