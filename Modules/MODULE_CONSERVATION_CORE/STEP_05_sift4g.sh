#!/usr/bin/env bash
#SBATCH --job-name=CONS_05_sift4g
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=slurm_logs/05_sift4g_%j.out
#SBATCH --error=slurm_logs/05_sift4g_%j.err
set -euo pipefail

# =============================================================================
# STEP 05 — Build SIFT4G Custom Database + Annotate
# =============================================================================
# A. Download UniRef90 (if not present)
# B. Build custom SIFT4G database from catfish genome + annotation
# C. Annotate missense variants with SIFT4G scores
#
# IMPORTANT: SIFT4G is ONE catfish-centered functional layer.
# It is built once from the reference genome/annotation, NOT per-sample.
# It only applies to protein-changing coding variants.
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config/pipeline.config.sh"

DB_BUILD="${DIR_SIFT4G}/db"
DB_OUT="${DIR_SIFT4G}/db/fClaHyb_Gar_LG"
UNIREF_DIR="${DIR_SIFT4G}/uniref90"

mkdir -p "${DB_BUILD}" "${UNIREF_DIR}" "${DIR_SIFT4G}/annotated" "${DIR_LOGS}"

# ─────────────────────────────────────────
# A. Download UniRef90 (alignment database for SIFT)
# ─────────────────────────────────────────
echo "=== A. Checking UniRef90 ==="

UNIREF_FA="${UNIREF_DIR}/uniref90.fasta"
if [[ -f "${UNIREF_FA}" ]]; then
    echo "  [SKIP] UniRef90 already downloaded"
else
    echo "  [DOWNLOAD] UniRef90 (~30 GB, this will take hours)..."
    cd "${UNIREF_DIR}"
    wget -c "https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz" \
        -O uniref90.fasta.gz
    echo "  Decompressing..."
    gunzip uniref90.fasta.gz
    echo "  UniRef90 ready: $(grep -c '^>' "${UNIREF_FA}") sequences"
fi

# Build BLAST database for UniRef90
UNIREF_DB="${UNIREF_DIR}/uniref90"
if [[ ! -f "${UNIREF_DB}.phr" && ! -f "${UNIREF_DB}.psq" ]]; then
    echo "  Building BLAST database for UniRef90..."
    makeblastdb -in "${UNIREF_FA}" -dbtype prot -out "${UNIREF_DB}" -parse_seqids
fi

# ─────────────────────────────────────────
# B. Build SIFT4G Custom Database
# ─────────────────────────────────────────
echo ""
echo "=== B. Building SIFT4G database for fClaHyb_Gar_LG ==="

SIFT4G_CREATE="${SIFT4G_DIR}/SIFT4G_Create_Genomic_DB"

if [[ -d "${DB_OUT}" && -f "${DB_OUT}/chr_list.txt" ]]; then
    echo "  [SKIP] SIFT4G database already built"
else
    # Create configuration file for SIFT4G database builder
    SIFT_CONFIG="${DB_BUILD}/sift4g_config.txt"
    cat > "${SIFT_CONFIG}" <<SIFTCFG
GENETIC_CODE_TABLE=1
GENETIC_CODE_TABLENAME=Standard
MITO_GENETIC_CODE_TABLE=2
MITO_GENETIC_CODE_TABLENAME=Vertebrate Mitochondrial

# Paths
PARENT_DIR=${DB_BUILD}
ORG=fClaHyb_Gar_LG
ORG_VERSION=1.0
DBSNP_VCF_FILE=

# Reference
GENE_DOWNLOAD_SITE=FILE
ENSEMBL_DOWNLOAD=NO

# Source files
FASTA_DIR=${DB_BUILD}/chr_fasta
GENE_ANNOTATION_DIR=${DB_BUILD}/gene_annotation

# SIFT4G executable
SIFT4G_PATH=${SIFT4G_BIN}

# Protein database for alignment
PROTEIN_DB=${UNIREF_DB}

# Subst file
SUBST_DIR=${DB_OUT}

# Number of threads
NUM_THREADS=16
SIFTCFG

    # Prepare chromosome FASTAs (SIFT4G needs individual chr files)
    CHR_DIR="${DB_BUILD}/chr_fasta"
    mkdir -p "${CHR_DIR}"
    echo "  Extracting per-chromosome FASTAs..."
    for CHR in "${CHROMS[@]}"; do
        if [[ ! -f "${CHR_DIR}/${CHR}.fa" ]]; then
            samtools faidx "${REF_FASTA}" "${CHR}" > "${CHR_DIR}/${CHR}.fa"
        fi
    done

    # Prepare gene annotation (GTF format required by SIFT4G builder)
    GENE_DIR="${DB_BUILD}/gene_annotation"
    mkdir -p "${GENE_DIR}"
    if [[ ! -f "${GENE_DIR}/fClaHyb_Gar_LG.gtf" ]]; then
        cp "${REF_GTF}" "${GENE_DIR}/fClaHyb_Gar_LG.gtf"
    fi

    # Run SIFT4G database builder
    echo "  Running SIFT4G_Create_Genomic_DB..."
    echo "  This may take 12-24 hours depending on genome size and UniRef90..."

    cd "${SIFT4G_CREATE}"
    perl scripts/SIFT4G_Create_Genomic_DB.pl \
        -config "${SIFT_CONFIG}" \
        2>&1 | tee "${DIR_LOGS}/05_sift4g_build.log"

    cd "${MODCONS}"
fi

# ─────────────────────────────────────────
# C. Annotate Variants with SIFT4G
# ─────────────────────────────────────────
echo ""
echo "=== C. Annotating variants with SIFT4G ==="

SIFT_ANNOTATOR="${SIFT4G_DIR}/SIFT4G_Annotator.jar"

for CHR in "${CHROMS[@]}"; do
    INPUT="${DIR_VARIANTS}/normalized/${CHR}.clair3.norm.vcf.gz"
    OUTPUT="${DIR_SIFT4G}/annotated/${CHR}.sift4g.vcf"
    TSV="${DIR_SIFT4G}/annotated/${CHR}.sift4g.tsv"

    if [[ ! -f "${INPUT}" ]]; then
        continue
    fi
    if [[ -f "${TSV}" ]]; then
        echo "  [SKIP] ${CHR} — already annotated"
        continue
    fi

    echo "  [SIFT4G] ${CHR}..."

    # Decompress VCF for SIFT4G annotator (it prefers plain VCF)
    PLAIN_VCF="${DIR_SIFT4G}/annotated/${CHR}.input.vcf"
    bcftools view "${INPUT}" > "${PLAIN_VCF}"

    # Run SIFT4G Annotator
    java -Xmx8g -jar "${SIFT_ANNOTATOR}" \
        -c \
        -i "${PLAIN_VCF}" \
        -d "${DB_OUT}" \
        -r "${DIR_SIFT4G}/annotated/${CHR}_results" \
        -t \
        2>"${DIR_LOGS}/05_sift4g_annotate_${CHR}.err" || {
            echo "    WARNING: SIFT4G annotation failed for ${CHR}"
            continue
        }

    # Parse SIFT4G output to TSV
    # SIFT4G adds SIFT_SCORE and SIFT_PRED to VCF INFO
    rm -f "${PLAIN_VCF}"
done

# ─────────────────────────────────────────
# D. Consolidate SIFT4G results → single TSV
# ─────────────────────────────────────────
echo ""
echo "=== D. Consolidating SIFT4G results ==="

python3 - "${DIR_SIFT4G}/annotated" "${DIR_SIFT4G}" <<'PYEOF'
"""Consolidate SIFT4G annotation results into a single TSV."""
import sys, os, glob, csv

annotated_dir = sys.argv[1]
outdir = sys.argv[2]

out_path = os.path.join(outdir, 'sift4g_scores.tsv')
header = ['chr', 'pos', 'ref', 'alt', 'var_key', 'gene_id', 'transcript_id',
          'aa_change', 'sift_score', 'sift_class', 'sift_median_info']

n = 0
with open(out_path, 'w') as out:
    out.write('\t'.join(header) + '\n')

    # Look for SIFT4G result files (typically *_SIFTpredictions.tsv or annotated VCFs)
    for result_dir in sorted(glob.glob(os.path.join(annotated_dir, '*_results'))):
        for sift_file in glob.glob(os.path.join(result_dir, '*SIFT*predictions*')):
            print(f"  Parsing {sift_file}...")
            with open(sift_file) as fh:
                for line in fh:
                    if line.startswith('#') or line.startswith('Coordinates'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) < 7:
                        continue
                    # SIFT4G output format varies; parse flexibly
                    try:
                        coord = parts[0]  # chr:pos
                        ref_alt = parts[1] if len(parts) > 1 else ''
                        gene = parts[2] if len(parts) > 2 else ''
                        tx = parts[3] if len(parts) > 3 else ''
                        aa = parts[4] if len(parts) > 4 else ''
                        score_str = parts[5] if len(parts) > 5 else ''
                        median = parts[6] if len(parts) > 6 else ''

                        # Parse score
                        try:
                            score = float(score_str)
                            sclass = 'DELETERIOUS' if score < 0.05 else 'TOLERATED'
                        except (ValueError, TypeError):
                            score = ''
                            sclass = ''

                        # Parse coordinates
                        if ':' in coord:
                            chrom, pos = coord.split(':')[:2]
                        else:
                            chrom, pos = coord, ''

                        ref, alt = '', ''
                        if '/' in ref_alt:
                            ref, alt = ref_alt.split('/')[:2]

                        var_key = f"{chrom}:{pos}:{ref}:{alt}"
                        row = [chrom, pos, ref, alt, var_key, gene, tx, aa,
                               str(score), sclass, median]
                        out.write('\t'.join(row) + '\n')
                        n += 1
                    except Exception:
                        continue

print(f"\n  Total SIFT4G annotations: {n}")
print(f"  Output: {out_path}")
PYEOF

echo "=== STEP 05 complete ==="
echo "  SIFT4G database: ${DB_OUT}"
echo "  SIFT4G scores:   ${DIR_SIFT4G}/sift4g_scores.tsv"
