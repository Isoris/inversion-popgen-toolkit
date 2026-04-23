#!/usr/bin/env bash
#SBATCH --job-name=CONS_00_setup
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=slurm_logs/00_setup_%j.out
#SBATCH --error=slurm_logs/00_setup_%j.err
set -euo pipefail

# =============================================================================
# STEP 00 — Project Setup & Catfish Reference Preparation
# =============================================================================
# Creates directory tree, converts GTF→GFF3, extracts CDS/protein/cDNA,
# builds gene→transcript→protein ID mapping, selects canonical transcripts.
#
# INTEGRATED from Module 1 scripts:
#   - Picard NormalizeFasta (line length standardization, header cleanup)
#   - seqkit fx2tab for scaffold/contig sizes
#   - Gap detection logic from detgaps
#   - MD5 checksums for reproducibility
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config/pipeline.config.sh"

mkdir -p "${MODCONS}/slurm_logs"

# --- Create full directory tree ---
echo "=== Creating directory structure ==="
for d in 00_config 01_reference 02_annotation 03_variants 04_snpeff 05_sift4g \
         06_orthofinder_core 06A_wgd_ploidy_check 06B_orthofinder_fish \
         07_comparative_genomes 08_cafe 09_duplication_modes \
         10_cactus 11_maf 12_gerp 13_phast \
         16_merged_variant_tables 17_burden_tables 18_roh_overlap \
         19_stats_inputs 20_results 21_logs; do
    mkdir -p "${MODCONS}/${d}"
done
mkdir -p "${DIR_VARIANTS}"/{clair3_merged,delly_merged,manta_merged,sv_combined,normalized}
mkdir -p "${DIR_SNPEFF}/db" "${DIR_SIFT4G}"/{db,annotated}
mkdir -p "${DIR_GERP}"/{core,fish,siluriform}
mkdir -p "${DIR_PHAST}"/{core,fish,siluriform}
echo "  Done"

# --- Validate inputs ---
echo "=== Validating input files ==="
for f in "${REF_FASTA}" "${REF_FAI}" "${REF_GTF}"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Missing: $f"; exit 1
    fi
done
echo "  Reference: ${REF_FASTA}"
echo "  GTF:       ${REF_GTF}"

# ─────────────────────────────────────────
# A. Normalize reference FASTA (from Module 1 logic)
# ─────────────────────────────────────────
echo ""
echo "=== A. Normalizing reference FASTA ==="
NORM_FA="${DIR_REF}/fClaHyb_Gar_LG.normalized.fa"

if [[ -f "${NORM_FA}" ]]; then
    echo "  [SKIP] Already normalized"
else
    if command -v picard &>/dev/null; then
        echo "  Running Picard NormalizeFasta..."
        picard NormalizeFasta \
            -I "${REF_FASTA}" \
            -O "${NORM_FA}" \
            -LINE_LENGTH 80 \
            -TRUNCATE_SEQUENCE_NAMES_AT_WHITESPACE true \
            -QUIET true \
            -VERBOSITY ERROR
        # Clean headers: remove everything after first space, strip pipe prefixes
        # (from module_1_prepare_the_genomes.sh line 302)
        sed -i -E 's/^>[^|]*\|/>/; s/^(>[^ ]+).*/\1/' "${NORM_FA}"
        samtools faidx "${NORM_FA}"
        echo "  Normalized: ${NORM_FA}"
    else
        echo "  Picard not found — using samtools faidx on original"
        ln -sf "${REF_FASTA}" "${NORM_FA}"
    fi
fi

# ─────────────────────────────────────────
# B. Scaffold/contig sizes (seqkit fx2tab logic from Module 1 step 6)
# ─────────────────────────────────────────
echo ""
echo "=== B. Computing scaffold sizes ==="
SIZES_FILE="${DIR_REF}/fClaHyb_Gar_LG.scaffolds.sizes"
if [[ ! -f "${SIZES_FILE}" ]]; then
    if command -v seqkit &>/dev/null; then
        seqkit fx2tab -nl "${REF_FASTA}" > "${SIZES_FILE}"
        echo "  Sizes: $(wc -l < "${SIZES_FILE}") sequences"
    else
        # Fallback: use FAI
        awk '{print $1"\t"$2}' "${REF_FAI}" > "${SIZES_FILE}"
        echo "  Sizes from FAI: $(wc -l < "${SIZES_FILE}") sequences"
    fi
fi

# ─────────────────────────────────────────
# C. Gap detection (from Module 1 detgaps logic)
# ─────────────────────────────────────────
echo ""
echo "=== C. Detecting gaps in reference ==="
GAPS_FILE="${DIR_REF}/fClaHyb_Gar_LG.scaffolds.gaps"
if [[ ! -f "${GAPS_FILE}" ]]; then
    if command -v detgaps &>/dev/null; then
        detgaps "${REF_FASTA}" > "${GAPS_FILE}"
        echo "  Gaps: $(wc -l < "${GAPS_FILE}") gap regions"
    else
        # Simple gap detection: find runs of N
        python3 -c "
import sys
with open('${REF_FASTA}') as f:
    name, pos = '', 0
    in_gap = False
    gap_start = 0
    for line in f:
        if line.startswith('>'):
            if in_gap: print(f'{name}\t{gap_start}\t{pos}')
            name = line[1:].split()[0]
            pos = 0; in_gap = False
        else:
            for c in line.strip():
                if c in 'Nn':
                    if not in_gap: gap_start = pos; in_gap = True
                else:
                    if in_gap: print(f'{name}\t{gap_start}\t{pos}'); in_gap = False
                pos += 1
    if in_gap: print(f'{name}\t{gap_start}\t{pos}')
" > "${GAPS_FILE}"
        echo "  Gaps (Python): $(wc -l < "${GAPS_FILE}") gap regions"
    fi
fi

# Gap count per chromosome (from detgaps_c.sh logic)
GAPCOUNTS="${DIR_REF}/fClaHyb_Gar_LG.gap_counts_per_chr.tsv"
if [[ -f "${GAPS_FILE}" && ! -f "${GAPCOUNTS}" ]]; then
    echo -e "chr\tn_gaps" > "${GAPCOUNTS}"
    awk '{chr[$1]++} END{for(c in chr) print c"\t"chr[c]}' "${GAPS_FILE}" | sort >> "${GAPCOUNTS}"
fi

# ─────────────────────────────────────────
# D. MD5 checksums (from Module 1 step 9)
# ─────────────────────────────────────────
echo ""
echo "=== D. MD5 checksums ==="
MD5_FILE="${DIR_REF}/reference_files.md5"
md5sum "${REF_FASTA}" "${REF_FAI}" "${REF_GTF}" > "${MD5_FILE}" 2>/dev/null || true
echo "  Checksums: ${MD5_FILE}"

# ─────────────────────────────────────────
# E. Convert GTF to GFF3
# ─────────────────────────────────────────
echo ""
echo "=== E. Converting GTF → GFF3 ==="
if [[ ! -f "${REF_GFF3}" ]]; then
    if command -v gffread &>/dev/null; then
        gffread "${REF_GTF}" -o "${REF_GFF3}" --force-exons -F
        echo "  Converted with gffread"
    elif command -v agat_convert_sp_gff2gff.pl &>/dev/null; then
        agat_convert_sp_gff2gff.pl --gff "${REF_GTF}" -o "${REF_GFF3}"
        echo "  Converted with AGAT"
    else
        echo "ERROR: Neither gffread nor AGAT available"
        exit 1
    fi
else
    echo "  [SKIP] GFF3 exists"
fi

n_gene=$(grep -c $'\tgene\t' "${REF_GFF3}" || true)
n_cds=$(grep -c $'\tCDS\t' "${REF_GFF3}" || true)
echo "  GFF3: ${n_gene} genes, ${n_cds} CDS features"
[[ "$n_cds" -eq 0 ]] && { echo "ERROR: No CDS features in GFF3"; exit 1; }

# ─────────────────────────────────────────
# F. Extract CDS / protein / cDNA (gffread)
# ─────────────────────────────────────────
echo ""
echo "=== F. Extracting CDS / protein / cDNA ==="
if command -v gffread &>/dev/null; then
    [[ ! -f "${REF_CDS}" ]]  && gffread "${REF_GFF3}" -g "${REF_FASTA}" -x "${REF_CDS}"
    [[ ! -f "${REF_PEP}" ]]  && gffread "${REF_GFF3}" -g "${REF_FASTA}" -y "${REF_PEP}"
    [[ ! -f "${REF_CDNA}" ]] && gffread "${REF_GFF3}" -g "${REF_FASTA}" -w "${REF_CDNA}"
    echo "  CDS:     $(grep -c '^>' "${REF_CDS}" 2>/dev/null || echo 0) sequences"
    echo "  Protein: $(grep -c '^>' "${REF_PEP}" 2>/dev/null || echo 0) sequences"
    echo "  cDNA:    $(grep -c '^>' "${REF_CDNA}" 2>/dev/null || echo 0) sequences"
else
    echo "  WARNING: gffread not found"
fi

# ─────────────────────────────────────────
# G. Gene → transcript → protein ID mapping + canonical selection
# ─────────────────────────────────────────
echo ""
echo "=== G. Building ID mapping table ==="

python3 - "${REF_GFF3}" "${GENE_MAP}" "${CANONICAL_TX}" <<'PYEOF'
import sys, re
from collections import defaultdict

gff3_path, map_out, canon_out = sys.argv[1], sys.argv[2], sys.argv[3]

genes = {}
tx_parent = {}
tx_info = {}
cds_lengths = defaultdict(int)

with open(gff3_path) as fh:
    for line in fh:
        if line.startswith('#'): continue
        parts = line.strip().split('\t')
        if len(parts) < 9: continue
        chrom, _, ftype, start, end, _, strand, _, attrs = parts
        attr_dict = {}
        for a in attrs.split(';'):
            if '=' in a:
                k, v = a.split('=', 1)
                attr_dict[k.strip()] = v.strip()

        if ftype == 'gene':
            gid = attr_dict.get('ID', '')
            genes[gid] = {'chr': chrom, 'start': int(start), 'end': int(end),
                          'strand': strand,
                          'biotype': attr_dict.get('gene_biotype', attr_dict.get('biotype', 'unknown'))}
        elif ftype in ('mRNA', 'transcript'):
            tid = attr_dict.get('ID', '')
            parent = attr_dict.get('Parent', '')
            tx_parent[tid] = parent
            tx_info[tid] = {'chr': chrom, 'start': int(start), 'end': int(end), 'cds_len': 0}
        elif ftype == 'CDS':
            parent = attr_dict.get('Parent', '')
            cds_lengths[parent] += int(end) - int(start) + 1

for tid, clen in cds_lengths.items():
    if tid in tx_info: tx_info[tid]['cds_len'] = clen

gene_canonical = {}
for tid, gid in tx_parent.items():
    clen = tx_info.get(tid, {}).get('cds_len', 0)
    if gid not in gene_canonical or clen > gene_canonical[gid][1]:
        gene_canonical[gid] = (tid, clen)

with open(map_out, 'w') as out:
    out.write("gene_id\ttranscript_id\tchr\tgene_start\tgene_end\tstrand\tbiotype\tcds_length\tis_canonical\n")
    for tid, gid in sorted(tx_parent.items()):
        gi = genes.get(gid, {})
        ti = tx_info.get(tid, {})
        is_canon = 'yes' if gene_canonical.get(gid, ('',))[0] == tid else 'no'
        out.write(f"{gid}\t{tid}\t{gi.get('chr','')}\t{gi.get('start','')}\t{gi.get('end','')}\t"
                  f"{gi.get('strand','')}\t{gi.get('biotype','')}\t{ti.get('cds_len',0)}\t{is_canon}\n")

with open(canon_out, 'w') as out:
    for gid, (tid, clen) in sorted(gene_canonical.items()):
        out.write(f"{tid}\n")

print(f"  Genes: {len(genes)}, Transcripts: {len(tx_parent)}, Canonical: {len(gene_canonical)}")
PYEOF

# ─────────────────────────────────────────
# H. Reference summary report (Module 1 style)
# ─────────────────────────────────────────
echo ""
echo "=== H. Reference genome summary ==="
REPORT="${DIR_REF}/reference_summary.txt"
{
    echo "======================================"
    echo "REFERENCE GENOME SUMMARY"
    echo "======================================"
    echo "Genome:    fClaHyb_Gar_LG (Clarias gariepinus Gar subgenome, extracted from the haplotype-resolved F1 hybrid assembly)"
    echo "FASTA:     ${REF_FASTA}"
    echo ""

    # Assembly stats
    TOTAL_BP=$(awk '{s+=$2}END{print s}' "${REF_FAI}")
    N_SEQ=$(wc -l < "${REF_FAI}")
    N_CHR=$(grep -c 'C_gar_LG' "${REF_FAI}" || true)
    LARGEST=$(sort -k2 -rn "${REF_FAI}" | head -1 | awk '{print $1"\t"$2}')

    echo "Total bp:       ${TOTAL_BP}"
    echo "Sequences:      ${N_SEQ}"
    echo "Chromosomes:    ${N_CHR}"
    echo "Largest:        ${LARGEST}"

    # Gaps
    if [[ -f "${GAPS_FILE}" ]]; then
        N_GAPS=$(wc -l < "${GAPS_FILE}")
        echo "Gaps:           ${N_GAPS} gap regions"
    fi

    # Annotation
    echo ""
    echo "GFF3:           ${REF_GFF3}"
    echo "Genes:          $(grep -c $'\tgene\t' "${REF_GFF3}" || echo 0)"
    echo "CDS features:   $(grep -c $'\tCDS\t' "${REF_GFF3}" || echo 0)"
    echo "Protein seqs:   $(grep -c '^>' "${REF_PEP}" 2>/dev/null || echo 0)"
    echo ""

    # Chromosome table
    echo "--- Chromosome sizes ---"
    awk '$1 ~ /^C_gar_LG/ {printf "%s\t%d\n", $1, $2}' "${REF_FAI}" | sort -V

} > "${REPORT}"

cat "${REPORT}"

echo ""
echo "=== STEP 00 complete ==="
