#!/usr/bin/env bash
#SBATCH --job-name=cargo_60_genebed
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/60_genebed_%j.out
#SBATCH --error=logs/60_genebed_%j.err
set -euo pipefail

# =============================================================================
# STEP_60 — Build gene_intervals.bed.gz from REF_GFF3
# =============================================================================
# One-time setup. Idempotent — skips if outputs exist.
# Produces:
#   ${GENE_BED}                  — sorted, bgzipped, tabix-indexed gene intervals
#   ${GENE_LENGTH_TSV}           — gene_id  length_bp  cds_length_bp  n_transcripts
#   ${TRANSCRIPT_GENE_MAP}       — transcript_id  gene_id (one row per transcript)
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_cargo_config.sh"

mkdir -p logs

# Skip if all outputs exist
if [[ -f "${GENE_BED}" && -f "${GENE_BED_TBI}" && \
      -f "${GENE_LENGTH_TSV}" && -f "${TRANSCRIPT_GENE_MAP}" ]]; then
  echo "[STEP_60] All outputs exist — skipping"
  echo "  ${GENE_BED}"
  echo "  ${GENE_LENGTH_TSV}"
  echo "  ${TRANSCRIPT_GENE_MAP}"
  exit 0
fi

echo "[STEP_60] Parsing GFF3: ${REF_GFF3}"

python3 - "${REF_GFF3}" "${GENE_BED%.gz}" "${GENE_LENGTH_TSV}" "${TRANSCRIPT_GENE_MAP}" <<'PYEOF'
"""Parse GFF3 → gene BED + length table + transcript→gene map.

Robust to the two common GFF3 attribute conventions:
  - Ensembl-style:  ID=gene:GENEID; biotype=protein_coding
  - NCBI-style:     ID=gene-GENEID; gene_biotype=protein_coding

Only protein-coding genes are emitted to the BED. Non-coding genes (lncRNA,
miRNA, etc.) are excluded because all downstream cargo analyses operate on
missense / LoF variants.
"""
import sys, gzip, re
from collections import defaultdict

gff_path, bed_path, len_path, txmap_path = sys.argv[1:5]


def parse_attrs(s):
    out = {}
    for kv in s.strip(';').split(';'):
        if '=' in kv:
            k, v = kv.split('=', 1)
            out[k.strip()] = v.strip()
    return out


def strip_prefix(gid):
    # ID=gene:ENSXX → ENSXX, ID=gene-XYZ → XYZ
    if ':' in gid:
        return gid.split(':', 1)[1]
    if gid.startswith('gene-') or gid.startswith('rna-'):
        return gid.split('-', 1)[1]
    return gid


genes = {}                # gene_id → (chrom, start_0, end, strand)
gene_biotype = {}         # gene_id → biotype
transcripts = {}          # tx_id → gene_id
cds_lengths = defaultdict(int)  # tx_id → sum of CDS lengths

opener = gzip.open if gff_path.endswith('.gz') else open
with opener(gff_path, 'rt') as fh:
    for line in fh:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 9:
            continue
        chrom, src, ftype, start, end, score, strand, phase, attrs_s = parts
        try:
            start_i = int(start) - 1   # GFF is 1-based inclusive → BED 0-based half-open
            end_i = int(end)
        except ValueError:
            continue
        attrs = parse_attrs(attrs_s)

        if ftype == 'gene':
            gid_raw = attrs.get('ID', '')
            if not gid_raw:
                continue
            gid = strip_prefix(gid_raw)
            biot = (attrs.get('biotype') or attrs.get('gene_biotype')
                    or attrs.get('type') or '')
            genes[gid] = (chrom, start_i, end_i, strand)
            gene_biotype[gid] = biot

        elif ftype in ('mRNA', 'transcript'):
            tx_raw = attrs.get('ID', '')
            par_raw = attrs.get('Parent', '')
            if not tx_raw or not par_raw:
                continue
            tx_id = strip_prefix(tx_raw)
            # Parent can be comma-separated; take the first
            par = par_raw.split(',')[0]
            gene_id = strip_prefix(par)
            transcripts[tx_id] = gene_id

        elif ftype == 'CDS':
            par_raw = attrs.get('Parent', '')
            if not par_raw:
                continue
            for par in par_raw.split(','):
                tx_id = strip_prefix(par)
                cds_lengths[tx_id] += (end_i - start_i)

# ── Filter genes to protein-coding ──
def is_protein_coding(gid):
    biot = gene_biotype.get(gid, '').lower()
    if biot == 'protein_coding':
        return True
    # Fallback: gene has at least one transcript with CDS
    for tx, g in transcripts.items():
        if g == gid and cds_lengths.get(tx, 0) > 0:
            return True
    return False

pc_genes = {g: v for g, v in genes.items() if is_protein_coding(g)}
print(f"  Total genes parsed: {len(genes)}", file=sys.stderr)
print(f"  Protein-coding kept: {len(pc_genes)}", file=sys.stderr)

# ── Write BED (sorted by chrom, start) ──
bed_rows = sorted(
    [(g, *v) for g, v in pc_genes.items()],
    key=lambda r: (r[1], r[2])
)
with open(bed_path, 'w') as out:
    for gid, chrom, start, end, strand in bed_rows:
        out.write(f"{chrom}\t{start}\t{end}\t{gid}\t0\t{strand}\n")
print(f"  BED rows: {len(bed_rows)} → {bed_path}", file=sys.stderr)

# ── Write length table ──
gene_n_tx = defaultdict(int)
gene_max_cds = defaultdict(int)
for tx, gid in transcripts.items():
    if gid in pc_genes:
        gene_n_tx[gid] += 1
        if cds_lengths.get(tx, 0) > gene_max_cds[gid]:
            gene_max_cds[gid] = cds_lengths[tx]

with open(len_path, 'w') as out:
    out.write("gene_id\tlength_bp\tcds_length_bp\tn_transcripts\n")
    for gid, chrom, start, end, strand in bed_rows:
        out.write(f"{gid}\t{end-start}\t{gene_max_cds.get(gid, 0)}\t{gene_n_tx.get(gid, 0)}\n")
print(f"  Length table → {len_path}", file=sys.stderr)

# ── Write transcript → gene map ──
with open(txmap_path, 'w') as out:
    out.write("transcript_id\tgene_id\tcds_length_bp\n")
    for tx, gid in sorted(transcripts.items()):
        if gid in pc_genes:
            out.write(f"{tx}\t{gid}\t{cds_lengths.get(tx, 0)}\n")
print(f"  Transcript map → {txmap_path}", file=sys.stderr)
PYEOF

# Compress + tabix-index
echo "[STEP_60] Compressing and indexing BED..."
bgzip -f "${GENE_BED%.gz}"
tabix -p bed -f "${GENE_BED}"

echo "[STEP_60] Done"
echo "  ${GENE_BED}             $(zcat "${GENE_BED}" | wc -l) genes"
echo "  ${GENE_LENGTH_TSV}      $(($(wc -l < "${GENE_LENGTH_TSV}") - 1)) rows"
echo "  ${TRANSCRIPT_GENE_MAP}  $(($(wc -l < "${TRANSCRIPT_GENE_MAP}") - 1)) transcripts"
