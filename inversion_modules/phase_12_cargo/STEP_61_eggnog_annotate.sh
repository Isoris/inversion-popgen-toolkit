#!/usr/bin/env bash
#SBATCH --job-name=cargo_61_eggnog
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/61_eggnog_%j.out
#SBATCH --error=logs/61_eggnog_%j.err
set -euo pipefail

# =============================================================================
# STEP_61 — eggNOG-mapper functional annotation (one-time, OPTIONAL)
# =============================================================================
# Required for STEP_C63 functional enrichment (6B). If you skip this step,
# enrichment will be limited to gene-family / paralog grouping only and GO/KEGG
# columns will be empty.
#
# Outputs:
#   ${EGGNOG_OUT}            — raw eggNOG-mapper output (.emapper.annotations)
#   ${GENE_FUNCTION_TSV}     — per-gene resolved table:
#                                gene_id  preferred_name  best_OG  description
#                                go_terms (semicolon-sep)  kegg_pathway  family
#
# Idempotent — skips if outputs exist.
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_cargo_config.sh"

mkdir -p logs

if [[ -f "${GENE_FUNCTION_TSV}" ]]; then
  echo "[STEP_61] ${GENE_FUNCTION_TSV} already exists — skipping"
  exit 0
fi

: "${REF_PEP:?REF_PEP must be set in 00_inversion_config.sh}"
: "${TRANSCRIPT_GENE_MAP:?Run STEP_60 first}"

# eggNOG-mapper paths — adjust if your install differs
EGGNOG_DATA_DIR="${EGGNOG_DATA_DIR:-${BASE}/13-programs/eggnog_data}"
EMAPPER_BIN="${EMAPPER_BIN:-emapper.py}"

if ! command -v "${EMAPPER_BIN}" >/dev/null 2>&1; then
  echo "[STEP_61] WARNING: ${EMAPPER_BIN} not found in PATH"
  echo "  Install via:  conda install -c bioconda eggnog-mapper"
  echo "  Then download DB:  download_eggnog_data.py --data_dir ${EGGNOG_DATA_DIR}"
  echo "[STEP_61] Skipping eggNOG run — STEP_C63 will run with family-only enrichment"
  # Write a stub gene_function table so downstream doesn't crash
  awk 'BEGIN{OFS="\t"; print "gene_id","preferred_name","best_OG","description","go_terms","kegg_pathway","family"}
       NR>1 {seen[$2]=1} END{for (g in seen) print g,"","","","","",""}' \
       "${TRANSCRIPT_GENE_MAP}" > "${GENE_FUNCTION_TSV}"
  exit 0
fi

WORKDIR="${CARGO_DIR}/eggnog_work"
mkdir -p "${WORKDIR}"

echo "[STEP_61] Running eggNOG-mapper on ${REF_PEP}"
echo "  Threads: ${SLURM_CPUS_PER_TASK:-32}"
echo "  Data dir: ${EGGNOG_DATA_DIR}"

"${EMAPPER_BIN}" \
  -i "${REF_PEP}" \
  --output cargo_eggnog \
  --output_dir "${WORKDIR}" \
  --data_dir "${EGGNOG_DATA_DIR}" \
  --cpu "${SLURM_CPUS_PER_TASK:-32}" \
  -m diamond \
  --override

cp "${WORKDIR}/cargo_eggnog.emapper.annotations" "${EGGNOG_OUT}"
echo "[STEP_61] eggNOG raw output → ${EGGNOG_OUT}"

# ── Resolve protein → gene → consensus annotation ──
echo "[STEP_61] Resolving per-gene annotations..."
python3 - "${EGGNOG_OUT}" "${TRANSCRIPT_GENE_MAP}" "${GENE_FUNCTION_TSV}" <<'PYEOF'
"""
Aggregate eggNOG protein-level annotations to gene level.

eggNOG outputs one row per query protein. Multiple transcripts of the same
gene may produce conflicting annotations; we resolve by:
  - preferred_name: take the most common non-empty value across transcripts
  - best_OG: take the most common non-empty value
  - go_terms: union across transcripts (semicolon-separated, sorted)
  - kegg_pathway: union across transcripts
  - family: parsed from PFAMs / best_OG description
"""
import sys
from collections import defaultdict, Counter

eggnog_path, txmap_path, out_path = sys.argv[1:4]

# ── Load transcript → gene ──
tx_gene = {}
with open(txmap_path) as f:
    next(f)
    for line in f:
        parts = line.rstrip('\n').split('\t')
        if len(parts) >= 2:
            tx_gene[parts[0]] = parts[1]

# ── Parse eggNOG annotations ──
# Columns (eggNOG-mapper v2.1+):
# query, seed_ortholog, evalue, score, eggNOG_OGs, max_annot_lvl, COG_category,
# Description, Preferred_name, GOs, EC, KEGG_ko, KEGG_Pathway, KEGG_Module,
# KEGG_Reaction, KEGG_rclass, BRITE, KEGG_TC, CAZy, BiGG_Reaction, PFAMs
gene_data = defaultdict(lambda: {
    'preferred_name': Counter(), 'best_OG': Counter(), 'description': Counter(),
    'go_terms': set(), 'kegg_pathway': set(), 'pfams': set()
})

header_idx = {}
with open(eggnog_path) as f:
    for line in f:
        line = line.rstrip('\n')
        if line.startswith('##'):
            continue
        if line.startswith('#query'):
            cols = line.lstrip('#').split('\t')
            header_idx = {c: i for i, c in enumerate(cols)}
            continue
        if not line.strip():
            continue
        parts = line.split('\t')
        if not header_idx:
            continue
        def get(col):
            i = header_idx.get(col, -1)
            return parts[i] if 0 <= i < len(parts) else ''
        query = get('query')
        gene_id = tx_gene.get(query)
        if not gene_id:
            # try stripping version suffix
            gene_id = tx_gene.get(query.split('.')[0])
        if not gene_id:
            continue
        gd = gene_data[gene_id]
        for col_src, key in [('Preferred_name', 'preferred_name'),
                             ('eggNOG_OGs', 'best_OG'),
                             ('Description', 'description')]:
            v = get(col_src).strip()
            if v and v != '-':
                # eggNOG_OGs is comma-separated; take the first (most specific)
                if key == 'best_OG':
                    v = v.split(',')[0]
                gd[key][v] += 1
        for col_src, key in [('GOs', 'go_terms'),
                             ('KEGG_Pathway', 'kegg_pathway'),
                             ('PFAMs', 'pfams')]:
            v = get(col_src).strip()
            if v and v != '-':
                for tok in v.split(','):
                    tok = tok.strip()
                    if tok:
                        gd[key].add(tok)

# ── Write resolved per-gene table ──
def top(counter):
    return counter.most_common(1)[0][0] if counter else ''

with open(out_path, 'w') as out:
    out.write("gene_id\tpreferred_name\tbest_OG\tdescription\tgo_terms\tkegg_pathway\tfamily\n")
    # Include all genes from the transcript map, even those eggNOG didn't annotate
    all_genes = set(tx_gene.values())
    for gid in sorted(all_genes):
        gd = gene_data.get(gid)
        if gd is None:
            out.write(f"{gid}\t\t\t\t\t\t\n")
            continue
        # Family heuristic: most common PFAM, or the OG name fragment
        family = ''
        if gd['pfams']:
            family = sorted(gd['pfams'])[0]
        elif gd['best_OG']:
            family = top(gd['best_OG']).split('@')[0]
        out.write('\t'.join([
            gid,
            top(gd['preferred_name']),
            top(gd['best_OG']),
            top(gd['description']),
            ';'.join(sorted(gd['go_terms'])),
            ';'.join(sorted(gd['kegg_pathway'])),
            family,
        ]) + '\n')

print(f"[STEP_61] Resolved {len(gene_data)} annotated genes / {len(set(tx_gene.values()))} total → {out_path}", file=sys.stderr)
PYEOF

echo "[STEP_61] Done"
echo "  ${GENE_FUNCTION_TSV}    $(($(wc -l < "${GENE_FUNCTION_TSV}") - 1)) gene rows"
