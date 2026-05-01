#!/usr/bin/env bash
# =============================================================================
# run_wfmash_gar_vs_mac.sh — pairwise alignment of the Gar and Mac haplotypes
#                            for cross-species breakpoint detection.
# -----------------------------------------------------------------------------
# Input:
#   $GAR_FA = /project/lt200308-agbsci/01-catfish_assembly/02-annot/fClaHyb_Gar_LG.fa
#   $MAC_FA = /project/lt200308-agbsci/01-catfish_assembly/02-annot/fClaHyb_Mac_LG.fa
# Output:
#   ${OUTDIR}/gar_vs_mac.paf       — primary PAF
#   ${OUTDIR}/wfmash.log           — stderr capture for diagnostics
#   ${OUTDIR}/run_metadata.json    — params + timing + checksums for provenance
#
# Convention: Gar is the QUERY, Mac is the TARGET. Downstream breakpoint code
# walks Gar coords and reports which Mac chrom/orient each syntenic block
# came from. Reverse run (Mac as query) is fine but produces a different
# breakpoint catalogue — don't mix them.
#
# wfmash defaults documented in the handoff:
#   -X     all-vs-all minus self (we want cross-species so this excludes the
#          assembly self-alignment that would otherwise dominate)
#   -p 90  ≥90% ANI minimum. Cgar↔Cmac ~10-15 Mya divergence is consistent
#          with 92-95% ANI in syntenic regions. Drop to 85 if LG ends look
#          sparse in the PAF (rerun with: -p 85).
#   -s 50000  segment length 50 kb — manuscript scale; matches the inversion
#             candidate sizes in the 226-sample cohort
#   -n 1   one-to-one mapping. The handoff offers -n 5 as fallback if some
#          inversions have multiple homologous segments. Try -n 1 first.
#   -t 16  threads. LANTA login node is fine for the test pair; full run
#          should go through SLURM.
#
# Usage:
#   bash run_wfmash_gar_vs_mac.sh                 # default params
#   P=85 S=25000 bash run_wfmash_gar_vs_mac.sh    # finer + looser
#   N=5  bash run_wfmash_gar_vs_mac.sh            # top-5 mappings
#
# Sanity-check workflow before full run:
#   1. samtools faidx fClaHyb_Gar_LG.fa LG28 > test_gar_LG28.fa
#      samtools faidx fClaHyb_Mac_LG.fa LG28 > test_mac_LG28.fa
#   2. Run with -t 4 on those two single-chrom FASTAs — confirm wfmash works
#      end-to-end and produces a sensible PAF (a few hundred rows for one
#      pair of chromosome-scale sequences with known inversion at LG28).
#   3. If the test passes, run this driver as-is on the full FASTAs.
#
# =============================================================================

set -euo pipefail

# -------- Inputs --------------------------------------------------------------
GAR_FA="${GAR_FA:-/project/lt200308-agbsci/01-catfish_assembly/02-annot/fClaHyb_Gar_LG.fa}"
MAC_FA="${MAC_FA:-/project/lt200308-agbsci/01-catfish_assembly/02-annot/fClaHyb_Mac_LG.fa}"

# -------- Params (overridable via env) ----------------------------------------
P="${P:-90}"
S="${S:-50000}"
N="${N:-1}"
T="${T:-16}"

# -------- Output --------------------------------------------------------------
OUTDIR="${OUTDIR:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/cross_species_breakpoints/01_wfmash}"
mkdir -p "$OUTDIR"

# -------- Pre-flight ----------------------------------------------------------
[[ -s "$GAR_FA" ]] || { echo "FATAL: GAR_FA not readable: $GAR_FA" >&2; exit 1; }
[[ -s "$MAC_FA" ]] || { echo "FATAL: MAC_FA not readable: $MAC_FA" >&2; exit 1; }

# wfmash needs samtools-indexed FASTAs (.fai) AND BWT-style indices for the
# k-mer mapper. The .fai is fast to make if it's missing; let's create it
# defensively so we don't fail a 30-min run on a missing aux file.
if [[ ! -s "${GAR_FA}.fai" ]]; then
  echo "[run_wfmash] indexing $GAR_FA"
  samtools faidx "$GAR_FA"
fi
if [[ ! -s "${MAC_FA}.fai" ]]; then
  echo "[run_wfmash] indexing $MAC_FA"
  samtools faidx "$MAC_FA"
fi

command -v wfmash >/dev/null || {
  echo "FATAL: wfmash not in PATH. Activate the conda env first:" >&2
  echo "  conda activate assembly  # or whichever env has wfmash" >&2
  echo "If wfmash isn't installed: mamba install -c bioconda wfmash" >&2
  exit 1
}

# -------- Run -----------------------------------------------------------------
PAF="${OUTDIR}/gar_vs_mac.paf"
LOG="${OUTDIR}/wfmash.log"
META="${OUTDIR}/run_metadata.json"

echo "[run_wfmash] start $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "[run_wfmash] params: -p $P -s $S -n $N -t $T"
echo "[run_wfmash] query  : $GAR_FA"
echo "[run_wfmash] target : $MAC_FA"
echo "[run_wfmash] output : $PAF"

# Important: wfmash treats the FIRST positional arg as the target and the
# SECOND as the query. We want Gar = query, so order is: MAC_FA (target)
# then GAR_FA (query). The -X flag excludes self-mapping, which here means
# "skip pairs where target_name == query_name" — useful because the Gar
# and Mac haplotypes share LG-prefix chrom names (LG01..LG28 in Gar, LG01..LG27
# in Mac). NOTE: if your assembly used distinct chrom prefixes (e.g.
# Cgar_LG01 vs Cmac_LG01), you can drop -X without affecting correctness.
T0=$(date +%s)
{
  wfmash \
    "$MAC_FA" \
    "$GAR_FA" \
    -X \
    -p "$P" \
    -s "$S" \
    -n "$N" \
    -t "$T" \
    > "$PAF"
} 2> "$LOG"
T1=$(date +%s)
ELAPSED=$((T1 - T0))

# -------- Provenance ----------------------------------------------------------
N_ROWS=$(wc -l < "$PAF" || echo 0)
GAR_SHA=$(sha256sum "$GAR_FA" | awk '{print $1}')
MAC_SHA=$(sha256sum "$MAC_FA" | awk '{print $1}')

cat > "$META" <<JSON
{
  "tool":            "wfmash",
  "tool_version":    "$(wfmash --version 2>&1 | head -n1)",
  "generated_at":    "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "elapsed_seconds": $ELAPSED,
  "query":  { "path": "$GAR_FA", "role": "query",  "sha256": "$GAR_SHA" },
  "target": { "path": "$MAC_FA", "role": "target", "sha256": "$MAC_SHA" },
  "params": { "p": $P, "s": $S, "n": $N, "t": $T, "X": true },
  "output": { "paf": "$PAF", "n_rows": $N_ROWS }
}
JSON

echo "[run_wfmash] done ($ELAPSED s, $N_ROWS PAF rows)"
echo "[run_wfmash] PAF  : $PAF"
echo "[run_wfmash] META : $META"
echo "[run_wfmash] LOG  : $LOG"
