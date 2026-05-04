#!/bin/bash
# =============================================================================
# scripts/build_ref_n_bed.sh — extract N / assembly-gap regions from a FASTA
# =============================================================================
# Emits a 3-column BED of every run of N characters in the reference, per
# chromosome. Runs shorter than MIN_GAP_BP (default 100) are skipped.
#
# Usage:
#   bash scripts/build_ref_n_bed.sh <reference.fasta> [min_gap_bp]
# Outputs:
#   ref_n.ALL.bed                  genome-wide
#   ref_n.<CHR>.bed                one per chromosome
#
# Uses awk (not Biostrings / Python) so it works on any HPC without R setup.
# =============================================================================
set -euo pipefail

FASTA="${1:-}"
MIN_GAP_BP="${2:-100}"
[[ -z "${FASTA}" ]] && { echo "Usage: $0 <reference.fasta> [min_gap_bp]" >&2; exit 1; }
[[ -f "${FASTA}" ]] || { echo "Not found: ${FASTA}" >&2; exit 1; }

OUT_DIR="$(dirname "${FASTA}")/ref_n_bed"
mkdir -p "${OUT_DIR}"
ALL_BED="${OUT_DIR}/ref_n.ALL.bed"

echo "[build_ref_n_bed] scanning ${FASTA} for N-runs >= ${MIN_GAP_BP} bp"

awk -v MINLEN="${MIN_GAP_BP}" -v OUTDIR="${OUT_DIR}" '
  BEGIN { chr = ""; pos = 0; in_gap = 0; gap_start = 0; total = 0 }
  /^>/ {
    # Flush any open gap on the previous chrom
    if (in_gap && (pos - gap_start) >= MINLEN) {
      print chr "\t" gap_start "\t" pos >> (OUTDIR "/ref_n.ALL.bed")
      print chr "\t" gap_start "\t" pos >> (OUTDIR "/ref_n." chr ".bed")
      total++
    }
    chr = substr($1, 2)
    pos = 0; in_gap = 0; gap_start = 0
    # Truncate per-chrom BED at start so re-runs are idempotent
    printf "" > (OUTDIR "/ref_n." chr ".bed")
    next
  }
  {
    line_len = length($0)
    for (i = 1; i <= line_len; i++) {
      c = substr($0, i, 1)
      if (c == "N" || c == "n") {
        if (!in_gap) { in_gap = 1; gap_start = pos + i - 1 }
      } else {
        if (in_gap) {
          gap_end = pos + i - 1
          if ((gap_end - gap_start) >= MINLEN) {
            print chr "\t" gap_start "\t" gap_end >> (OUTDIR "/ref_n.ALL.bed")
            print chr "\t" gap_start "\t" gap_end >> (OUTDIR "/ref_n." chr ".bed")
            total++
          }
          in_gap = 0
        }
      }
    }
    pos += line_len
  }
  END {
    if (in_gap && (pos - gap_start) >= MINLEN) {
      print chr "\t" gap_start "\t" pos >> (OUTDIR "/ref_n.ALL.bed")
      print chr "\t" gap_start "\t" pos >> (OUTDIR "/ref_n." chr ".bed")
      total++
    }
    print "[build_ref_n_bed] found " total " N-runs >= " MINLEN " bp" > "/dev/stderr"
  }
' < <(
  if [[ "${FASTA}" == *.gz ]]; then zcat "${FASTA}"; else cat "${FASTA}"; fi
)

# Also emit an empty-file-guard: if no per-chrom BED was created for a chrom,
# create one as an empty BED (so Q04 can cheaply check existence).
# We do that by listing chroms from a .fai if present.
if [[ -f "${FASTA}.fai" ]]; then
  awk '{ print $1 }' "${FASTA}.fai" | while read -r c; do
    [[ -f "${OUT_DIR}/ref_n.${c}.bed" ]] || : > "${OUT_DIR}/ref_n.${c}.bed"
  done
fi

echo "[build_ref_n_bed] BED files in ${OUT_DIR}/"
wc -l "${ALL_BED}" 2>/dev/null || true
