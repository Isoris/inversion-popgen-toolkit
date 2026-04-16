#!/usr/bin/env bash
# =============================================================================
# prepare_gaps_agp.sh — Generate gaps BED + AGP for PHASE_01C block_detect
#
# Reads the reference FASTA, detects N-runs (gaps between contigs), and
# produces two files:
#
#   1. <ref>.gaps.bed  — 3-column BED (chrom, start, end) of all N-runs
#                        0-based half-open coordinates, no header
#
#   2. <ref>.agp       — Standard AGP with W (contig) and N (gap) rows
#                        01C reads W rows for contig boundary positions
#
# Usage:
#   bash prepare_gaps_agp.sh /path/to/fClaHyb_Gar_LG.fa
#   bash prepare_gaps_agp.sh /path/to/fClaHyb_Gar_LG.fa [min_gap_size]
#
# The min_gap_size (default: 1) filters gaps shorter than N bp.
# For T2T assemblies with no gaps, both files will be empty/trivial
# and 01C gracefully skips assembly-error annotation.
#
# Dependencies: awk, standard POSIX tools (no python/seqtk needed)
# =============================================================================

set -euo pipefail

REF="${1:?Usage: prepare_gaps_agp.sh <reference.fa> [min_gap_size]}"
MIN_GAP="${2:-1}"

if [[ ! -f "${REF}" ]]; then
  echo "[ERROR] Reference not found: ${REF}" >&2
  exit 1
fi

BASENAME="${REF%.fa}"
BASENAME="${BASENAME%.fasta}"
GAPS_BED="${BASENAME}.gaps.bed"
AGP_FILE="${BASENAME}.agp"

echo "[prepare] Reference: ${REF}"
echo "[prepare] Min gap size: ${MIN_GAP}"
echo "[prepare] Output gaps: ${GAPS_BED}"
echo "[prepare] Output AGP:  ${AGP_FILE}"

# =============================================================================
# STEP 1: Detect N-runs → gaps BED
# =============================================================================
# Linearise FASTA, scan for N/n runs, output BED (0-based, half-open).
# Works on multi-line FASTA of any line width.

echo "[prepare] Detecting N-gaps..."

awk -v min_gap="${MIN_GAP}" '
BEGIN { OFS="\t" }
/^>/ {
  # Flush previous chromosome
  if (chr != "" && in_gap && (pos - gap_start) >= min_gap)
    print chr, gap_start, pos
  chr = substr($1, 2)
  pos = 0
  in_gap = 0
  next
}
{
  for (i = 1; i <= length($0); i++) {
    c = substr($0, i, 1)
    if (c == "N" || c == "n") {
      if (!in_gap) { gap_start = pos; in_gap = 1 }
    } else {
      if (in_gap) {
        gap_len = pos - gap_start
        if (gap_len >= min_gap) print chr, gap_start, pos
        in_gap = 0
      }
    }
    pos++
  }
}
END {
  if (chr != "" && in_gap && (pos - gap_start) >= min_gap)
    print chr, gap_start, pos
}
' "${REF}" > "${GAPS_BED}"

N_GAPS=$(wc -l < "${GAPS_BED}")
echo "[prepare] Found ${N_GAPS} gap(s)"

if [[ "${N_GAPS}" -eq 0 ]]; then
  echo "[prepare] No gaps found — assembly appears T2T"
  echo "[prepare] ${GAPS_BED} is empty (01C will skip gap annotation)"
fi

# =============================================================================
# STEP 2: Build AGP from FASTA + gaps
# =============================================================================
# Standard AGP format (https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/)
#
# For each chromosome:
#   - Split into contigs at N-gaps
#   - W rows = contigs, N rows = gaps
#   - Coordinates are 1-based inclusive
#
# If no gaps: single W row spanning the whole chromosome.

echo "[prepare] Building AGP..."

# First, get chromosome lengths from .fai or compute on the fly
FAI="${REF}.fai"
if [[ ! -f "${FAI}" ]]; then
  echo "[prepare] No .fai found, computing lengths from FASTA..."
  # Compute lengths inline
  awk '/^>/ { if (chr) print chr"\t"len; chr=substr($1,2); len=0; next }
       { len += length($0) }
       END { if (chr) print chr"\t"len }' "${REF}" > "${BASENAME}.chrlen.tmp"
  LENFILE="${BASENAME}.chrlen.tmp"
else
  # Use fai (cols: name, length, ...)
  awk -F'\t' '{ print $1"\t"$2 }' "${FAI}" > "${BASENAME}.chrlen.tmp"
  LENFILE="${BASENAME}.chrlen.tmp"
fi

# Now build AGP: merge chrlen + gaps
awk -F'\t' -v OFS='\t' -v gaps_file="${GAPS_BED}" '
BEGIN {
  # Load gaps into arrays keyed by chrom
  while ((getline line < gaps_file) > 0) {
    split(line, f, "\t")
    chr = f[1]; gs = f[2]+0; ge = f[3]+0  # 0-based
    n_gaps[chr]++
    gap_starts[chr, n_gaps[chr]] = gs
    gap_ends[chr, n_gaps[chr]] = ge
  }
  close(gaps_file)
}
{
  chr = $1; chrlen = $2+0
  part = 0
  prev_end = 0  # 0-based end of last gap (= start of next contig)
  ng = n_gaps[chr]+0
  contig_num = 0

  if (ng == 0) {
    # No gaps: single contig = whole chromosome
    part++
    contig_num++
    print chr, 1, chrlen, part, "W", chr "_ctg" contig_num, 1, chrlen, "+"
  } else {
    for (g = 1; g <= ng; g++) {
      gs = gap_starts[chr, g]  # 0-based start
      ge = gap_ends[chr, g]    # 0-based end

      # Contig before this gap
      if (gs > prev_end) {
        part++
        contig_num++
        ctg_start_1 = prev_end + 1       # 1-based
        ctg_end_1 = gs                     # 1-based (gs is 0-based = last base +1, so this is correct as end)
        ctg_len = ctg_end_1 - ctg_start_1 + 1
        print chr, ctg_start_1, ctg_end_1, part, "W", chr "_ctg" contig_num, 1, ctg_len, "+"
      }

      # Gap
      part++
      gap_start_1 = gs + 1  # 1-based
      gap_end_1 = ge         # 1-based
      gap_len = gap_end_1 - gap_start_1 + 1
      print chr, gap_start_1, gap_end_1, part, "N", gap_len, "scaffold", "yes", "na"

      prev_end = ge  # 0-based end of gap
    }

    # Contig after last gap
    if (prev_end < chrlen) {
      part++
      contig_num++
      ctg_start_1 = prev_end + 1
      ctg_end_1 = chrlen
      ctg_len = ctg_end_1 - ctg_start_1 + 1
      print chr, ctg_start_1, ctg_end_1, part, "W", chr "_ctg" contig_num, 1, ctg_len, "+"
    }
  }
}
' "${LENFILE}" > "${AGP_FILE}"

# Cleanup
rm -f "${BASENAME}.chrlen.tmp"

N_CONTIGS=$(awk -F'\t' '$5=="W"' "${AGP_FILE}" | wc -l)
N_AGP_GAPS=$(awk -F'\t' '$5=="N"' "${AGP_FILE}" | wc -l)
N_CHRS=$(awk -F'\t' '$5=="W"' "${AGP_FILE}" | cut -f1 | sort -u | wc -l)

echo "[prepare] AGP: ${N_CONTIGS} contigs, ${N_AGP_GAPS} gaps across ${N_CHRS} chromosomes"

# =============================================================================
# STEP 3: Validate
# =============================================================================

echo ""
echo "[prepare] === SUMMARY ==="
echo "  Gaps BED:  ${GAPS_BED} (${N_GAPS} gaps)"
echo "  AGP:       ${AGP_FILE} (${N_CONTIGS} contigs, ${N_AGP_GAPS} gaps)"
echo ""

if [[ "${N_GAPS}" -gt 0 ]]; then
  echo "[prepare] Top gaps by size:"
  awk -F'\t' '{ print $1, $3-$2, "bp" }' "${GAPS_BED}" | sort -k2,2nr | head -10
  echo ""
fi

echo "[prepare] Usage with 01C:"
echo "  Rscript PHASE_01C_block_detect.R <precomp_dir> <outdir> \\"
echo "    --gaps ${GAPS_BED} \\"
echo "    --agp ${AGP_FILE} \\"
echo "    --mode hatchery"
echo ""
echo "[prepare] Or via launcher (auto-detected if files are at expected paths):"
echo "  sbatch LAUNCH_01C_block_detect.slurm"
echo ""
echo "[DONE]"
