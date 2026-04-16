#!/usr/bin/env bash
# =============================================================================
# inspect_triangle_v2.sh — Quick inspection of triangle insulation results
#
# Usage:
#   bash inspect_triangle_v2.sh [triangle_dir]
# =============================================================================

BASE="${BASE:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"
TRI_DIR="${1:-${BASE}/inversion_localpca_v7/06_mds_candidates/snake_regions_multiscale/triangles_v2}"

echo "================================================================"
echo "  Triangle v2 Results Inspection"
echo "  Dir: ${TRI_DIR}"
echo "================================================================"

if [[ ! -d "${TRI_DIR}" ]]; then
  echo "[ERROR] Directory not found: ${TRI_DIR}"
  echo "Try: ls ${BASE}/inversion_localpca_v7/06_mds_candidates/snake_regions_multiscale/triangles*"
  exit 1
fi

echo ""
echo "── Files ──"
ls -lh "${TRI_DIR}/"

# ── Intervals summary ──
INTERVALS="${TRI_DIR}/triangle_intervals.tsv.gz"
if [[ -f "${INTERVALS}" ]]; then
  echo ""
  echo "── Intervals ──"
  n_total=$(zcat "${INTERVALS}" | tail -n +2 | wc -l)
  echo "  Total intervals: ${n_total}"

  echo ""
  echo "  By type:"
  zcat "${INTERVALS}" | tail -n +2 | awk -F'\t' '{print $NF}' | sort | uniq -c | sort -rn
  # NF = last column is usually interval_type

  echo ""
  echo "  Header:"
  zcat "${INTERVALS}" | head -1 | tr '\t' '\n' | nl

  echo ""
  echo "  Top 20 by squareness (strong_triangle intervals):"
  # Find squareness column index
  sq_col=$(zcat "${INTERVALS}" | head -1 | tr '\t' '\n' | grep -n "squareness" | head -1 | cut -d: -f1)
  type_col=$(zcat "${INTERVALS}" | head -1 | tr '\t' '\n' | grep -n "interval_type" | head -1 | cut -d: -f1)
  chr_col=$(zcat "${INTERVALS}" | head -1 | tr '\t' '\n' | grep -n "chrom" | head -1 | cut -d: -f1)
  start_col=$(zcat "${INTERVALS}" | head -1 | tr '\t' '\n' | grep -n "start_mb\|start_bp" | head -1 | cut -d: -f1)
  end_col=$(zcat "${INTERVALS}" | head -1 | tr '\t' '\n' | grep -n "end_mb\|end_bp" | head -1 | cut -d: -f1)

  if [[ -n "${sq_col}" ]]; then
    echo "  (squareness col=${sq_col}, type col=${type_col})"
    zcat "${INTERVALS}" | head -1
    zcat "${INTERVALS}" | tail -n +2 | sort -t$'\t' -k${sq_col} -rn | head -20
  else
    echo "  (squareness column not found, showing first 20 rows)"
    zcat "${INTERVALS}" | head -21
  fi

  echo ""
  echo "  Squareness distribution:"
  if [[ -n "${sq_col}" ]]; then
    zcat "${INTERVALS}" | tail -n +2 | awk -F'\t' -v c="${sq_col}" '{
      v=$c;
      if (v >= 0.6) bin=">=0.60 (strong)"
      else if (v >= 0.4) bin="0.40-0.59 (moderate)"
      else if (v >= 0.2) bin="0.20-0.39 (weak)"
      else bin="<0.20 (noise)"
      print bin
    }' | sort | uniq -c | sort -rn
  fi
fi

# ── Sub-regimes ──
SUBREG="${TRI_DIR}/triangle_subregimes.tsv.gz"
if [[ -f "${SUBREG}" ]]; then
  echo ""
  echo "── Sub-regimes ──"
  n_sub=$(zcat "${SUBREG}" | tail -n +2 | wc -l)
  echo "  Total sub-regimes: ${n_sub}"
  echo "  Header:"
  zcat "${SUBREG}" | head -1 | tr '\t' '\n' | nl
  echo "  First 5:"
  zcat "${SUBREG}" | head -6
fi

# ── Bridges ──
BRIDGES="${TRI_DIR}/triangle_bridges.tsv.gz"
if [[ -f "${BRIDGES}" ]]; then
  echo ""
  echo "── Bridges ──"
  n_br=$(zcat "${BRIDGES}" | tail -n +2 | wc -l)
  echo "  Total bridges: ${n_br}"

  echo "  By strength:"
  zcat "${BRIDGES}" | tail -n +2 | awk -F'\t' '{print $NF}' | sort | uniq -c | sort -rn | head -5
fi

# ── Window PA ──
PA="${TRI_DIR}/triangle_window_pa.tsv.gz"
if [[ -f "${PA}" ]]; then
  echo ""
  echo "── Window PA matrix ──"
  n_pa=$(zcat "${PA}" | tail -n +2 | wc -l)
  echo "  Total window rows: ${n_pa}"
  echo "  Header:"
  zcat "${PA}" | head -1 | tr '\t' '\n' | nl
fi

echo ""
echo "================================================================"
echo "  Done. Inspect specific intervals with:"
echo "    zcat ${INTERVALS} | awk -F'\\t' '\$NF==\"strong_triangle\"' | less -S"
echo "================================================================"
