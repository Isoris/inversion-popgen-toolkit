#!/usr/bin/env bash
###############################################################################
# extract_module1_review_metrics.sh
#
# Run this on the HPC to extract all quantitative metrics requested by the
# Module 1 peer review. Produces a single TSV that can be pasted directly
# into the revised methods.
#
# Usage:
#   bash extract_module1_review_metrics.sh > module1_review_metrics.tsv
#
# Requires: the Module 1 outputs to exist at their expected paths.
###############################################################################
set -euo pipefail

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
FASTP_LOGDIR="${BASE}/00-samples/fastp/logs"
FASTP_DIR="${BASE}/00-samples/fastp"
MERGED_DIR="${BASE}/02-merged_per_sample"
FILTERED_DIR="${BASE}/03-popgen_filtered"     # adjust if different
BAMROOT="${BASE}/01-minimap2-bams"
MASH_DIR="${BASE}/00-samples/species_assign_mash_fastp"
GENOME_SIZE=963905721

echo "============================================================"
echo "MODULE 1 PEER REVIEW — QUANTITATIVE METRICS EXTRACTION"
echo "Date: $(date '+%F %T')"
echo "============================================================"
echo ""

# =========================================================================
# M1 §1.1: Fastp QC summary
# =========================================================================
echo "=== §1.1 FASTP QC ==="

# Count total fastp log files (= run-lane units processed)
N_LOGS=$(ls -1 "${FASTP_LOGDIR}"/*.fastp.log 2>/dev/null | wc -l)
echo "run_lane_units_processed: ${N_LOGS}"

# If fastp_qc_summary.tsv exists (from S03), extract key metrics
QC_SUM="${FASTP_LOGDIR}/fastp_qc_summary.tsv"
if [[ -s "$QC_SUM" ]]; then
  echo "fastp_qc_summary: ${QC_SUM}"
  cat "$QC_SUM"
else
  echo "[TODO] Run S03_extract_fastp_qc.sh first to produce fastp_qc_summary.tsv"
fi

# Count EOF redo cases (from S02)
REDO_DIR="${BASE}/00-samples/fastp/gzip_check_report_4139502_20260116_151011"
if [[ -d "$REDO_DIR" ]]; then
  REDO_TSV=$(find "$REDO_DIR" -name "fastp_unexpected_eof_redo.tsv" | head -1)
  if [[ -s "$REDO_TSV" ]]; then
    N_REDO=$(tail -n +2 "$REDO_TSV" | wc -l)
    echo "eof_redo_samples: ${N_REDO}"
  fi
fi

echo ""

# =========================================================================
# m1: Actual read length — check from fastp logs
# =========================================================================
echo "=== m1 READ LENGTH ==="
FIRST_LOG=$(ls -1 "${FASTP_LOGDIR}"/*.fastp.log 2>/dev/null | head -1)
if [[ -s "$FIRST_LOG" ]]; then
  grep -i "read.length\|read1_mean_length\|Read1 before\|mean length" "$FIRST_LOG" | head -5
fi
echo ""

# =========================================================================
# m3: Run-lane units per sample
# =========================================================================
echo "=== m3 RUN-LANE UNITS PER SAMPLE ==="
ls -1 "${FASTP_DIR}"/*.R1.fastp.fq.gz 2>/dev/null \
  | sed 's#.*/##; s/\.[^.]*\.[^.]*\.R1\.fastp\.fq\.gz$//' \
  | sort | uniq -c | sort -rn \
  | awk '
    BEGIN { n=0; sum=0 }
    { counts[NR]=$1; n++; sum+=$1 }
    END {
      asort(counts)
      printf "n_samples: %d\n", n
      printf "min_runlanes: %d\n", counts[1]
      printf "max_runlanes: %d\n", counts[n]
      printf "mean_runlanes: %.1f\n", sum/n
      mid = int(n/2)+1
      printf "median_runlanes: %d\n", counts[mid]
    }
  '
echo ""

# =========================================================================
# M1 §1.2: Species assignment results
# =========================================================================
echo "=== §1.2 SPECIES ASSIGNMENT ==="
MASH_METRIC="${MASH_DIR}/S04_species_assign_mash.metric.tsv"
if [[ -s "$MASH_METRIC" ]]; then
  echo "mash_metric_file: ${MASH_METRIC}"
  # Count assignments
  awk -F'\t' 'NR>1 {print $NF}' "$MASH_METRIC" | sort | uniq -c | sort -rn
  echo ""
  # Count per-CGA consensus
  echo "Per-CGA consensus:"
  awk -F'\t' 'NR>1 {cga[$2]=$NF} END {for(c in cga) print cga[c]}' "$MASH_METRIC" \
    | sort | uniq -c | sort -rn
else
  echo "[TODO] Mash metric file not found at ${MASH_METRIC}"
fi
echo ""

# =========================================================================
# M1 §1.3: Mapping rates from flagstat
# =========================================================================
echo "=== §1.3 MAPPING RATES ==="
N_FLAGSTAT=0
TOTAL_MAPPED_PCT=0
for fs in "${MERGED_DIR}"/*/*.flagstat.txt; do
  [[ -s "$fs" ]] || continue
  pct=$(grep "mapped" "$fs" | head -1 | grep -oP '\([\d.]+%' | tr -d '(%')
  if [[ -n "$pct" ]]; then
    TOTAL_MAPPED_PCT=$(awk "BEGIN{print $TOTAL_MAPPED_PCT + $pct}")
    N_FLAGSTAT=$((N_FLAGSTAT + 1))
  fi
done
if (( N_FLAGSTAT > 0 )); then
  MEAN_MAP=$(awk "BEGIN{printf \"%.2f\", $TOTAL_MAPPED_PCT / $N_FLAGSTAT}")
  echo "n_samples_with_flagstat: ${N_FLAGSTAT}"
  echo "mean_mapping_rate_pct: ${MEAN_MAP}"
fi
echo ""

# =========================================================================
# M1 §1.4: Duplicate rates from markdup stderr
# =========================================================================
echo "=== §1.4 DUPLICATE RATES ==="
N_DUP=0
TOTAL_DUP_PCT=0
for md in "${MERGED_DIR}"/*/*.markdup.stderr.txt; do
  [[ -s "$md" ]] || continue
  # samtools markdup -s prints "READ X WRITTEN Y EXCLUDED Z (A%)"
  pct=$(grep -oP 'DUPLICATE.*?[\d.]+%' "$md" | grep -oP '[\d.]+(?=%)' | head -1)
  if [[ -n "$pct" ]]; then
    TOTAL_DUP_PCT=$(awk "BEGIN{print $TOTAL_DUP_PCT + $pct}")
    N_DUP=$((N_DUP + 1))
  fi
done
if (( N_DUP > 0 )); then
  MEAN_DUP=$(awk "BEGIN{printf \"%.2f\", $TOTAL_DUP_PCT / $N_DUP}")
  echo "n_samples_with_dup_stats: ${N_DUP}"
  echo "mean_duplicate_rate_pct: ${MEAN_DUP}"
else
  # Alternative: from S15 QC table
  QC_TABLE="${MERGED_DIR}/bam_qc_summary/per_sample_qc.tsv"
  if [[ -s "$QC_TABLE" ]]; then
    echo "Extracting from S15 QC table..."
    awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)if($i~"dup")col=i} NR>1&&col{sum+=$col;n++} END{printf "mean_dup_rate: %.2f%%\nn_samples: %d\n", sum/n, n}' "$QC_TABLE"
  else
    echo "[TODO] No markdup stderr or S15 QC table found"
  fi
fi
echo ""

# =========================================================================
# M1 §1.5: TLEN p99 scope — was it from one sample or all?
# =========================================================================
echo "=== §1.5 TLEN P99 SCOPE ==="
# Check if S11 was run on all or one
TLEN_FILES=$(find "${MERGED_DIR}" -name "*tlen*" -o -name "*S11*" 2>/dev/null | head -5)
if [[ -n "$TLEN_FILES" ]]; then
  echo "TLEN output files found:"
  echo "$TLEN_FILES"
  N_TLEN=$(echo "$TLEN_FILES" | wc -l)
  echo "n_tlen_files: ${N_TLEN}"
else
  echo "[INFO] Check manually: was S11 run on one sample or all 226?"
fi
echo ""

# =========================================================================
# M1 §1.7: Depth distribution from S15 or mosdepth
# =========================================================================
echo "=== §1.7 DEPTH DISTRIBUTION ==="
QC_TABLE_S15=$(find "${MERGED_DIR}" -name "per_sample_qc.tsv" -o -name "bam_qc_table.tsv" 2>/dev/null | head -1)
if [[ -s "$QC_TABLE_S15" ]]; then
  echo "S15 QC table: ${QC_TABLE_S15}"
  # Extract depth column
  awk -F'\t' '
    NR==1 { for(i=1;i<=NF;i++) if($i ~ /depth/) dcol=i; next }
    dcol && $dcol+0>0 { depths[++n]=$dcol+0; sum+=$dcol }
    END {
      if(n==0){print "[WARN] No depth data"; exit}
      asort(depths)
      printf "n_samples: %d\n", n
      printf "mean_depth: %.1fx\n", sum/n
      printf "median_depth: %.1fx\n", depths[int(n/2)+1]
      printf "min_depth: %.1fx\n", depths[1]
      printf "max_depth: %.1fx\n", depths[n]
      printf "p10_depth: %.1fx\n", depths[int(n*0.1)+1]
      printf "p90_depth: %.1fx\n", depths[int(n*0.9)]
    }' "$QC_TABLE_S15"
else
  echo "[TODO] Run S15_summarize_bam_qc.sh first"
fi
echo ""

# =========================================================================
# m8: Final sample count entering Module 2
# =========================================================================
echo "=== m8 FINAL SAMPLE COUNT ==="
SAMPLE_LIST="${BASE}/popstruct_thin/list_of_samples_one_per_line_same_bamfile_list.tsv"
if [[ -s "$SAMPLE_LIST" ]]; then
  N_FINAL=$(wc -l < "$SAMPLE_LIST")
  echo "final_samples_entering_module2: ${N_FINAL}"
else
  echo "[TODO] Sample list not found"
fi

BAMLIST="${BASE}/bamlist.pp.samechr.tlenP99.filtered.txt"
if [[ -s "$BAMLIST" ]]; then
  N_BAMS=$(wc -l < "$BAMLIST")
  echo "final_filtered_bams: ${N_BAMS}"
fi
echo ""

# =========================================================================
# M3: MAPQ 60 impact — fraction of reads lost
# =========================================================================
echo "=== M3 MAPQ60 IMPACT ==="
echo "[INFO] To compute: take one representative merged BAM and count:"
echo "  total_mapped=\$(samtools view -c -F 4 sample.bam)"
echo "  mapq60=\$(samtools view -c -F 4 -q 60 sample.bam)"
echo "  mapq30=\$(samtools view -c -F 4 -q 30 sample.bam)"
echo "  echo \"mapq60_retention: \$((mapq60 * 100 / total_mapped))%\""
echo "  echo \"mapq30_retention: \$((mapq30 * 100 / total_mapped))%\""
echo ""
echo "Or use the S14 MQ0/MQ30 tracking per chromosome:"
MQ_TABLE=$(find "${BASE}" -name "*mapq*" -o -name "*mq_dist*" 2>/dev/null | head -1)
if [[ -n "$MQ_TABLE" ]]; then
  echo "Found: $MQ_TABLE"
fi
echo ""

# =========================================================================
# Tool versions
# =========================================================================
echo "=== TOOL VERSIONS ==="
for tool in fastp minimap2 samtools mosdepth bedtools mash meryl; do
  v=$(command -v $tool 2>/dev/null && $tool --version 2>&1 | head -1 || echo "not found")
  echo "${tool}: ${v}"
done
# BamUtil
BAM=$(command -v bam 2>/dev/null || echo "")
if [[ -n "$BAM" ]]; then
  echo "BamUtil: $($BAM 2>&1 | head -1 || echo 'version unknown')"
fi
echo ""

echo "============================================================"
echo "DONE. Fill [TODO] placeholders in methods with these values."
echo "============================================================"
