#!/usr/bin/env bash
###############################################################################
# STEP_A08_check_bam_bai_pairs.sh
#
# Checks that each .sr.bam has a matching .sr.bam.bai and vice versa.
# Reports missing BAM or BAI files.
#
# Usage:
#   cd /path/to/minimap2-bams && bash STEP_A08_check_bam_bai_pairs.sh
#
# Outputs:
#   bam_bai_check.tsv
#   STEP_A08_check_bam_bai_pairs.arg
#   STEP_A08_check_bam_bai_pairs.results
###############################################################################
set -euo pipefail

STEP="STEP_A08_check_bam_bai_pairs"
timestamp(){ date '+%F %T'; }

OUT="bam_bai_check.tsv"
ARGFILE="${STEP}.arg"
RESULTSFILE="${STEP}.results"

{
  echo -e "key\tvalue"
  echo -e "step\t${STEP}"
  echo -e "script_name\t$(basename "$0")"
  echo -e "datetime\t$(timestamp)"
  echo -e "host\t$(hostname)"
  echo -e "working_dir\t$(pwd)"
} > "$ARGFILE"

echo -e "status\tpath" > "$OUT"

n_missing_bai=0
n_missing_bam=0

find batch*/ -type f -name "*.sr.bam" ! -name "*.tmp.*" -print \
| sed 's/\.sr\.bam$//' \
| while read -r p; do
    if [[ ! -s "${p}.sr.bam.bai" ]]; then
      echo -e "MISSING_BAI\t${p}.sr.bam"
      n_missing_bai=$((n_missing_bai+1))
    fi
  done >> "$OUT"

find batch*/ -type f -name "*.sr.bam.bai" -print \
| sed 's/\.sr\.bam\.bai$//' \
| while read -r p; do
    if [[ ! -s "${p}.sr.bam" ]]; then
      echo -e "MISSING_BAM\t${p}.sr.bam.bai"
      n_missing_bam=$((n_missing_bam+1))
    fi
  done >> "$OUT"

n_issues=$(( $(wc -l < "$OUT") - 1 ))

{
  echo -e "key\tvalue\tdescription"
  echo -e "arg_file\t${ARGFILE}\tParameter record"
  echo -e "check_table\t${OUT}\tBAM/BAI consistency check"
  echo -e "n_issues\t${n_issues}\tTotal missing BAM or BAI files"
} > "$RESULTSFILE"

if [[ "$n_issues" -eq 0 ]]; then
  echo "[OK] All BAM/BAI pairs consistent"
else
  echo "[WARN] Found $n_issues issues — see $OUT"
fi
