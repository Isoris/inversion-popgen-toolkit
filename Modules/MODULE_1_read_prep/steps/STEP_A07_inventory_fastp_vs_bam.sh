#!/usr/bin/env bash
###############################################################################
# STEP_A07_inventory_fastp_vs_bam.sh
#
# Compares fastp output pairs against generated .sr.bam files and their
# indexes, then reports whether each sample is complete, missing, or
# inconsistent.
#
# Usage:
#   bash STEP_A07_inventory_fastp_vs_bam.sh
#
# Outputs:
#   fastp_vs_bam.tsv
#   STEP_A07_inventory_fastp_vs_bam.arg
#   STEP_A07_inventory_fastp_vs_bam.results
###############################################################################
set -euo pipefail

STEP="STEP_A07_inventory_fastp_vs_bam"
timestamp(){ date '+%F %T'; }

FASTP_DIR="/scratch/lt200308-agbsci/Quentin_project/00-samples/fastp"
BAMROOT="/scratch/lt200308-agbsci/Quentin_project/01-minimap2-bams"
OUT="fastp_vs_bam.tsv"
ARGFILE="${STEP}.arg"
RESULTSFILE="${STEP}.results"

# ---- Write .arg ----
{
  echo -e "key\tvalue"
  echo -e "step\t${STEP}"
  echo -e "script_name\t$(basename "$0")"
  echo -e "datetime\t$(timestamp)"
  echo -e "host\t$(hostname)"
  echo -e "fastp_dir\t${FASTP_DIR}"
  echo -e "bam_root\t${BAMROOT}"
  echo -e "output\t${OUT}"
} > "$ARGFILE"

echo -e "id\tfastp_R1\tfastp_R2\tfastp_pair_ok\tbam_path\tbam_ok\tbai_ok\tstatus" > "$OUT"

# Collect fastp IDs (must have R1)
tmp_fastp=$(mktemp)
tmp_bam=$(mktemp)
tmp_all=$(mktemp)
trap 'rm -f "$tmp_fastp" "$tmp_bam" "$tmp_all"' EXIT

find "$FASTP_DIR" -maxdepth 1 -type f -name "*.R1.fastp.fq.gz" | sort \
  | sed 's#.*/##; s/\.R1\.fastp\.fq\.gz$//' > "$tmp_fastp"

find "$BAMROOT"/batch*/ -type f -name "*.sr.bam" ! -name "*.tmp.*" | sort \
  | sed 's#.*/##; s/\.sr\.bam$//' > "$tmp_bam"

cat "$tmp_fastp" "$tmp_bam" | sort -u > "$tmp_all"

n_ok=0; n_issues=0
while read -r id; do
  r1="${FASTP_DIR}/${id}.R1.fastp.fq.gz"
  r2="${FASTP_DIR}/${id}.R2.fastp.fq.gz"
  r1_ok=0; [[ -s "$r1" ]] && r1_ok=1
  r2_ok=0; [[ -s "$r2" ]] && r2_ok=1
  pair_ok=0; (( r1_ok == 1 && r2_ok == 1 )) && pair_ok=1

  bam="$(find "$BAMROOT"/batch*/ -type f -name "${id}.sr.bam" ! -name "*.tmp.*" | head -n 1)"
  bam_ok=0; [[ -n "$bam" && -s "$bam" ]] && bam_ok=1
  bai_ok=0; (( bam_ok == 1 )) && [[ -s "${bam}.bai" ]] && bai_ok=1

  status="OK"
  if (( pair_ok == 1 && bam_ok == 0 )); then status="FASTP_OK_BUT_NO_BAM"
  elif (( pair_ok == 0 && bam_ok == 1 )); then status="BAM_OK_BUT_NO_FASTP_PAIR"
  elif (( bam_ok == 1 && bai_ok == 0 )); then status="BAM_OK_BUT_MISSING_BAI"
  elif (( pair_ok == 0 && bam_ok == 0 )); then status="NEITHER_PRESENT"
  fi

  [[ "$status" == "OK" ]] && n_ok=$((n_ok+1)) || n_issues=$((n_issues+1))

  echo -e "${id}\t${r1}\t${r2}\t${pair_ok}\t${bam:-NA}\t${bam_ok}\t${bai_ok}\t${status}"
done < "$tmp_all" >> "$OUT"

# ---- Write .results ----
{
  echo -e "key\tvalue\tdescription"
  echo -e "arg_file\t${ARGFILE}\tParameter record"
  echo -e "inventory\t${OUT}\tFastp vs BAM inventory comparison"
  echo -e "n_total\t$(wc -l < "$tmp_all")\tTotal unique IDs"
  echo -e "n_ok\t${n_ok}\tIDs with both fastp pair and BAM+BAI"
  echo -e "n_issues\t${n_issues}\tIDs with missing or inconsistent files"
} > "$RESULTSFILE"

echo "[DONE] Wrote: $OUT (OK=$n_ok  issues=$n_issues)"
