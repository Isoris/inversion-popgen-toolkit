#!/usr/bin/env bash
###############################################################################
# STEP_A02_fastp_eof_redo_table.sh
#
# Finds fastp outputs with unexpected EOF errors, extracts affected sample IDs,
# searches uploadku for matching raw read pairs, and writes redo/missing tables.
#
# Usage:
#   bash STEP_A02_fastp_eof_redo_table.sh
#
# Outputs:
#   ${REPORT_DIR}/fastp_unexpected_eof_redo.tsv
#   ${REPORT_DIR}/fastp_unexpected_eof_notfound_in_uploadku.tsv
#   ${REPORT_DIR}/STEP_A02_fastp_eof_redo_table.arg
#   ${REPORT_DIR}/STEP_A02_fastp_eof_redo_table.results
###############################################################################
set -euo pipefail

STEP="STEP_A02_fastp_eof_redo_table"
timestamp(){ date '+%F %T'; }

# ---- Paths ----
REPORT_DIR="/scratch/lt200308-agbsci/Quentin_project/00-samples/fastp/gzip_check_report_4139502_20260116_151011"
GZIP_TSV="${REPORT_DIR}/gzip_check_all.tsv"
FASTP_DIR="/scratch/lt200308-agbsci/Quentin_project/00-samples/fastp"
UPLOADKU_GLOB="/project/*agb*/uploadku"

OUT_TSV="${REPORT_DIR}/fastp_unexpected_eof_redo.tsv"
MISS_TSV="${REPORT_DIR}/fastp_unexpected_eof_notfound_in_uploadku.tsv"
ARGFILE="${REPORT_DIR}/${STEP}.arg"
RESULTSFILE="${REPORT_DIR}/${STEP}.results"

# ---- Sanity ----
[[ -s "$GZIP_TSV" ]] || { echo "[ERROR] Missing $GZIP_TSV" >&2; exit 1; }
[[ -d "$FASTP_DIR" ]] || { echo "[ERROR] Missing $FASTP_DIR" >&2; exit 1; }

# ---- Write .arg ----
{
  echo -e "key\tvalue"
  echo -e "step\t${STEP}"
  echo -e "script_name\t$(basename "$0")"
  echo -e "datetime\t$(timestamp)"
  echo -e "host\t$(hostname)"
  echo -e "gzip_tsv\t${GZIP_TSV}"
  echo -e "fastp_dir\t${FASTP_DIR}"
  echo -e "uploadku_glob\t${UPLOADKU_GLOB}"
  echo -e "out_tsv\t${OUT_TSV}"
  echo -e "miss_tsv\t${MISS_TSV}"
} > "$ARGFILE"

# ---- Collect uploadku dirs ----
mapfile -t UPDIRS < <(ls -d ${UPLOADKU_GLOB} 2>/dev/null || true)
if [[ "${#UPDIRS[@]}" -eq 0 ]]; then
  echo "[ERROR] No uploadku dirs found via: ${UPLOADKU_GLOB}" >&2
  exit 1
fi
echo "[INFO] Found ${#UPDIRS[@]} uploadku dirs"

# ---- Helper: find raw R1/R2 by ID ----
find_raw_pair() {
  local id="$1" r1="" r2=""
  for d in "${UPDIRS[@]}"; do
    r1=$(find "$d" -type f \
        \( -name "*${id}*R1*.f*q.gz" -o -name "*${id}*_R1*.f*q.gz" \
           -o -name "*${id}*1.f*q.gz" -o -name "*${id}*_1.f*q.gz" \) \
        2>/dev/null | head -n 1 || true)
    r2=$(find "$d" -type f \
        \( -name "*${id}*R2*.f*q.gz" -o -name "*${id}*_R2*.f*q.gz" \
           -o -name "*${id}*2.f*q.gz" -o -name "*${id}*_2.f*q.gz" \) \
        2>/dev/null | head -n 1 || true)
    if [[ -n "$r1" && -n "$r2" ]]; then
      echo -e "${r1}\t${r2}"
      return 0
    fi
  done
  echo -e "\t"
}

# ---- Extract EOF IDs ----
mapfile -t IDS < <(
  grep -F "unexpected end of file" "$GZIP_TSV" \
  | awk '{print $1}' \
  | xargs -n1 basename \
  | sed -E 's/\.R[12]\.fastp\.fq\.gz$//' \
  | sort -u
)
echo "[INFO] Unexpected-EOF IDs: ${#IDS[@]}"

# ---- Write outputs ----
echo -e "id\tfastp_R1\tfastp_R2\traw_R1_uploadku\traw_R2_uploadku\tstatus" > "$OUT_TSV"
echo -e "id\tfastp_R1\tfastp_R2\tnote" > "$MISS_TSV"

ok=0; miss=0
for id in "${IDS[@]}"; do
  f1="${FASTP_DIR}/${id}.R1.fastp.fq.gz"
  f2="${FASTP_DIR}/${id}.R2.fastp.fq.gz"
  raw_pair="$(find_raw_pair "$id")"
  raw1="$(cut -f1 <<<"$raw_pair")"
  raw2="$(cut -f2 <<<"$raw_pair")"

  if [[ -n "$raw1" && -n "$raw2" ]]; then
    echo -e "${id}\t${f1}\t${f2}\t${raw1}\t${raw2}\tFOUND_RAW_PAIR" >> "$OUT_TSV"
    ok=$((ok+1))
  else
    echo -e "${id}\t${f1}\t${f2}\tRAW_NOT_FOUND_IN_UPLOADKU" >> "$OUT_TSV"
    echo -e "${id}\t${f1}\t${f2}\tRAW_NOT_FOUND_IN_UPLOADKU" >> "$MISS_TSV"
    miss=$((miss+1))
  fi
done

# ---- Write .results ----
{
  echo -e "key\tvalue\tdescription"
  echo -e "arg_file\t${ARGFILE}\tParameter record"
  echo -e "redo_tsv\t${OUT_TSV}\tSamples needing fastp redo with raw read paths"
  echo -e "missing_tsv\t${MISS_TSV}\tSamples where raw reads were not found"
  echo -e "n_eof_ids\t${#IDS[@]}\tTotal IDs with unexpected EOF"
  echo -e "n_found_raw\t${ok}\tIDs with raw pairs found in uploadku"
  echo -e "n_missing_raw\t${miss}\tIDs with no raw pairs in uploadku"
} > "$RESULTSFILE"

echo "[DONE] Wrote: $OUT_TSV"
echo "[DONE] Missing: $MISS_TSV"
echo "[SUMMARY] FOUND_RAW_PAIR=$ok  RAW_NOT_FOUND=$miss"
