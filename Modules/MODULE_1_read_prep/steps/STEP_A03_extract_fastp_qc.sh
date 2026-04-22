#!/usr/bin/env bash
###############################################################################
# STEP_A03_extract_fastp_qc.sh
#
# Parses fastp log files into a per-sample QC table and a compact numeric
# summary (min/mean/median/max) for key read-quality metrics.
#
# Usage:
#   bash STEP_A03_extract_fastp_qc.sh
#
# Outputs:
#   ${LOGDIR}/fastp_qc_table.tsv
#   ${LOGDIR}/fastp_qc_summary.tsv
#   ${LOGDIR}/STEP_A03_extract_fastp_qc.arg
#   ${LOGDIR}/STEP_A03_extract_fastp_qc.results
###############################################################################
set -euo pipefail

STEP="STEP_A03_extract_fastp_qc"
timestamp(){ date '+%F %T'; }

LOGDIR="/scratch/lt200308-agbsci/Quentin_project/00-samples/fastp/logs"
OUTTABLE="${LOGDIR}/fastp_qc_table.tsv"
OUTSUM="${LOGDIR}/fastp_qc_summary.tsv"
ARGFILE="${LOGDIR}/${STEP}.arg"
RESULTSFILE="${LOGDIR}/${STEP}.results"

cd "$LOGDIR"

NLOGFILES=$(ls -1 *.fastp.log 2>/dev/null | wc -l)
if [[ "$NLOGFILES" -eq 0 ]]; then
  echo "[ERROR] No *.fastp.log files found in $LOGDIR" >&2
  exit 1
fi

# ---- Write .arg ----
{
  echo -e "key\tvalue"
  echo -e "step\t${STEP}"
  echo -e "script_name\t$(basename "$0")"
  echo -e "datetime\t$(timestamp)"
  echo -e "host\t$(hostname)"
  echo -e "logdir\t${LOGDIR}"
  echo -e "n_log_files\t${NLOGFILES}"
  echo -e "output_table\t${OUTTABLE}"
  echo -e "output_summary\t${OUTSUM}"
} > "$ARGFILE"

# ---- 1) Build per-log table ----
awk -v OFS="\t" '
function pct_from_line(line,   m){
  if (match(line, /\(([0-9.]+)%\)/, m)) return m[1];
  return "";
}
BEGIN{
  print "ID","Sample","Run","Lane",
        "Adapter_R1","Adapter_R2",
        "R1_reads_before","R2_reads_before",
        "R1_Q20pct_before","R1_Q30pct_before","R2_Q20pct_before","R2_Q30pct_before",
        "R1_reads_after","R2_reads_after",
        "R1_Q20pct_after","R1_Q30pct_after","R2_Q20pct_after","R2_Q30pct_after",
        "Reads_passed","Fail_lowqual","Fail_toomanyN","Fail_tooshort",
        "Reads_adapter_trimmed","Bases_adapter_trimmed",
        "Reads_polyX","Bases_polyX",
        "Dup_rate_pct","Insert_peak"
}
FNR==1{
  id=FILENAME; sub(/\.fastp\.log$/,"",id)
  sample=run=lane=""
  split(id,a,"."); sample=a[1]; run=a[2]; lane=a[3]

  ad1=ad2=""
  r1rb=r2rb=r1ra=r2ra=""
  r1q20b=r1q30b=r2q20b=r2q30b=""
  r1q20a=r1q30a=r2q20a=r2q30a=""
  pass=flq=fn=fts=0
  r_ad=0; b_ad=0
  r_px=0; b_px=0
  dup=""; ins=""
}
/^No adapter detected for read1/ { ad1="No" }
/^No adapter detected for read2/ { ad2="No" }
/^Detected adapter for read1/     { ad1="Yes" }
/^Detected adapter for read2/     { ad2="Yes" }

/^Read1 before filtering:/ { mode="before"; which="R1"; next }
/^Read2 before filtering:/ { mode="before"; which="R2"; next }
/^Read1 after filtering:/  { mode="after";  which="R1"; next }
/^Read2 after filtering:/  { mode="after";  which="R2"; next }

/^[[:space:]]*total reads:/{
  gsub(/,/,"",$3)
  if (mode=="before" && which=="R1") r1rb=$3
  if (mode=="before" && which=="R2") r2rb=$3
  if (mode=="after"  && which=="R1") r1ra=$3
  if (mode=="after"  && which=="R2") r2ra=$3
}

/^[[:space:]]*Q20 bases:/{
  p=pct_from_line($0)
  if (mode=="before" && which=="R1") r1q20b=p
  if (mode=="before" && which=="R2") r2q20b=p
  if (mode=="after"  && which=="R1") r1q20a=p
  if (mode=="after"  && which=="R2") r2q20a=p
}
/^[[:space:]]*Q30 bases:/{
  p=pct_from_line($0)
  if (mode=="before" && which=="R1") r1q30b=p
  if (mode=="before" && which=="R2") r2q30b=p
  if (mode=="after"  && which=="R1") r1q30a=p
  if (mode=="after"  && which=="R2") r2q30a=p
}

/^reads passed filter:/                { gsub(/,/,"",$4); pass=$4 }
/^reads failed due to low quality:/    { gsub(/,/,"",$6); flq=$6 }
/^reads failed due to too many N:/     { gsub(/,/,"",$7); fn=$7 }
/^reads failed due to too short:/      { gsub(/,/,"",$6); fts=$6 }
/^reads with adapter trimmed:/         { gsub(/,/,"",$5); r_ad=$5 }
/^bases trimmed due to adapters:/      { gsub(/,/,"",$6); b_ad=$6 }
/^reads with polyX in 3'\'' end:/      { gsub(/,/,"",$7); r_px=$7 }
/^bases trimmed in polyX tail:/        { gsub(/,/,"",$7); b_px=$7 }
/^Duplication rate:/                   { dup=$3; gsub(/%/,"",dup) }
/^Insert size peak/                    { ins=$NF }

ENDFILE{
  if (ad1=="") ad1="NA"
  if (ad2=="") ad2="NA"
  print id,sample,run,lane,
        ad1,ad2,
        r1rb,r2rb,
        r1q20b,r1q30b,r2q20b,r2q30b,
        r1ra,r2ra,
        r1q20a,r1q30a,r2q20a,r2q30a,
        pass,flq,fn,fts,
        r_ad,b_ad,
        r_px,b_px,
        dup,ins
}
' *.fastp.log > "$OUTTABLE"

echo "[OK] Wrote table: $OUTTABLE ($(( $(wc -l < "$OUTTABLE") - 1 )) rows)"

# ---- 2) Summary statistics ----
python3 - "$OUTTABLE" "$OUTSUM" <<'PY'
import sys
import pandas as pd

table = pd.read_csv(sys.argv[1], sep="\t")
outpath = sys.argv[2]

cols = [
    "R1_Q30pct_before", "R2_Q30pct_before",
    "R1_Q30pct_after",  "R2_Q30pct_after",
    "Fail_toomanyN",
    "Reads_adapter_trimmed",
    "Dup_rate_pct",
    "Insert_peak",
]

for c in cols:
    table[c] = pd.to_numeric(table[c], errors="coerce")

rows = []
for c in cols:
    s = table[c].dropna()
    if s.empty:
        rows.append({"metric": c, "n": 0, "min": None, "mean": None, "median": None, "max": None})
    else:
        rows.append({
            "metric": c,
            "n": int(s.shape[0]),
            "min": float(s.min()),
            "mean": float(s.mean()),
            "median": float(s.median()),
            "max": float(s.max()),
        })

out = pd.DataFrame(rows)
out.to_csv(outpath, sep="\t", index=False)
print(f"[OK] Wrote summary: {outpath}")
PY

# ---- Write .results ----
{
  echo -e "key\tvalue\tdescription"
  echo -e "arg_file\t${ARGFILE}\tParameter record"
  echo -e "qc_table\t${OUTTABLE}\tPer-sample fastp QC metrics"
  echo -e "qc_summary\t${OUTSUM}\tSummary statistics across samples"
  echo -e "n_samples\t$(( $(wc -l < "$OUTTABLE") - 1 ))\tNumber of samples in QC table"
} > "$RESULTSFILE"

echo "[DONE] Wrote: $ARGFILE  $RESULTSFILE"
