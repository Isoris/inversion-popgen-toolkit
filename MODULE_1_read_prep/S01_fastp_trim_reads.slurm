#!/usr/bin/env bash
###############################################################################
# S01_fastp_trim_reads.slurm
#
# Re-runs fastp on a TSV-defined set of read pairs (typically redo cases
# such as unexpected EOF), using parallel per-sample jobs on a full node.
#
# Usage:
#   sbatch --export=ALL,TSV=/path/to/redo.tsv S01_fastp_trim_reads.slurm
#
# Input TSV columns (tab-separated, with header):
#   Sample  Run  Lane  R1  R2
#
# Outputs per sample:
#   ${OUTFQDIR}/<Sample>.<Run>.<Lane>.R1.fastp.fq.gz
#   ${OUTFQDIR}/<Sample>.<Run>.<Lane>.R2.fastp.fq.gz
#   ${LOGDIR}/<Sample>.<Run>.<Lane>.fastp.redo.log
#
# Sidecar files:
#   ${OUTFQDIR}/S01_fastp_trim_reads.arg
#   ${OUTFQDIR}/S01_fastp_trim_reads.results
#   ${OUTFQDIR}/S01_fastp_trim_reads.metric.tsv
###############################################################################
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -c 127
#SBATCH --mem=237GB
#SBATCH -t 1-00:00:00
#SBATCH -J fastp_redo
#SBATCH -o fastp_redo.%j.out
#SBATCH -e fastp_redo.%j.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

# ---- Paths ----
FASTP="$(command -v fastp)"
BASE="/scratch/lt200308-agbsci/Quentin_project"
INDIR="${BASE}/00-samples"
OUTFQDIR="${INDIR}/fastp"
LOGDIR="${OUTFQDIR}/logs"
mkdir -p "$OUTFQDIR" "$LOGDIR"

: "${TSV:?ERROR: set TSV in sbatch --export (TSV=/path/to/file.tsv)}"

# ---- Parallelism ----
CPUS="${SLURM_CPUS_PER_TASK:-1}"
TFASTP="${TFASTP:-8}"
NJOBS="${NJOBS:-$(( CPUS / TFASTP ))}"
(( NJOBS < 1 )) && NJOBS=1

# ---- fastp parameters ----
# These are the exact trimming parameters used for all WGS read prep.
# Changing any of these requires re-running the entire downstream pipeline.
FASTP_ARGS=(
  --detect_adapter_for_pe
  --trim_poly_g --trim_poly_x
  -5 -3 -q 20
  -f 5 -F 5 -t 5 -T 5
  -n 0 -l 30
)

timestamp(){ date '+%F %T'; }

# ---- Sidecar setup ----
STEP="S01_fastp_trim_reads"
ARGFILE="${OUTFQDIR}/${STEP}.arg"
RESULTSFILE="${OUTFQDIR}/${STEP}.results"
METRICFILE="${OUTFQDIR}/${STEP}.metric.tsv"
METRIC_TMPDIR="${OUTFQDIR}/.${STEP}.metric.parts"

SCRIPT_NAME="$(basename "$0")"
HOSTNAME_NOW="$(hostname)"
JOB_ID="${SLURM_JOB_ID:-NA}"
NOW="$(timestamp)"
FASTP_VERSION="$("$FASTP" --version 2>&1 | head -n 1 || echo NA)"

mkdir -p "$METRIC_TMPDIR"
rm -f "${METRIC_TMPDIR}"/*.tsv 2>/dev/null || true

# ---- Write .arg ----
{
  echo -e "key\tvalue"
  echo -e "step\t${STEP}"
  echo -e "script_name\t${SCRIPT_NAME}"
  echo -e "datetime\t${NOW}"
  echo -e "host\t${HOSTNAME_NOW}"
  echo -e "slurm_job_id\t${JOB_ID}"
  echo -e "tool_name\tfastp"
  echo -e "tool_path\t${FASTP}"
  echo -e "tool_version\t${FASTP_VERSION}"
  echo -e "base\t${BASE}"
  echo -e "input_dir\t${INDIR}"
  echo -e "output_fastq_dir\t${OUTFQDIR}"
  echo -e "logdir\t${LOGDIR}"
  echo -e "tsv\t${TSV}"
  echo -e "total_cpus\t${CPUS}"
  echo -e "threads_per_fastp\t${TFASTP}"
  echo -e "njobs\t${NJOBS}"
  echo -e "detect_adapter_for_pe\tyes"
  echo -e "trim_poly_g\tyes"
  echo -e "trim_poly_x\tyes"
  echo -e "quality_trim_5prime\tyes"
  echo -e "quality_trim_3prime\tyes"
  echo -e "qualified_quality_phred\t20"
  echo -e "trim_front_r1\t5"
  echo -e "trim_front_r2\t5"
  echo -e "trim_tail_r1\t5"
  echo -e "trim_tail_r2\t5"
  echo -e "max_n_bases\t0"
  echo -e "min_read_length\t30"
  echo -e "mode\tfastp_redo_from_tsv"
  echo -e "command_template\t${FASTP} -i <R1> -I <R2> -o <O1> -O <O2> ${FASTP_ARGS[*]} -w ${TFASTP}"
} > "$ARGFILE"

echo "[$(timestamp)] [INFO] host=${HOSTNAME_NOW} JOB_ID=${JOB_ID}"
echo "[$(timestamp)] [INFO] FASTP=$FASTP ($FASTP_VERSION)"
echo "[$(timestamp)] [INFO] TSV=$TSV"
echo "[$(timestamp)] [INFO] NJOBS=$NJOBS TFASTP=$TFASTP total_threads=$((NJOBS*TFASTP))"
echo "[$(timestamp)] [INFO] Wrote arg file: $ARGFILE"

# ---- Helpers ----
safe_size_bytes() {
  local f="$1"
  if [[ -e "$f" ]]; then stat -c %s "$f" 2>/dev/null || echo "NA"; else echo "NA"; fi
}

write_metric_row() {
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$1" "$2" "$3" "$4" "$5" "$6"
}

# ---- Per-sample worker ----
fastp_one() {
  local SAMPLE="$1" RUN="$2" LANE="$3" R1="$4" R2="$5"
  local ID="${SAMPLE}.${RUN}.${LANE}"
  local OUT1="${OUTFQDIR}/${ID}.R1.fastp.fq.gz"
  local OUT2="${OUTFQDIR}/${ID}.R2.fastp.fq.gz"
  local LOG="${LOGDIR}/${ID}.fastp.redo.log"
  local PART_METRIC="${METRIC_TMPDIR}/${ID}.tsv"
  local STATUS="FAILED"

  local IN1_BYTES IN2_BYTES
  IN1_BYTES="$(safe_size_bytes "$R1")"
  IN2_BYTES="$(safe_size_bytes "$R2")"

  {
    echo "[$(timestamp)] [START] $ID"
    echo "[$(timestamp)] [IN]  R1=$R1"
    echo "[$(timestamp)] [IN]  R2=$R2"
    echo "[$(timestamp)] [OUT] O1=$OUT1"
    echo "[$(timestamp)] [OUT] O2=$OUT2"

    [[ -s "$R1" ]] || { echo "[$(timestamp)] [ERROR] missing R1"; exit 2; }
    [[ -s "$R2" ]] || { echo "[$(timestamp)] [ERROR] missing R2"; exit 2; }

    rm -f "$OUT1" "$OUT2"

    "$FASTP" \
      -i "$R1" -I "$R2" \
      -o "$OUT1" -O "$OUT2" \
      "${FASTP_ARGS[@]}" \
      -w "$TFASTP"

    echo "[$(timestamp)] [DONE] OK"
  } > "$LOG" 2>&1 && STATUS="OK"

  local OUT1_BYTES OUT2_BYTES
  OUT1_BYTES="$(safe_size_bytes "$OUT1")"
  OUT2_BYTES="$(safe_size_bytes "$OUT2")"

  {
    write_metric_row "$ID" "fastp_redo" "input_r1_bytes" "$IN1_BYTES" "bytes" "R1_input_file_size"
    write_metric_row "$ID" "fastp_redo" "input_r2_bytes" "$IN2_BYTES" "bytes" "R2_input_file_size"
    write_metric_row "$ID" "fastp_redo" "output_r1_bytes" "$OUT1_BYTES" "bytes" "OUT1_fastp_file_size"
    write_metric_row "$ID" "fastp_redo" "output_r2_bytes" "$OUT2_BYTES" "bytes" "OUT2_fastp_file_size"
    write_metric_row "$ID" "fastp_redo" "output_r1_exists" "$([[ -s "$OUT1" ]] && echo yes || echo no)" "yesno" "test_-s"
    write_metric_row "$ID" "fastp_redo" "output_r2_exists" "$([[ -s "$OUT2" ]] && echo yes || echo no)" "yesno" "test_-s"
    write_metric_row "$ID" "fastp_redo" "redo_status" "$STATUS" "class" "fastp_exit_status"
  } > "$PART_METRIC"
}

export -f fastp_one timestamp safe_size_bytes write_metric_row
export FASTP TFASTP OUTFQDIR LOGDIR METRIC_TMPDIR
export FASTP_ARGS

# ---- Run ----
NROWS=$(( $(wc -l < "$TSV") - 1 ))
echo "[$(timestamp)] [INFO] rows(excl header)=$NROWS"

i=0
tail -n +2 "$TSV" | while IFS=$'\t' read -r SAMPLE RUN LANE R1 R2; do
  i=$((i+1))
  while (( $(jobs -rp | wc -l) >= NJOBS )); do sleep 1; done
  echo "[$(timestamp)] [QUEUE] ($i/$NROWS) ${SAMPLE}.${RUN}.${LANE}"
  fastp_one "$SAMPLE" "$RUN" "$LANE" "$R1" "$R2" &
done
wait

# ---- Merge metrics ----
{
  echo -e "sample\tmodule\tmetric\tvalue\tunit\tderived_from"
  find "$METRIC_TMPDIR" -maxdepth 1 -type f -name "*.tsv" -print0 \
    | sort -z | xargs -0 cat
} > "$METRICFILE"

# ---- Write .results ----
{
  echo -e "key\tvalue\tdescription"
  echo -e "arg_file\t${ARGFILE}\tParameter record"
  echo -e "metric_file\t${METRICFILE}\tPer-sample metric sidecar"
  echo -e "output_fastq_dir\t${OUTFQDIR}\tDirectory containing trimmed FASTQ pairs"
  echo -e "log_dir\t${LOGDIR}\tPer-sample fastp logs"
  echo -e "n_samples_processed\t${NROWS}\tTotal read pairs processed"
} > "$RESULTSFILE"

echo "[$(timestamp)] [DONE] Wrote: $ARGFILE  $METRICFILE  $RESULTSFILE"
