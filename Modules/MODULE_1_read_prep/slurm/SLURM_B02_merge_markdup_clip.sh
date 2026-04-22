#!/usr/bin/env bash
###############################################################################
# SLURM_B02_merge_markdup_clip.slurm
#
# SLURM array job: merges per-run BAMs into one BAM per sample, validates
# RG/SM consistency, performs fixmate + markdup, optionally runs clipOverlap,
# generates QC outputs (flagstat, samtools stats), and appends a summary row.
#
# Usage:
#   # First run once to build sample list, then:
#   sbatch --array=1-<N> SLURM_B02_merge_markdup_clip.slurm
#
# Outputs per sample:
#   ${OUTBASE}/<Sample>/<Sample>.merged.markdup.clip.bam
#   ${OUTBASE}/<Sample>/<Sample>.flagstat.txt
#   ${OUTBASE}/<Sample>/<Sample>.samtools.stats.txt
#   ${OUTBASE}/results/QC_ROWS.tsv  (one row appended per sample)
#
# Sidecar (written once if missing):
#   ${OUTBASE}/SLURM_B02_merge_markdup_clip.arg
###############################################################################
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=120GB
#SBATCH -t 2-00:00:00
#SBATCH -J merge_md_clip
#SBATCH -o merge_md_clip.%A_%a.out
#SBATCH -e merge_md_clip.%A_%a.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

STEP="SLURM_B02_merge_markdup_clip"
timestamp(){ date '+%F %T'; }

# ---- Config ----
BAMBASE="/scratch/lt200308-agbsci/Quentin_project/01-minimap2-bams"
OUTBASE="/scratch/lt200308-agbsci/Quentin_project/02-merged_per_sample"
REF="/scratch/lt200308-agbsci/Quentin_project/00-samples/fClaHyb_Gar_LG.fa"

THREADS="${SLURM_NTASKS:-16}"
DO_CLIP="${DO_CLIP:-1}"       # 1=clipOverlap on, 0=off

SAMTOOLS="$(command -v samtools)"
BAMUTIL="$(command -v bam)"   # BamUtil clipOverlap

mkdir -p "$OUTBASE"/{lists,logs,tmp,results}

# ---- Safe sort memory ----
mem_mb="${SLURM_MEM_PER_NODE:-0}"
(( mem_mb == 0 )) && mem_mb=$(awk '/MemTotal/ {print int($2/1024)}' /proc/meminfo)
usable_mb=$(( mem_mb * 80 / 100 ))
m_per_thread=$(( usable_mb / THREADS ))
(( m_per_thread < 256 )) && m_per_thread=256
(( m_per_thread > 1500 )) && m_per_thread=1500
SORT_MEM_PER_THREAD="${m_per_thread}M"

# ---- Build input lists (once) ----
ALLBAMS="${OUTBASE}/lists/all_run_bams.txt"
if [[ ! -s "$ALLBAMS" ]]; then
  find "$BAMBASE"/batch* -type f -name "*.sr.bam" ! -name "*.tmp.*" | sort > "$ALLBAMS"
fi

SAMPLES="${OUTBASE}/lists/samples.txt"
if [[ ! -s "$SAMPLES" ]]; then
  awk -F'/' '{print $NF}' "$ALLBAMS" \
    | sed 's/\.sr\.bam$//' \
    | awk -F'.' '{print $1}' \
    | sort -u > "$SAMPLES"
fi

NSAMP=$(wc -l < "$SAMPLES")
[[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || { echo "[ERROR] Run as: sbatch --array=1-${NSAMP} $0"; exit 1; }

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES")
[[ -n "$SAMPLE" ]] || { echo "[ERROR] Empty SAMPLE for task ${SLURM_ARRAY_TASK_ID}"; exit 1; }

# ---- Write .arg (once) ----
ARGFILE="${OUTBASE}/${STEP}.arg"
if [[ ! -s "$ARGFILE" ]]; then
  {
    echo -e "key\tvalue"
    echo -e "step\t${STEP}"
    echo -e "script_name\t$(basename "$0")"
    echo -e "datetime\t$(timestamp)"
    echo -e "host\t$(hostname)"
    echo -e "slurm_job_id\t${SLURM_JOB_ID:-NA}"
    echo -e "samtools\t${SAMTOOLS} ($("$SAMTOOLS" --version | head -n 1))"
    echo -e "bamutil\t${BAMUTIL}"
    echo -e "reference\t${REF}"
    echo -e "bam_base\t${BAMBASE}"
    echo -e "out_base\t${OUTBASE}"
    echo -e "threads\t${THREADS}"
    echo -e "sort_mem_per_thread\t${SORT_MEM_PER_THREAD}"
    echo -e "do_clip\t${DO_CLIP}"
    echo -e "n_samples\t${NSAMP}"
    echo -e "pipeline\tmerge -> namesort -> fixmate -> coordsort -> markdup -> clipOverlap -> index -> flagstat -> stats"
  } > "$ARGFILE"
fi

# ---- Per-sample file set ----
SDIR="${OUTBASE}/${SAMPLE}"
mkdir -p "$SDIR"

RUNLIST="${SDIR}/${SAMPLE}.run_bams.txt"
grep -F "/${SAMPLE}." "$ALLBAMS" > "$RUNLIST" || true
[[ -s "$RUNLIST" ]] || { echo "[ERROR] No run BAMs for $SAMPLE"; exit 2; }

TMPD="${OUTBASE}/tmp/${SAMPLE}.tmp_${SLURM_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPD"
export TMPDIR="$TMPD"

LOG="${SDIR}/${SAMPLE}.pipeline.log"
FINAL_BAM="${SDIR}/${SAMPLE}.merged.markdup${DO_CLIP:+.clip}.bam"
FINAL_BAI="${FINAL_BAM}.bai"
QC_ROW="${OUTBASE}/results/QC_ROWS.tsv"

# Skip if done
if [[ -s "$FINAL_BAM" && -s "$FINAL_BAI" ]]; then
  echo "[$(timestamp)] [SKIP] $SAMPLE already finished" | tee -a "$LOG"
  exit 0
fi

echo "[$(timestamp)] [INFO] SAMPLE=$SAMPLE task=${SLURM_ARRAY_TASK_ID}/${NSAMP}" | tee -a "$LOG"
echo "[$(timestamp)] [INFO] THREADS=$THREADS sort_mem=$SORT_MEM_PER_THREAD $(wc -l < "$RUNLIST") run BAMs" | tee -a "$LOG"

# ---- 0) Validate RG SM tags ----
echo "[$(timestamp)] [STEP] validate @RG SM tags" | tee -a "$LOG"
BADSM=0
while read -r bam; do
  sm_list=$("$SAMTOOLS" view -H "$bam" | awk -F'\t' '
    $1=="@RG"{ for(i=1;i<=NF;i++) if($i~/^SM:/) {sub(/^SM:/,"",$i); print $i} }' | sort -u)
  [[ -z "$sm_list" ]] && { echo "[ERROR] No @RG SM in $bam" | tee -a "$LOG"; BADSM=1; continue; }
  while read -r sm; do
    [[ "$sm" != "$SAMPLE" ]] && { echo "[ERROR] SM mismatch: $sm != $SAMPLE in $bam" | tee -a "$LOG"; BADSM=1; }
  done <<< "$sm_list"
done < "$RUNLIST"
(( BADSM != 0 )) && { echo "[FAIL] Fix RG first" | tee -a "$LOG"; exit 3; }

# ---- 1) Merge ----
MERGED="${TMPD}/${SAMPLE}.merged.bam"
echo "[$(timestamp)] [STEP] merge" | tee -a "$LOG"
"$SAMTOOLS" merge -@ "$THREADS" -f -b "$RUNLIST" -o "$MERGED"

# ---- 2) namesort -> fixmate -> coordsort ----
NSORT="${TMPD}/${SAMPLE}.namesort.bam"
FIXM="${TMPD}/${SAMPLE}.fixmate.bam"
CSORT="${TMPD}/${SAMPLE}.possort.bam"

echo "[$(timestamp)] [STEP] name-sort" | tee -a "$LOG"
"$SAMTOOLS" sort -n -@ "$THREADS" -m "$SORT_MEM_PER_THREAD" -T "${TMPD}/nsort" -o "$NSORT" "$MERGED"

echo "[$(timestamp)] [STEP] fixmate" | tee -a "$LOG"
"$SAMTOOLS" fixmate -m -@ "$THREADS" "$NSORT" "$FIXM"

echo "[$(timestamp)] [STEP] coord-sort" | tee -a "$LOG"
"$SAMTOOLS" sort -@ "$THREADS" -m "$SORT_MEM_PER_THREAD" -T "${TMPD}/csort" -o "$CSORT" "$FIXM"

# ---- 3) markdup ----
DUP="${TMPD}/${SAMPLE}.markdup.bam"
echo "[$(timestamp)] [STEP] markdup (mark only)" | tee -a "$LOG"
"$SAMTOOLS" markdup -@ "$THREADS" -s "$CSORT" "$DUP" 2> "${SDIR}/${SAMPLE}.markdup.stderr.txt"

# ---- 4) clipOverlap ----
if [[ "$DO_CLIP" -eq 1 ]]; then
  echo "[$(timestamp)] [STEP] clipOverlap" | tee -a "$LOG"
  "$BAMUTIL" clipOverlap --in "$DUP" --out "$FINAL_BAM" --storeOrig OC --stats --unmapped \
    > "${SDIR}/${SAMPLE}.clipOverlap.stdout.txt" 2> "${SDIR}/${SAMPLE}.clipOverlap.stats.txt"
else
  mv "$DUP" "$FINAL_BAM"
fi

# ---- 5) QC ----
echo "[$(timestamp)] [STEP] index + quickcheck + flagstat + stats" | tee -a "$LOG"
"$SAMTOOLS" index -@ "$THREADS" "$FINAL_BAM"

QC_STATUS="OK"
"$SAMTOOLS" quickcheck -v "$FINAL_BAM" > "${SDIR}/${SAMPLE}.quickcheck.txt" 2>&1 || QC_STATUS="BAD"
"$SAMTOOLS" flagstat -@ "$THREADS" "$FINAL_BAM" > "${SDIR}/${SAMPLE}.flagstat.txt"
"$SAMTOOLS" stats -@ "$THREADS" "$FINAL_BAM" > "${SDIR}/${SAMPLE}.samtools.stats.txt"

# Parse quick metrics
total_reads=$(awk '/in total/ {print $1; exit}' "${SDIR}/${SAMPLE}.flagstat.txt" || echo NA)
mapped_reads=$(awk '/ mapped \(/ && $0 !~ /primary/ {print $1; exit}' "${SDIR}/${SAMPLE}.flagstat.txt" || echo NA)
pct_mapped=$(awk '/ mapped \(/ && $0 !~ /primary/ {gsub(/[()%]/,""); print $5; exit}' "${SDIR}/${SAMPLE}.flagstat.txt" || echo NA)
ins_mean=$(awk -F'\t' '$1=="SN" && $2~/(insert size average)/ {print $3; exit}' "${SDIR}/${SAMPLE}.samtools.stats.txt" || echo NA)
err_rate=$(awk -F'\t' '$1=="SN" && $2~/(error rate)/ {print $3; exit}' "${SDIR}/${SAMPLE}.samtools.stats.txt" || echo NA)

# Append QC row (thread-safe via flock)
mkdir -p "$(dirname "$QC_ROW")"
{
  flock -x 9
  if [[ ! -s "$QC_ROW" ]]; then
    printf "bam\tsample\tn_run_bams\tstatus\tthreads\tsort_mem\tn_reads\tn_mapped\tpct_mapped\tinsert_mean\terror_rate\n" > "$QC_ROW"
  fi
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$FINAL_BAM" "$SAMPLE" "$(wc -l < "$RUNLIST")" "$QC_STATUS" \
    "$THREADS" "$SORT_MEM_PER_THREAD" \
    "$total_reads" "$mapped_reads" "$pct_mapped" "$ins_mean" "$err_rate" >> "$QC_ROW"
} 9>>"${OUTBASE}/results/.lock"

echo "[$(timestamp)] [DONE] $SAMPLE status=$QC_STATUS" | tee -a "$LOG"
rm -rf "$TMPD"
