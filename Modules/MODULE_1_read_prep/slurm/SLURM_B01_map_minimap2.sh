#!/usr/bin/env bash
###############################################################################
# SLURM_B01_map_minimap2.slurm
#
# Maps fastp-processed paired-end reads with minimap2 in SR mode, adds
# read-group information, sorts and indexes each per-run BAM.
#
# Usage:
#   sbatch --export=ALL,BATCH=00 SLURM_B01_map_minimap2.slurm
#   sbatch --export=ALL,BATCH=01,TSV=/path/to/batch01.tsv SLURM_B01_map_minimap2.slurm
#
# Input TSV columns (tab-separated, with header):
#   Sample  Run  Lane  R1  R2
#
# Outputs per run-lane:
#   ${OUTDIR}/<Sample>.<Run>.<Lane>.sr.bam
#   ${OUTDIR}/<Sample>.<Run>.<Lane>.sr.bam.bai
#   ${LOGDIR}/<Sample>.<Run>.<Lane>.log
#
# Sidecar files:
#   ${OUTDIR}/SLURM_B01_map_minimap2.arg
#   ${OUTDIR}/SLURM_B01_map_minimap2.results
###############################################################################
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 127
#SBATCH --mem=237GB
#SBATCH -t 5-00:00:00
#SBATCH -J map_mm2
#SBATCH -o map_mm2.%j.out
#SBATCH -e map_mm2.%j.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

STEP="SLURM_B01_map_minimap2"
timestamp(){ date '+%F %T'; }

# ---- Paths ----
MINIMAP2="$(command -v minimap2)"
SAMTOOLS="$(command -v samtools)"

BASE="/scratch/lt200308-agbsci/Quentin_project"
INDIR="${BASE}/00-samples"
REF="${INDIR}/fClaHyb_Gar_LG.fa"

: "${BATCH:=00}"
TSV="${TSV:-${INDIR}/cga_fastq.batch.${BATCH}.tsv}"
OUTDIR="${BASE}/01-minimap2-bams/batch${BATCH}"
LOGDIR="${OUTDIR}/logs"
mkdir -p "$OUTDIR" "$LOGDIR"

# ---- Parallelism ----
NJOBS="${NJOBS:-8}"
TMAP="${TMAP:-16}"

# ---- Minimap2 parameters ----
# -ax sr     : short-read preset
# -Y         : use soft clipping for supplementary
# -c         : output CIGAR in PAF (ignored in SAM but harmless)
# --secondary-seq  : output sequence for secondary alignments
# --sam-hit-only   : omit unmapped reads from SAM output
MM2_ARGS=(-ax sr -Y -c --secondary-seq --sam-hit-only)

ARGFILE="${OUTDIR}/${STEP}.arg"
RESULTSFILE="${OUTDIR}/${STEP}.results"

MM2_VERSION="$("$MINIMAP2" --version 2>/dev/null || echo NA)"
SAM_VERSION="$("$SAMTOOLS" --version | head -n 1 || echo NA)"

# ---- Write .arg ----
{
  echo -e "key\tvalue"
  echo -e "step\t${STEP}"
  echo -e "script_name\t$(basename "$0")"
  echo -e "datetime\t$(timestamp)"
  echo -e "host\t$(hostname)"
  echo -e "slurm_job_id\t${SLURM_JOB_ID:-NA}"
  echo -e "minimap2_path\t${MINIMAP2}"
  echo -e "minimap2_version\t${MM2_VERSION}"
  echo -e "samtools_path\t${SAMTOOLS}"
  echo -e "samtools_version\t${SAM_VERSION}"
  echo -e "reference\t${REF}"
  echo -e "batch\t${BATCH}"
  echo -e "tsv\t${TSV}"
  echo -e "outdir\t${OUTDIR}"
  echo -e "njobs\t${NJOBS}"
  echo -e "threads_per_map\t${TMAP}"
  echo -e "mm2_preset\tsr"
  echo -e "mm2_args\t${MM2_ARGS[*]}"
} > "$ARGFILE"

echo "[$(timestamp)] [INFO] BATCH=$BATCH NJOBS=$NJOBS TMAP=$TMAP"
echo "[$(timestamp)] [INFO] REF=$REF"
echo "[$(timestamp)] [INFO] TSV=$TSV"

# ---- Build index once ----
if [[ ! -s "${REF}.mmi" ]]; then
  echo "[$(timestamp)] Building minimap2 index..."
  "$MINIMAP2" -d "${REF}.mmi" "$REF"
fi

# ---- Per-run worker ----
map_one() {
  local SAMPLE="$1" RUN="$2" LANE="$3" R1="$4" R2="$5"

  # Use fastp outputs (override raw paths from TSV)
  local FASTP_DIR="${INDIR}/fastp"
  R1="${FASTP_DIR}/${SAMPLE}.${RUN}.${LANE}.R1.fastp.fq.gz"
  R2="${FASTP_DIR}/${SAMPLE}.${RUN}.${LANE}.R2.fastp.fq.gz"

  local SM="$SAMPLE"
  local ID="${SAMPLE}.${RUN}.${LANE}"
  local LB="$SAMPLE"
  local PU="${RUN}.${LANE}"
  local OUTBAM="${OUTDIR}/${ID}.sr.bam"
  local ROWLOG="${LOGDIR}/${ID}.log"

  {
    echo "[$(timestamp)] [START] ID=$ID"
    echo "[$(timestamp)] [RG] @RG\\tID:${ID}\\tSM:${SM}\\tLB:${LB}\\tPU:${PU}\\tPL:ILLUMINA"

    if [[ -s "$OUTBAM" && -s "${OUTBAM}.bai" ]]; then
      echo "[$(timestamp)] [SKIP] BAM+BAI already exist"
      exit 0
    fi

    "$MINIMAP2" "${MM2_ARGS[@]}" \
      -t "$TMAP" \
      -R "@RG\\tID:${ID}\\tSM:${SM}\\tLB:${LB}\\tPU:${PU}\\tPL:ILLUMINA" \
      "${REF}.mmi" "$R1" "$R2" \
    | "$SAMTOOLS" sort -@ "$TMAP" -o "$OUTBAM" -

    "$SAMTOOLS" index -@ "$TMAP" "$OUTBAM"
    echo "[$(timestamp)] [DONE] OK"
  } > "$ROWLOG" 2>&1
}

export -f map_one timestamp
export MINIMAP2 SAMTOOLS REF OUTDIR LOGDIR TMAP INDIR
export MM2_ARGS

NROWS=$(( $(wc -l < "$TSV") - 1 ))
echo "[$(timestamp)] TSV rows=$NROWS"

i=0
tail -n +2 "$TSV" | while IFS=$'\t' read -r SAMPLE RUN LANE R1 R2; do
  i=$((i+1))
  while (( $(jobs -rp | wc -l) >= NJOBS )); do sleep 3; done
  echo "[$(timestamp)] [QUEUE] ($i/$NROWS) ${SAMPLE}.${RUN}.${LANE}"
  map_one "$SAMPLE" "$RUN" "$LANE" "$R1" "$R2" &
done
wait

# ---- Write .results ----
n_bams=$(find "$OUTDIR" -name "*.sr.bam" ! -name "*.tmp.*" | wc -l)
{
  echo -e "key\tvalue\tdescription"
  echo -e "arg_file\t${ARGFILE}\tParameter record"
  echo -e "bam_dir\t${OUTDIR}\tPer-run BAM output directory"
  echo -e "log_dir\t${LOGDIR}\tPer-run mapping logs"
  echo -e "n_bams_produced\t${n_bams}\tTotal BAMs in output dir"
  echo -e "n_input_rows\t${NROWS}\tRows in input TSV"
} > "$RESULTSFILE"

echo "[$(timestamp)] [DONE] All mappings finished for batch${BATCH} ($n_bams BAMs)"
