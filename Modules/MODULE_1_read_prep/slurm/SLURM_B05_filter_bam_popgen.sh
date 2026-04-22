#!/usr/bin/env bash
###############################################################################
# SLURM_B05_filter_bam_popgen.slurm
#
# SLURM array job: filters merged markdup+clip BAMs for population genomics.
# Applies: MAPQ >= 60, proper pair, exclude unmapped/secondary/dup/supp,
# same-chromosome mate, absolute TLEN within empirical p99 band.
#
# Usage:
#   sbatch --array=0-<N-1> SLURM_B05_filter_bam_popgen.slurm
#
# Outputs per sample:
#   <sample>.filtered.bam + .bai
#   <sample>.filtered.flagstat.txt
#
# Sidecar:
#   ${RESULTS_DIR}/SLURM_B05_filter_bam_popgen.arg
###############################################################################
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=64GB
#SBATCH -t 1-00:00:00
#SBATCH -J bam_filter
#SBATCH -o bam_filter.%A_%a.out
#SBATCH -e bam_filter.%A_%a.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

STEP="SLURM_B05_filter_bam_popgen"
timestamp(){ date '+%F %T'; }

SAMTOOLS="$(command -v samtools)"
BASE="/scratch/lt200308-agbsci/Quentin_project"
IN_GLOB="${BASE}/02-merged_per_sample/*/*.merged.markdup.clip.bam"
RESULTS_DIR="${BASE}/02-merged_per_sample/results"
THREADS="${SLURM_NTASKS:-16}"

# ---- Filtering parameters ----
# These are the exact popgen BAM filters used for all downstream ANGSD work.
MINMAPQ=60
REMOVE_MASK=0xF0C               # unmapped + mate-unmapped + secondary + QCfail + dup + supp
REQUIRE_PROPER_PAIRS=1
ENFORCE_SAME_CHR=1
ENFORCE_TLEN=1
MIN_TLEN=0
MAX_TLEN=514                    # empirical p99(abs(TLEN)) = 514 bp

mkdir -p "$RESULTS_DIR"

# ---- Write .arg (once) ----
ARGFILE="${RESULTS_DIR}/${STEP}.arg"
if [[ ! -s "$ARGFILE" ]]; then
  {
    echo -e "key\tvalue"
    echo -e "step\t${STEP}"
    echo -e "script_name\t$(basename "$0")"
    echo -e "datetime\t$(timestamp)"
    echo -e "host\t$(hostname)"
    echo -e "slurm_job_id\t${SLURM_JOB_ID:-NA}"
    echo -e "samtools\t${SAMTOOLS} ($("$SAMTOOLS" --version | head -n1))"
    echo -e "min_mapq\t${MINMAPQ}"
    echo -e "remove_mask\t${REMOVE_MASK}"
    echo -e "require_proper_pairs\t${REQUIRE_PROPER_PAIRS}"
    echo -e "enforce_same_chr\t${ENFORCE_SAME_CHR}"
    echo -e "enforce_tlen\t${ENFORCE_TLEN}"
    echo -e "min_tlen\t${MIN_TLEN}"
    echo -e "max_tlen\t${MAX_TLEN}"
    echo -e "note_tlen\tEmpirical p99 from S11; p95=446"
    echo -e "threads\t${THREADS}"
  } > "$ARGFILE"
fi

# ---- Select this task's BAM ----
mapfile -t BAMS < <(ls -1 $IN_GLOB 2>/dev/null | sort)
NBAM="${#BAMS[@]}"
(( NBAM == 0 )) && { echo "[ERROR] No BAMs found: $IN_GLOB" >&2; exit 2; }

BAM="${BAMS[$SLURM_ARRAY_TASK_ID]}"
SDIR="$(dirname "$BAM")"
BN="$(basename "$BAM")"
SAMPLE="${BN%%.*}"

OUT="${SDIR}/${SAMPLE}.merged.markdup.clip.pp.samechr.tlenP99.filtered.bam"
LOG="${SDIR}/${SAMPLE}.filter.pipeline.log"

echo "[$(timestamp)] SAMPLE=$SAMPLE" | tee "$LOG"
echo "[$(timestamp)] BAM=$BAM" | tee -a "$LOG"
echo "[$(timestamp)] MAPQ>=${MINMAPQ} PP=${REQUIRE_PROPER_PAIRS} SAMECHR=${ENFORCE_SAME_CHR} TLEN=[${MIN_TLEN},${MAX_TLEN}]" | tee -a "$LOG"

[[ -s "$BAM" ]] || { echo "[ERROR] missing BAM" | tee -a "$LOG"; exit 3; }

# ---- Filter ----
VIEW_ARGS=(-@ "$THREADS" -h -F "$REMOVE_MASK" -q "$MINMAPQ")
(( REQUIRE_PROPER_PAIRS == 1 )) && VIEW_ARGS+=(-f 0x2)

"$SAMTOOLS" view "${VIEW_ARGS[@]}" "$BAM" \
| awk -v samechr="$ENFORCE_SAME_CHR" -v dotlen="$ENFORCE_TLEN" -v min="$MIN_TLEN" -v max="$MAX_TLEN" '
  BEGIN{OFS="\t"}
  /^@/ {print; next}
  {
    if (samechr==1 && !($7=="=" || $7==$3)) next
    if (dotlen==1) {
      t=$9; if (t<0) t=-t
      if (t<min || t>max) next
    }
    print
  }' \
| "$SAMTOOLS" view -@ "$THREADS" -bh -o "$OUT" - 2>>"$LOG"

"$SAMTOOLS" index -@ "$THREADS" "$OUT" >> "$LOG" 2>&1
"$SAMTOOLS" flagstat -@ "$THREADS" "$OUT" > "${SDIR}/${SAMPLE}.filtered.flagstat.txt"

echo "[$(timestamp)] [DONE] $SAMPLE" | tee -a "$LOG"
