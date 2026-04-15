#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --mem=180G
#SBATCH -t 3-00:00:00
#SBATCH -A lt200308
#SBATCH -J clair3_disc
#SBATCH -o logs/clair3_disc.%A_%a.out
#SBATCH -e logs/clair3_disc.%A_%a.err

set -euo pipefail

# ============================================================
# 01_run_clair3_discovery.sh
#
# Clair3 variant discovery for one sample on one chromosome.
# Designed for SLURM array submission: one job per sample.
#
# Clair3 parameters (identical to the validated chr28 pilot):
#   --qual=20
#   --snp_min_af=0.08
#   --indel_min_af=0.08
#   --min_mq=20
#   --min_coverage=2
#   --gvcf
#   --include_all_ctgs
#   --use_whatshap_for_intermediate_phasing   ← phase blocks
#   --use_whatshap_for_final_output_phasing   ← phased GT + PS tags
#   --platform=ilmn
#   --model_path=/opt/models/ilmn
#
# Phase output:
#   WhatsHap produces phased genotypes (pipe-separated 0|1) and PS tags
#   in the final VCF. These feed directly into STEP02B of post-processing.
#
# Usage (per-chromosome array):
#   sbatch --array=1-226 01_run_clair3_discovery.sh --chrom C_gar_LG28
#
# Usage (single sample):
#   bash 01_run_clair3_discovery.sh --chrom C_gar_LG28 --line 5
#
# Usage (explicit BAM):
#   bash 01_run_clair3_discovery.sh --chrom C_gar_LG28 \
#       --bam /path/to/sample.bam --sample SAMPLENAME
# ============================================================

THREADS="${SLURM_NTASKS:-64}"

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
ROOT="${BASE}/MODULE_4A_SNP_INDEL50_Clair3"
CLAIR3_SIF="${ROOT}/clair3_latest.sif"

REF="${BASE}/00-samples/fClaHyb_Gar_LG.fa"
MANIFEST="${ROOT}/meta/clair3_sample_manifest.tsv"
CHROM=""
BAM=""
SAMPLE=""
LINE_NO=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --chrom|--ctg)  CHROM="$2"; shift 2 ;;
        --bam)          BAM="$(realpath "$2")"; shift 2 ;;
        --sample)       SAMPLE="$2"; shift 2 ;;
        --manifest)     MANIFEST="$(realpath "$2")"; shift 2 ;;
        --line)         LINE_NO="$2"; shift 2 ;;
        --ref)          REF="$(realpath "$2")"; shift 2 ;;
        --threads)      THREADS="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

# ── Validate chromosome ──────────────────────────────────────────────────────
[[ -n "$CHROM" ]] || { echo "[ERROR] --chrom is required" >&2; exit 1; }

BED="${ROOT}/beds/${CHROM}.bed"
[[ -s "$BED" ]] || { echo "[ERROR] Missing BED for $CHROM: $BED" >&2; exit 1; }

# ── Resolve sample/BAM ──────────────────────────────────────────────────────
if [[ -z "$BAM" ]]; then
    [[ -s "$MANIFEST" ]] || { echo "[ERROR] Need --bam or valid --manifest" >&2; exit 1; }

    if [[ -z "$LINE_NO" ]]; then
        LINE_NO="${SLURM_ARRAY_TASK_ID:-}"
    fi
    [[ -n "$LINE_NO" ]] || { echo "[ERROR] Need --line or SLURM_ARRAY_TASK_ID" >&2; exit 1; }

    ROW=$(( LINE_NO + 1 ))  # +1 for header
    SAMPLE="$(awk -v r="$ROW" 'BEGIN{FS="\t"} NR==r{print $1}' "$MANIFEST")"
    BAM="$(awk -v r="$ROW" 'BEGIN{FS="\t"} NR==r{print $2}' "$MANIFEST")"

    [[ -n "$SAMPLE" ]] || { echo "[ERROR] No sample at manifest line $LINE_NO" >&2; exit 1; }
    [[ -n "$BAM" ]]    || { echo "[ERROR] No BAM at manifest line $LINE_NO" >&2; exit 1; }
fi

[[ -s "$BAM" ]] || { echo "[ERROR] BAM not found: $BAM" >&2; exit 1; }

if [[ -z "$SAMPLE" ]]; then
    SAMPLE="$(basename "$BAM" .bam | sed 's/\.merged.*//; s/\.markdup.*//')"
fi

# ── Environment setup ────────────────────────────────────────────────────────
source ~/.bashrc
mamba activate /lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly
module load Apptainer

# ── Validate inputs ──────────────────────────────────────────────────────────
[[ -s "$CLAIR3_SIF" ]] || { echo "[ERROR] Missing Clair3 container: $CLAIR3_SIF" >&2; exit 1; }
[[ -s "$REF" ]]        || { echo "[ERROR] Missing REF: $REF" >&2; exit 1; }
[[ -f "${REF}.fai" ]]  || samtools faidx "$REF"

# Ensure BAM index
if [[ ! -f "${BAM}.bai" && ! -f "${BAM%.bam}.bai" ]]; then
    echo "[INFO] Indexing BAM: $BAM"
    samtools index "$BAM"
fi

# ── Output directories ───────────────────────────────────────────────────────
mkdir -p "${ROOT}/logs" "${ROOT}/runs" \
         "${ROOT}/vcf/${CHROM}" "${ROOT}/gvcf/${CHROM}" "${ROOT}/meta"

RUN_DIR="${ROOT}/runs/${SAMPLE}_${CHROM}_clair3"
rm -rf "$RUN_DIR"
mkdir -p "$RUN_DIR"

# ── Build Clair3 command ─────────────────────────────────────────────────────
CMD=(
  /opt/bin/run_clair3.sh
  --bam_fn="$BAM"
  --ref_fn="$REF"
  --threads="$THREADS"
  --platform="ilmn"
  --model_path="/opt/models/ilmn"
  --output="$RUN_DIR"

  # ── Quality / sensitivity thresholds ──
  --qual=20
  --snp_min_af=0.08
  --indel_min_af=0.08
  --min_mq=20
  --min_coverage=2

  # ── Output modes ──
  --gvcf

  # ── Contig handling ──
  --include_all_ctgs
  --bed_fn="$BED"
  --ctg_name="$CHROM"

  # ── Phasing ──
  --use_whatshap_for_intermediate_phasing
  --use_whatshap_for_final_output_phasing
)

# ── Run ──────────────────────────────────────────────────────────────────────
echo "════════════════════════════════════════════════════════"
echo "[CLAIR3] SAMPLE  = $SAMPLE"
echo "[CLAIR3] CHROM   = $CHROM"
echo "[CLAIR3] BAM     = $BAM"
echo "[CLAIR3] REF     = $REF"
echo "[CLAIR3] BED     = $BED"
echo "[CLAIR3] THREADS = $THREADS"
echo "[CLAIR3] RUN_DIR = $RUN_DIR"
echo "════════════════════════════════════════════════════════"

START_TIME=$(date +%s)

apptainer exec \
  -B "$(dirname "$BAM"):$(dirname "$BAM")" \
  -B "$(dirname "$REF"):$(dirname "$REF")" \
  -B "$ROOT:$ROOT" \
  -B "$(dirname "$BED"):$(dirname "$BED")" \
  "$CLAIR3_SIF" \
  "${CMD[@]}"

END_TIME=$(date +%s)
ELAPSED=$(( END_TIME - START_TIME ))

echo "[CLAIR3] Runtime: ${ELAPSED}s ($(( ELAPSED / 60 ))m)"

# ── Locate and copy final VCF ────────────────────────────────────────────────
# Prefer phased merged output if present, else fall back in strict order.

FINAL_VCF=""
FINAL_VCF_TYPE=""

for candidate in \
    "${RUN_DIR}/phased_merge_output.vcf.gz:PHASED_MERGED" \
    "${RUN_DIR}/merge_output.vcf.gz:MERGED" \
    "${RUN_DIR}/full_alignment.vcf.gz:FULL_ALIGNMENT" \
    "${RUN_DIR}/pileup.vcf.gz:PILEUP"
do
    FILE="${candidate%%:*}"
    TYPE="${candidate##*:}"
    if [[ -s "$FILE" ]]; then
        FINAL_VCF="$FILE"
        FINAL_VCF_TYPE="$TYPE"
        break
    fi
done

if [[ -z "$FINAL_VCF" ]]; then
    FALLBACK_VCF="$(find "$RUN_DIR" -maxdepth 2 -type f -name '*.vcf.gz' \
        ! -name '*.tbi' ! -name '*.csi' \
        | grep -E 'phased_merge_output|merge_output|full_alignment|pileup' \
        | head -n 1 2>/dev/null || true)"
    if [[ -n "$FALLBACK_VCF" && -s "$FALLBACK_VCF" ]]; then
        FINAL_VCF="$FALLBACK_VCF"
        FINAL_VCF_TYPE="FALLBACK"
    fi
fi

[[ -n "$FINAL_VCF" ]] || {
    echo "[ERROR] No final VCF found in $RUN_DIR" >&2
    find "$RUN_DIR" -type f \( -name "*.vcf.gz" -o -name "*.gvcf.gz" \) | sort >&2
    exit 1
}

OUT_VCF="${ROOT}/vcf/${CHROM}/${SAMPLE}.${CHROM}.clair3.discovery.vcf.gz"
cp -f "$FINAL_VCF" "$OUT_VCF"
bcftools index -f "$OUT_VCF"

echo "[CLAIR3] VCF source type → $FINAL_VCF_TYPE"
echo "[CLAIR3] VCF source file → $FINAL_VCF"
echo "[CLAIR3] VCF copied to   → $OUT_VCF"

# ── Phase audit for copied VCF ───────────────────────────────────────────────
HAS_PS_HEADER=0
N_PHASED_GT=0
N_TOTAL_RECORDS=0
FOUND_PHASED_MERGED=0

if [[ "$FINAL_VCF_TYPE" == "PHASED_MERGED" ]]; then
    FOUND_PHASED_MERGED=1
fi

if bcftools view -h "$OUT_VCF" | grep -q '^##FORMAT=<ID=PS,'; then
    HAS_PS_HEADER=1
fi

N_PHASED_GT=$(bcftools view -H "$OUT_VCF" | awk '$10 ~ /\|/' | wc -l)
N_TOTAL_RECORDS=$(bcftools view -H "$OUT_VCF" | wc -l)

PHASE_AUDIT_TMP="${ROOT}/meta/phase_audit.discovery.${SAMPLE}.${CHROM}.tmp.tsv"
echo -e "SAMPLE\tCHROM\tFOUND_PHASED_MERGED\tVCF_SOURCE_TYPE\tHAS_PS_HEADER\tN_PHASED_GT\tN_TOTAL_RECORDS\tVCF_PATH" > "$PHASE_AUDIT_TMP"
echo -e "${SAMPLE}\t${CHROM}\t${FOUND_PHASED_MERGED}\t${FINAL_VCF_TYPE}\t${HAS_PS_HEADER}\t${N_PHASED_GT}\t${N_TOTAL_RECORDS}\t${OUT_VCF}" >> "$PHASE_AUDIT_TMP"

echo "[CLAIR3] Phase audit: found_phased_merge=${FOUND_PHASED_MERGED}, PS_header=${HAS_PS_HEADER}, phased_GT=${N_PHASED_GT}, total_records=${N_TOTAL_RECORDS}"

# ── Locate and copy final GVCF ───────────────────────────────────────────────
# Clair3 usually provides merge_output.gvcf.gz. Prefer phased gVCF if it
# ever exists, else fall back to merged gVCF.

FINAL_GVCF=""
FINAL_GVCF_TYPE=""

for candidate in \
    "${RUN_DIR}/phased_merge_output.gvcf.gz:PHASED_MERGED_GVCF" \
    "${RUN_DIR}/phased_merge_output.gvcf.vcf.gz:PHASED_MERGED_GVCF_ALT" \
    "${RUN_DIR}/merge_output.gvcf.gz:MERGED_GVCF" \
    "${RUN_DIR}/merge_output.gvcf.vcf.gz:MERGED_GVCF_ALT"
do
    FILE="${candidate%%:*}"
    TYPE="${candidate##*:}"
    if [[ -s "$FILE" ]]; then
        FINAL_GVCF="$FILE"
        FINAL_GVCF_TYPE="$TYPE"
        break
    fi
done

if [[ -z "$FINAL_GVCF" ]]; then
    FALLBACK_GVCF="$(find "$RUN_DIR" -maxdepth 2 -type f \
        \( -name '*.gvcf.gz' -o -name '*gvcf.vcf.gz' \) \
        ! -name '*.tbi' ! -name '*.csi' \
        | grep -E 'phased_merge_output|merge_output|gvcf' \
        | head -n 1 2>/dev/null || true)"
    if [[ -n "$FALLBACK_GVCF" && -s "$FALLBACK_GVCF" ]]; then
        FINAL_GVCF="$FALLBACK_GVCF"
        FINAL_GVCF_TYPE="FALLBACK_GVCF"
    fi
fi

if [[ -n "$FINAL_GVCF" ]]; then
    OUT_GVCF="${ROOT}/gvcf/${CHROM}/${SAMPLE}.${CHROM}.clair3.discovery.gvcf.gz"
    cp -f "$FINAL_GVCF" "$OUT_GVCF"
    bcftools index -f "$OUT_GVCF"
    echo "[CLAIR3] GVCF source type → $FINAL_GVCF_TYPE"
    echo "[CLAIR3] GVCF source file → $FINAL_GVCF"
    echo "[CLAIR3] GVCF copied to   → $OUT_GVCF"
fi

# ── Cleanup intermediate files (optional, saves disk) ───────────────────────
# Uncomment to remove the working directory after successful copy:
# rm -rf "$RUN_DIR"

echo ""
echo "[DONE] Clair3 discovery complete"
echo "[DONE] SAMPLE=$SAMPLE  CHROM=$CHROM"
echo "[DONE] VCF:  $OUT_VCF"
[[ -n "$FINAL_GVCF" ]] && echo "[DONE] GVCF: $OUT_GVCF"
echo "[DONE] Runtime: ${ELAPSED}s"
