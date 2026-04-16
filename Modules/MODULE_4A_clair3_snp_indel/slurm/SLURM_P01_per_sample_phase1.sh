#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=32G
#SBATCH -t 1-00:00:00
#SBATCH -A lt200308
#SBATCH -J c3post_v2
#SBATCH -o logs/c3post_v2.%A_%a.out
#SBATCH -e logs/c3post_v2.%A_%a.err
#SBATCH --array=1-226

set -euo pipefail

# ============================================================
# Per-sample STEP01-04 + STEP02B (phase) + STEP08
# v2: includes phase-aware blocks
# ============================================================

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
ROOT="${BASE}/MODULE_4A_SNP_INDEL50_Clair3"

SCRIPTS="${ROOT}/postprocess_scripts"
VCF_DIR="${ROOT}/vcf"
REF="${BASE}/00-samples/fClaHyb_Gar_LG.fa"
REF_FAI="${REF}.fai"
BED="${ROOT}/beds/C_gar_LG28.bed"
BAM_MANIFEST="${BASE}/pa_roary_results/00_manifests/sample_bam_minimap2_vs_P99TLENMAPQ30.tsv"
OUTBASE="${ROOT}/postprocess_results"
CHROM="C_gar_LG28"
LINE_NO=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --chrom)        CHROM="$2"; shift 2 ;;
        --vcf_dir)      VCF_DIR="$2"; shift 2 ;;
        --bed)          BED="$2"; shift 2 ;;
        --scripts)      SCRIPTS="$2"; shift 2 ;;
        --outbase)      OUTBASE="$2"; shift 2 ;;
        --bam_manifest) BAM_MANIFEST="$2"; shift 2 ;;
        --line)         LINE_NO="$2"; shift 2 ;;
        *) echo "Unknown: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "$LINE_NO" ]]; then
    LINE_NO="${SLURM_ARRAY_TASK_ID:-}"
fi
[[ -n "$LINE_NO" ]] || { echo "[ERROR] Need --line or SLURM_ARRAY_TASK_ID" >&2; exit 1; }

ROW=$(( LINE_NO + 1 ))
SAMPLE=$(awk -v r="$ROW" 'BEGIN{FS="\t"} NR==r{print $1}' "$BAM_MANIFEST")
BAM=$(awk -v r="$ROW" 'BEGIN{FS="\t"} NR==r{print $3}' "$BAM_MANIFEST")

[[ -n "$SAMPLE" ]] || { echo "[ERROR] No sample at line $LINE_NO" >&2; exit 1; }

echo "════════════════════════════════════════"
echo "[RUN v2] SAMPLE=$SAMPLE  CHROM=$CHROM  LINE=$LINE_NO"
echo "════════════════════════════════════════"

# Locate VCF
VCF=""
for pattern in \
    "${VCF_DIR}/${SAMPLE}.${CHROM}.clair3.discovery.vcf.gz" \
    "${VCF_DIR}/${SAMPLE}.chr28.clair3.discovery.vcf.gz" \
    "${VCF_DIR}/${SAMPLE}.clair3.vcf.gz" \
; do
    if [[ -s "$pattern" ]]; then VCF="$pattern"; break; fi
done

[[ -n "$VCF" ]] || { echo "[ERROR] No VCF for $SAMPLE in $VCF_DIR" >&2; exit 1; }
[[ -s "$BAM" ]] || { echo "[ERROR] BAM not found: $BAM" >&2; exit 1; }

OUTDIR="${OUTBASE}/${CHROM}/${SAMPLE}"
mkdir -p "$OUTDIR" "${ROOT}/logs"

echo "[STEP01] Parse & annotate …"
python "${SCRIPTS}/STEP01_parse_and_annotate_clair3_vcf.py" \
    --vcf "$VCF" --outdir "$OUTDIR" --sample "$SAMPLE" \
    --ref_fai "$REF_FAI" --bed "$BED"

echo "[STEP02] Local event blocks …"
python "${SCRIPTS}/STEP02_detect_local_event_blocks.py" \
    --annotated "${OUTDIR}/all_variants_annotated.tsv" --outdir "$OUTDIR"

echo "[STEP02B] Phase-aware blocks …"
python "${SCRIPTS}/STEP02B_phase_aware_blocks.py" \
    --vcf "$VCF" \
    --annotated "${OUTDIR}/all_variants_with_blocks.tsv" \
    --bam "$BAM" \
    --outdir "$OUTDIR" \
    --sample "$SAMPLE"

# Use phase-annotated input for rescue if available
INPUT03="${OUTDIR}/all_variants_with_blocks.tsv"
[[ -f "${OUTDIR}/all_variants_with_phase.tsv" ]] && INPUT03="${OUTDIR}/all_variants_with_phase.tsv"

echo "[STEP03] Strong single-sample rescue …"
python "${SCRIPTS}/STEP03_rescue_strong_single_sample.py" \
    --annotated "$INPUT03" --outdir "$OUTDIR"

echo "[STEP04] Weak indel candidates …"
python "${SCRIPTS}/STEP04_export_weak_indel_candidates.py" \
    --classified "${OUTDIR}/all_classified_step03.tsv" \
    --bam "$BAM" --ref "$REF" --outdir "$OUTDIR" --sample_id "$SAMPLE"

echo "[DONE v2] $SAMPLE → $OUTDIR"
