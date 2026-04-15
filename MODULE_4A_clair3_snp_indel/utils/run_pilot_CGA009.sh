#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# UPGRADED pilot run: CGA009 chr28
# Runs STEP01-04 + STEP02B (phase) + STEP08 (classification)
#      + STEP09 (marker package) + STEP10 (figures)
# ============================================================

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
ROOT="${BASE}/MODULE_4A_SNP_INDEL50_Clair3"
SCRIPTS="${ROOT}/postprocess_scripts"

VCF="${ROOT}/vcf/CGA009.chr28.clair3.discovery.vcf.gz"
BAM="${BASE}/02-merged_per_sample/CGA009/CGA009.merged.markdup.clip.pp.samechr.tlenP99.filtered.bam"
REF="${BASE}/00-samples/fClaHyb_Gar_LG.fa"
REF_FAI="${REF}.fai"
BED="${ROOT}/beds/C_gar_LG28.bed"
SAMPLE="CGA009"
CHROM="C_gar_LG28"

OUTDIR="${ROOT}/postprocess_results/${CHROM}/${SAMPLE}"
FIGDIR="${OUTDIR}/figures"
ZOODIR="${OUTDIR}/01_Use_for_markers_multiplex"
mkdir -p "$OUTDIR" "$FIGDIR" "$ZOODIR"

echo "════════════════════════════════════════════════════════"
echo "UPGRADED PILOT: $SAMPLE  $CHROM"
echo "════════════════════════════════════════════════════════"

# ─── STEP01: Parse & annotate VCF ───
echo ""
echo "──── STEP01: Parse & annotate VCF ────"
python "${SCRIPTS}/STEP01_parse_and_annotate_clair3_vcf.py" \
    --vcf "$VCF" \
    --outdir "$OUTDIR" \
    --sample "$SAMPLE" \
    --ref_fai "$REF_FAI" \
    --bed "$BED"

# ─── STEP02: Local event blocks ───
echo ""
echo "──── STEP02: Local event blocks ────"
python "${SCRIPTS}/STEP02_detect_local_event_blocks.py" \
    --annotated "${OUTDIR}/all_variants_annotated.tsv" \
    --outdir "$OUTDIR"

# ─── STEP02B: Phase-aware blocks (NEW) ───
echo ""
echo "──── STEP02B: Phase-aware local blocks ────"
python "${SCRIPTS}/STEP02B_phase_aware_blocks.py" \
    --vcf "$VCF" \
    --annotated "${OUTDIR}/all_variants_with_blocks.tsv" \
    --bam "$BAM" \
    --outdir "$OUTDIR" \
    --sample "$SAMPLE"

# ─── STEP03: Strong single-sample rescue ───
echo ""
echo "──── STEP03: Strong single-sample rescue ────"
# Use phase-annotated input if available, else fall back
if [[ -f "${OUTDIR}/all_variants_with_phase.tsv" ]]; then
    INPUT_FOR_STEP03="${OUTDIR}/all_variants_with_phase.tsv"
else
    INPUT_FOR_STEP03="${OUTDIR}/all_variants_with_blocks.tsv"
fi
python "${SCRIPTS}/STEP03_rescue_strong_single_sample.py" \
    --annotated "$INPUT_FOR_STEP03" \
    --outdir "$OUTDIR"

# ─── STEP04: Weak indel candidates + read evidence ───
echo ""
echo "──── STEP04: Weak indel candidates + read evidence ────"
python "${SCRIPTS}/STEP04_export_weak_indel_candidates.py" \
    --classified "${OUTDIR}/all_classified_step03.tsv" \
    --bam "$BAM" \
    --ref "$REF" \
    --outdir "$OUTDIR" \
    --sample_id "$SAMPLE"

# ─── STEP08: Final classification (single-sample, no population yet) ───
echo ""
echo "──── STEP08: Final classification ────"
python "${SCRIPTS}/STEP08_final_classification.py" \
    --step03 "${OUTDIR}/all_classified_step03.tsv" \
    --weak "${OUTDIR}/weak_indel_candidates.tsv" \
    --blocks "${OUTDIR}/local_overlap_blocks.tsv" \
    --outdir "$OUTDIR" \
    --sample_id "$SAMPLE"

# ─── STEP09: Marker handoff package (NEW) ───
echo ""
echo "──── STEP09: Marker handoff package ────"
PHASE_ARG=""
if [[ -f "${OUTDIR}/all_variants_with_phase.tsv" ]]; then
    PHASE_ARG="--phase ${OUTDIR}/all_variants_with_phase.tsv"
fi

python "${SCRIPTS}/STEP09_marker_handoff_package.py" \
    --classified "${OUTDIR}/final_variant_classification.tsv" \
    $PHASE_ARG \
    --ref "$REF" \
    --outdir "$ZOODIR" \
    --sample_id "$SAMPLE" \
    --flank_bp 150

# ─── STEP10: Publication figure (NEW) ───
echo ""
echo "──── STEP10: Publication figure ────"
# Get chromosome length from FAI
CHR_LEN=$(awk -v c="$CHROM" '$1==c{print $2}' "$REF_FAI")
[[ -n "$CHR_LEN" ]] || CHR_LEN=33300000

Rscript "${SCRIPTS}/STEP10_publication_figure.R" \
    "${OUTDIR}/final_variant_classification.tsv" \
    "${OUTDIR}/weak_indel_candidates.tsv" \
    "$FIGDIR" \
    "$CHR_LEN" \
    "$SAMPLE"

echo ""
echo "════════════════════════════════════════════════════════"
echo "[UPGRADED PILOT DONE]"
echo ""
echo "Core outputs:"
ls -lh "${OUTDIR}"/*.tsv 2>/dev/null | head -20
echo ""
echo "Phase outputs:"
ls -lh "${OUTDIR}"/phase_*.tsv 2>/dev/null
echo ""
echo "Marker package:"
ls -lhR "$ZOODIR" 2>/dev/null
echo ""
echo "Figures:"
ls -lh "${FIGDIR}"/*.png "${FIGDIR}"/*.pdf 2>/dev/null
echo ""
echo "════════════════════════════════════════════════════════"
