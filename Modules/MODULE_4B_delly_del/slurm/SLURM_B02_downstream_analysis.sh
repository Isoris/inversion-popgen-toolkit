#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=64G
#SBATCH -t 0-04:00:00
#SBATCH -J delly_downstream
#SBATCH -o logs/downstream_analysis.%j.out
#SBATCH -e logs/downstream_analysis.%j.err
# =============================================================================
# 07_downstream_analysis.sh — Build all downstream DEL analysis tables
#
# Runs:
#   07_build_master_annotation.py   → master DEL annotation table (226 + 81)
#   08_build_per_sample_summary.py  → per-sample burden by size/class/repeat + PCA
#   09_select_markers.py            → tier 1/2/3 marker selection
#   10_gene_summary_tables.py       → gene overlap, per-chromosome, top genes
#
# Requires: upstream pipeline complete (steps 01–06)
# =============================================================================
set -euo pipefail
source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${SCRIPT_DIR}/../00_module4b_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
set -a
source "${CONFIG}"
set +a

dv_init_dirs
dv_log "=== DOWNSTREAM DEL ANALYSIS ==="

# ── Resolve paths ──────────────────────────────────────────────────────────
FINAL_226_VCF="${DIR_FINAL}/catalog_226.DEL.vcf.gz"
GT_MATRIX_226="${DIR_FINAL}/catalog_226.DEL.GT_matrix.tsv"
FUNC_226="${DIR_ANNOT}/catalog_226.functional_class.tsv"
REPEATS_IN="${DIR_ANNOT}/catalog_226.DELs_in_repeats.bed"
DEPTH="${DIR_DEPTH}/depth_support_226.tsv"
MATE_QC="${DIR_MATDIST}/mate_distance_qc_226.tsv"
BINARY_GT="${DIR_SUMMARY}/DEL_binary_genotype_matrix.tsv"

# Check all inputs exist
for f in "${FINAL_226_VCF}" "${GT_MATRIX_226}" "${FUNC_226}" "${REPEATS_IN}" \
         "${DEPTH}" "${MATE_QC}" "${BINARY_GT}" "${SAMPLES_UNRELATED}" "${REF_FAI}"; do
  dv_check_file "$f" "$(basename "$f")"
done

# ── Resolve ancestry files ─────────────────────────────────────────────────
# Try to find the best qopt and matching sample list
QOPT=""
QOPT_SAMPLES=""

if [[ -f "${PA_BEST_SEED_TABLE}" ]]; then
  # Auto-detect best K
  BEST_K=$(awk -F'\t' 'NR>1 { if ($5 > max || NR==2) { max=$5; k=$1 } } END { print k }' \
           "${PA_BEST_SEED_TABLE}")
  BEST_K_PAD=$(printf "%02d" "${BEST_K}")
  BEST_PREFIX=$(awk -F'\t' -v k="${BEST_K}" 'NR>1 && $1==k { print $4 }' "${PA_BEST_SEED_TABLE}")

  if [[ -n "${BEST_PREFIX}" ]]; then
    CANDIDATE_QOPT="${PA_NGSADMIX_DIR}/${BEST_PREFIX}.qopt"
    if [[ -f "${CANDIDATE_QOPT}" ]]; then
      QOPT="${CANDIDATE_QOPT}"
      dv_log "  Ancestry: K=${BEST_K}, qopt=${QOPT}"
    fi
  fi

  # Sample list: the global sample list used for NGSadmix
  # This should match row order of qopt — typically samples_all_226.txt
  QOPT_SAMPLES="${SAMPLES_ALL}"
fi

ANCESTRY_ARGS=""
if [[ -n "${QOPT}" && -n "${QOPT_SAMPLES}" ]]; then
  ANCESTRY_ARGS="--qopt ${QOPT} --qopt_samples ${QOPT_SAMPLES}"
  dv_log "  Ancestry args: ${ANCESTRY_ARGS}"
else
  dv_log "  No ancestry data found. Skipping ancestry-aware features."
fi

# ── STEP 7: Master annotation table ───────────────────────────────────────
dv_log "--- 07: Master DEL annotation table ---"
python3 "${SCRIPT_DIR}/07_build_master_annotation.py" \
  --vcf_226 "${FINAL_226_VCF}" \
  --gt_matrix_226 "${GT_MATRIX_226}" \
  --functional_226 "${FUNC_226}" \
  --repeats_in_bed "${REPEATS_IN}" \
  --depth_support "${DEPTH}" \
  --mate_distance_qc "${MATE_QC}" \
  --samples_unrelated "${SAMPLES_UNRELATED}" \
  --outdir "${DIR_SUMMARY}" \
  ${ANCESTRY_ARGS}

MASTER_226="${DIR_SUMMARY}/master_DEL_annotation_226.tsv"
dv_check_file "${MASTER_226}" "Master annotation 226"
dv_log "  Master annotation: $(tail -n +2 "${MASTER_226}" | wc -l) DELs"

# ── STEP 8: Per-sample summary ────────────────────────────────────────────
dv_log "--- 08: Per-sample DEL summary ---"
python3 "${SCRIPT_DIR}/08_build_per_sample_summary.py" \
  --master_annot "${MASTER_226}" \
  --gt_matrix "${GT_MATRIX_226}" \
  --binary_gt "${BINARY_GT}" \
  --samples_unrelated "${SAMPLES_UNRELATED}" \
  --outdir "${DIR_SUMMARY}" \
  ${ANCESTRY_ARGS}

dv_log "  Per-sample summary: $(wc -l < "${DIR_SUMMARY}/per_sample_DEL_summary.tsv") lines"

# ── STEP 9: Marker selection ──────────────────────────────────────────────
dv_log "--- 09: Marker selection ---"
MARKER_ARGS=""
if [[ -f "${DIR_SUMMARY}/per_sample_DEL_summary.tsv" ]]; then
  MARKER_ARGS="--per_sample_summary ${DIR_SUMMARY}/per_sample_DEL_summary.tsv --gt_matrix ${GT_MATRIX_226}"
fi

python3 "${SCRIPT_DIR}/09_select_markers.py" \
  --master_annot "${MASTER_226}" \
  --outdir "${DIR_SUMMARY}" \
  ${MARKER_ARGS}

for tier in tier1_strict tier2_gene tier3_group_informative; do
  f="${DIR_SUMMARY}/markers_${tier}.tsv"
  if [[ -f "$f" ]]; then
    n=$(tail -n +2 "$f" | wc -l)
    dv_log "  Tier ${tier}: ${n} markers"
  fi
done

# ── STEP 10: Gene and chromosome summaries ─────────────────────────────────
dv_log "--- 10: Gene and chromosome summaries ---"
python3 "${SCRIPT_DIR}/10_gene_summary_tables.py" \
  --master_annot "${MASTER_226}" \
  --ref_fai "${REF_FAI}" \
  --outdir "${DIR_SUMMARY}"

for f in gene_DEL_overlap_table gene_summary_collapsed top_genes_report per_chromosome_summary; do
  ff="${DIR_SUMMARY}/${f}.tsv"
  if [[ -f "$ff" ]]; then
    dv_log "  ${f}: $(tail -n +2 "$ff" | wc -l) rows"
  fi
done

# ── Final listing ──────────────────────────────────────────────────────────
dv_log ""
dv_log "=== DOWNSTREAM ANALYSIS COMPLETE ==="
dv_log "New files in ${DIR_SUMMARY}:"
ls -lh "${DIR_SUMMARY}"/master_DEL_annotation_*.tsv \
       "${DIR_SUMMARY}"/per_sample_DEL_summary.tsv \
       "${DIR_SUMMARY}"/markers_*.tsv \
       "${DIR_SUMMARY}"/gene_*.tsv \
       "${DIR_SUMMARY}"/top_genes_report.tsv \
       "${DIR_SUMMARY}"/per_chromosome_summary.tsv \
       "${DIR_SUMMARY}"/frequency_class_summary.tsv \
       "${DIR_SUMMARY}"/size_class_summary.tsv \
       "${DIR_SUMMARY}"/marker_selection_summary.tsv \
       2>/dev/null | while read line; do dv_log "  ${line}"; done
