#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=64G
#SBATCH -t 0-04:00:00
#SBATCH -J delly_dup_downstream
#SBATCH -o logs/downstream_analysis.%j.out
#SBATCH -e logs/downstream_analysis.%j.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${SCRIPT_DIR}/../00_module4c_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }

set -a
source "${CONFIG}"
set +a

DIR_ANNOT="${OUTDIR}/08_annotation"
DIR_DEPTH="${OUTDIR}/09_depth_support"
DIR_MATDIST="${OUTDIR}/10_mate_distance_qc"
DIR_SUMMARY="${OUTDIR}/11_summary"
DIR_PLOTS="${OUTDIR}/12_plots"

mkdir -p "${DIR_SUMMARY}" "${DIR_PLOTS}"

dv_log "=== DOWNSTREAM DUP ANALYSIS ==="

FINAL_226_VCF="${DIR_FINAL}/catalog_226.DUP.vcf.gz"
GT_MATRIX_226="${DIR_FINAL}/catalog_226.DUP.GT_matrix.tsv"
FUNC_226="${DIR_ANNOT}/catalog_226.functional_class.tsv"
REPEATS_IN="${DIR_ANNOT}/catalog_226.DUPs_in_repeats.bed"
DEPTH="${DIR_DEPTH}/depth_support_226.tsv"
MATE_QC="${DIR_MATDIST}/span_qc_226.tsv"
BINARY_GT="${DIR_SUMMARY}/catalog_226.DUP.binary_gt_for_pca.tsv"

for f in "${FINAL_226_VCF}" "${GT_MATRIX_226}" "${FUNC_226}" "${REPEATS_IN}" "${DEPTH}" "${MATE_QC}" "${SAMPLES_UNRELATED}" "${REF_FAI}"; do
  [[ -f "$f" ]] || { echo "Missing required input: $f" >&2; exit 1; }
done

python3 "${SCRIPT_DIR}/make_dup_pa_tables.py"

QOPT=""
QOPT_SAMPLES=""
if [[ -f "${DELLY_PROJECT}/popstruct_thin/05_ngsadmix_global/best_seed_by_K.tsv" ]]; then
  PA_BEST_SEED_TABLE="${DELLY_PROJECT}/popstruct_thin/05_ngsadmix_global/best_seed_by_K.tsv"
  PA_NGSADMIX_DIR="${DELLY_PROJECT}/popstruct_thin/05_ngsadmix_global/runs_thin500"
  BEST_K=$(awk -F'\t' 'NR>1 { if ($5 > max || NR==2) { max=$5; k=$1 } } END { print k }' "${PA_BEST_SEED_TABLE}")
  BEST_PREFIX=$(awk -F'\t' -v k="${BEST_K}" 'NR>1 && $1==k { print $4 }' "${PA_BEST_SEED_TABLE}")
  CANDIDATE_QOPT="${PA_NGSADMIX_DIR}/${BEST_PREFIX}.qopt"
  if [[ -f "${CANDIDATE_QOPT}" ]]; then
    QOPT="${CANDIDATE_QOPT}"
    QOPT_SAMPLES="${SAMPLES_ALL}"
  fi
fi

ANCESTRY_ARGS=""
if [[ -n "${QOPT}" && -n "${QOPT_SAMPLES}" ]]; then
  ANCESTRY_ARGS="--qopt ${QOPT} --qopt_samples ${QOPT_SAMPLES}"
fi

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

python3 "${SCRIPT_DIR}/08_build_per_sample_summary.py" \
  --master_annot "${DIR_SUMMARY}/master_DUP_annotation_226.tsv" \
  --gt_matrix "${GT_MATRIX_226}" \
  --binary_gt "${BINARY_GT}" \
  --samples_unrelated "${SAMPLES_UNRELATED}" \
  --outdir "${DIR_SUMMARY}" \
  ${ANCESTRY_ARGS}

MARKER_ARGS=""
if [[ -f "${DIR_SUMMARY}/per_sample_DUP_summary.tsv" ]]; then
  MARKER_ARGS="--per_sample_summary ${DIR_SUMMARY}/per_sample_DUP_summary.tsv --gt_matrix ${GT_MATRIX_226}"
fi

python3 "${SCRIPT_DIR}/09_select_markers.py" \
  --master_annot "${DIR_SUMMARY}/master_DUP_annotation_226.tsv" \
  --outdir "${DIR_SUMMARY}" \
  ${MARKER_ARGS}

python3 "${SCRIPT_DIR}/10_gene_summary_tables.py" \
  --master_annot "${DIR_SUMMARY}/master_DUP_annotation_226.tsv" \
  --ref_fai "${REF_FAI}" \
  --outdir "${DIR_SUMMARY}"

dv_log "=== DOWNSTREAM DUP ANALYSIS COMPLETE ==="
