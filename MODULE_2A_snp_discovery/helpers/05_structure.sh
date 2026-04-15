#!/usr/bin/env bash
###############################################################################
# helpers/05_structure.sh — Unified structure analysis
#
# Runs PCAngsd, NGSadmix, evalAdmix, best-seed-by-K selection, and canonical
# exports — all in one script. Same script for any sample set.
#
# Usage:
#   bash helpers/05_structure.sh --samples all_samples.txt
#   bash helpers/05_structure.sh --samples pruned_samples.txt
#
# The --samples flag determines:
#   - which sample list is used for NGSadmix/PCAngsd
#   - the sample_set label in output directories and manifests
#   - the tag in canonical best-seed output files
#
# What this script does:
#   1) Determine sample_set label from filename
#   2) Print sbatch commands for PCAngsd (per LG × thin × K)
#   3) Run or print sbatch for global NGSadmix (all K × seeds)
#   4) Run or print sbatch for evalAdmix (all K × seeds)
#   5) Run best-seed-by-K selection (R script)
#   6) Write canonical exports:
#      - all_seed_metrics_by_K.tsv
#      - best_seed_by_K.tsv
#      - best_seed_copied_files.tsv
#      - cluster_palette_by_K.tsv
#      - sample_order_reference.tsv
#      - sample_main_ancestry_by_K.tsv
#
# The practical run order for the full pipeline is:
#   1) bash helpers/05_structure.sh --samples all_samples.txt      (structure_all)
#   2) bash helpers/06_relatedness.sh                               (relatedness)
#   3) bash helpers/05_structure.sh --samples pruned_samples.txt   (structure_pruned)
#   4) bash helpers/merge_structure_summaries.sh                    (merge both)
#
# Called by: run_step1.sh structure_all / structure_pruned
###############################################################################
set -euo pipefail
source "$(dirname "$0")/../config.sh"

# ---- Parse args ----
SAMPLES=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --samples) SAMPLES="$2"; shift 2;;
    *) echo "Usage: $0 --samples <sample_list.txt>"; exit 1;;
  esac
done

[[ -n "$SAMPLES" && -s "$SAMPLES" ]] || { echo "[ERROR] --samples required and must exist" >&2; exit 1; }

# Derive sample_set label from filename (e.g. "all" from all_samples.txt, "pruned" from pruned_samples.txt)
SAMPLE_SET="$(basename "$SAMPLES" | sed 's/_samples\.\(txt\|tsv\)$//; s/\.txt$//; s/\.tsv$//')"
N_SAMP=$(wc -l < "$SAMPLES")

echo "[$(timestamp)] 05_structure: sample_set=${SAMPLE_SET} (${N_SAMP} samples)"
echo "[$(timestamp)] samples_file=${SAMPLES}"

# ---- Output directories scoped by sample_set ----
STRUCT_BASE="${THIN_DIR}/05_structure_${SAMPLE_SET}"
NGSADMIX_DIR="${STRUCT_BASE}/ngsadmix_global/runs_thin${RELATE_THIN}"
EVALADMIX_DIR="${STRUCT_BASE}/ngsadmix_global/evaladmix_thin${RELATE_THIN}"
PCANGSD_DIR="${STRUCT_BASE}/pcangsd_byLG"
EXPORTS_DIR="${STRUCT_BASE}/canonical_exports"

mkdir -p "${NGSADMIX_DIR}" "${EVALADMIX_DIR}" "${PCANGSD_DIR}/logs" "${EXPORTS_DIR}"

# ---- Write .arg ----
ARGFILE="${STRUCT_BASE}/05_structure.arg"
{
  echo -e "key\tvalue"
  echo -e "step\t05_structure"
  echo -e "sample_set\t${SAMPLE_SET}"
  echo -e "samples_file\t${SAMPLES}"
  echo -e "n_samples\t${N_SAMP}"
  echo -e "datetime\t$(timestamp)"
  echo -e "host\t$(hostname)"
  echo -e "relate_thin\t${RELATE_THIN}"
  echo -e "K_range\t${K_MIN}-${K_MAX}"
  echo -e "seeds\t$(IFS=,; echo "${SEEDS[*]}")"
  echo -e "best_seed_rule\t${BEST_SEED_RULE}"
  echo -e "palette_name\t${PALETTE_NAME}"
} > "$ARGFILE"

# ===========================================================================
# 1) PCAngsd (per LG × thin × K)
# ===========================================================================
echo ""
echo "[$(timestamp)] === PCAngsd ==="

for W in "${THIN_FINE[@]}"; do
  LISTF="${PCANGSD_DIR}/beagle_LG_thin_${W}.list"
  BEAGLE_SRC="${THIN_DIR}/04_beagle_byRF_majmin/thin_${W}"
  if [[ -d "$BEAGLE_SRC" ]]; then
    find "$BEAGLE_SRC" -name "*.beagle.gz" | sort -V > "$LISTF" 2>/dev/null || true
    echo "[INFO] thin_${W}: $(wc -l < "$LISTF") beagle files"
  fi
done

echo "[INFO] PCAngsd submit command:"
echo "  sbatch --array=1-<TOTAL>%32 helpers/slurm_pcangsd.sh"
echo "  (using SAMPLE_LIST=${SAMPLES} for --tree-samples)"

# ===========================================================================
# 2) NGSadmix global (all K × seeds on thin-${RELATE_THIN} whole-genome)
# ===========================================================================
echo ""
echo "[$(timestamp)] === NGSadmix global (thin-${RELATE_THIN}) ==="

BEAGLE="${THIN_DIR}/03_merged_beagle/catfish.wholegenome.byRF.thin_${RELATE_THIN}.beagle.gz"
if [[ ! -s "$BEAGLE" ]]; then
  echo "[WARN] Missing global BEAGLE: ${BEAGLE}"
  echo "[WARN] Run make_beagles first"
else
  P="${DEFAULT_THREADS}"
  echo "[INFO] Running NGSadmix K=${K_MIN}..${K_MAX} × ${#SEEDS[@]} seeds..."

  for SEED in "${SEEDS[@]}"; do
    for K in $(seq "$K_MIN" "$K_MAX"); do
      PREFIX="${NGSADMIX_DIR}/thin${RELATE_THIN}_K$(printf "%02d" "$K")_seed${SEED}"
      if [[ -s "${PREFIX}.qopt" && "${FORCE:-0}" -eq 0 ]]; then
        continue
      fi
      echo "[RUN] K=$K seed=$SEED -> $(basename $PREFIX)"
      NGSadmix -likes "$BEAGLE" -K "$K" -minMaf ${NGSADMIX_MINMAF} -P "$P" -seed "$SEED" -o "$PREFIX" \
        > "${PREFIX}.stdout.log" 2>&1 || echo "[WARN] NGSadmix failed: K=$K seed=$SEED"
    done
  done
  echo "[$(timestamp)] NGSadmix done."
fi

# ===========================================================================
# 3) EvalAdmix (all K × seeds)
# ===========================================================================
echo ""
echo "[$(timestamp)] === EvalAdmix ==="

if [[ -s "$BEAGLE" && -d "${EVALADMIX_BIN}" ]]; then
  P="${DEFAULT_THREADS}"
  n_eval=0
  for K in $(seq "$K_MIN" "$K_MAX"); do
    KPAD=$(printf "%02d" "$K")
    for SEED in "${SEEDS[@]}"; do
      PREFIX="thin${RELATE_THIN}_K${KPAD}_seed${SEED}"
      FOPT="${NGSADMIX_DIR}/${PREFIX}.fopt.gz"
      QOPT="${NGSADMIX_DIR}/${PREFIX}.qopt"
      COROUT="${EVALADMIX_DIR}/${PREFIX}.corres.txt"

      if [[ ! -s "$FOPT" || ! -s "$QOPT" ]]; then continue; fi
      if [[ -s "$COROUT" && "${FORCE:-0}" -eq 0 ]]; then continue; fi

      "${EVALADMIX_BIN}/evalAdmix" \
        -beagle "$BEAGLE" -fname "$FOPT" -qname "$QOPT" \
        -P "$P" -nIts ${EVALADMIX_NITS} -misTol ${EVALADMIX_MISTOL} -minMaf ${EVALADMIX_MINMAF} \
        -o "$COROUT" 2>/dev/null || echo "[WARN] evalAdmix failed: ${PREFIX}"
      n_eval=$((n_eval+1))
    done
  done
  echo "[$(timestamp)] evalAdmix: ${n_eval} runs completed."
else
  echo "[WARN] Skipping evalAdmix (missing BEAGLE or evalAdmix binary)"
fi

# ===========================================================================
# 4) Best-seed-by-K selection + canonical exports
# ===========================================================================
echo ""
echo "[$(timestamp)] === Best-seed-by-K selection ==="

BEST_SEED_R="$(dirname "$0")/select_best_seed_by_K.R"
[[ -f "$BEST_SEED_R" ]] || { echo "[ERROR] Missing R script: $BEST_SEED_R" >&2; exit 1; }

Rscript "$BEST_SEED_R" \
  --base-dir "${STRUCT_BASE}/ngsadmix_global" \
  --run-dir "${NGSADMIX_DIR}" \
  --eval-dir "${EVALADMIX_DIR}" \
  --sample-file "${SAMPLES}" \
  --thin-label "${RELATE_THIN}" \
  --k-min "${K_MIN}" \
  --k-max "${K_MAX}" \
  --seeds "$(IFS=,; echo "${SEEDS[*]}")" \
  --palette-name "${PALETTE_NAME}" \
  --sample-set "${SAMPLE_SET}" \
  --out-dir "${EXPORTS_DIR}" \
  || echo "[WARN] Best-seed selection had issues"

# ===========================================================================
# 5) Write .results
# ===========================================================================
RESULTSFILE="${STRUCT_BASE}/05_structure.results"
{
  echo -e "key\tvalue\tdescription"
  echo -e "arg_file\t${ARGFILE}\tParameter record"
  echo -e "sample_set\t${SAMPLE_SET}\tSample set label"
  echo -e "n_samples\t${N_SAMP}\tSamples in this run"
  echo -e "ngsadmix_dir\t${NGSADMIX_DIR}\tNGSadmix run outputs"
  echo -e "evaladmix_dir\t${EVALADMIX_DIR}\tevalAdmix residual matrices"
  echo -e "exports_dir\t${EXPORTS_DIR}\tCanonical best-seed exports"
  for f in all_seed_metrics_by_K.tsv best_seed_by_K.tsv best_seed_copied_files.tsv \
           cluster_palette_by_K.tsv sample_order_reference.tsv sample_main_ancestry_by_K.tsv; do
    p="${EXPORTS_DIR}/${f}"
    [[ -s "$p" ]] && echo -e "export\t${p}\t${f}"
  done
} > "$RESULTSFILE"

echo ""
echo "[$(timestamp)] [DONE] 05_structure sample_set=${SAMPLE_SET}"
echo "[$(timestamp)] Canonical exports: ${EXPORTS_DIR}/"
ls -lh "${EXPORTS_DIR}/"*.tsv 2>/dev/null || true
