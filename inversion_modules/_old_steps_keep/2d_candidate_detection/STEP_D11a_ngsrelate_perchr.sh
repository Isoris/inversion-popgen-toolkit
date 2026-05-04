#!/bin/bash
#SBATCH --job-name=peel_ngsrel
#SBATCH --account=lt200308-agbsci
#SBATCH --partition=compute
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-28
#SBATCH --output=logs/peel_ngsrel_%a.out
#SBATCH --error=logs/peel_ngsrel_%a.err

# =============================================================================
# 11_run_ngsrelate_perchr.sh — Per-chromosome ngsRelate + NaToRA pruning
# =============================================================================
#
# Runs ngsRelate on EACH chromosome separately using per-chr BEAGLE files.
# Then applies NaToRA-style greedy pruning per chromosome.
# Output: per_chr_pruned.tsv (consumed by STEP_C01n --pruned_list)
#
# WHY PER-CHROMOSOME:
#   Genome-wide ngsRelate averages kinship over 28 chromosomes.
#   Two fish can be unrelated genome-wide but share a whole LG by IBD.
#   Per-chromosome ngsRelate catches this.
#
# INVERSION CONFOUND (honest warning):
#   ngsRelate on a chromosome WITH a big inversion will show inflated
#   kinship between inversion carriers. This is NOT a bug for the peeling
#   diagnostic — it's actually what we want. If ngsRelate says "fish A and
#   fish B are related on LG01" and we peel them and the block disappears,
#   that tells us those fish were essential to the signal. Whether it's
#   family or inversion is what L2/L3 help disambiguate.
#
#   For a "clean" relatedness estimate you'd need to mask inversions first.
#   But that's circular (need to know where inversions are to clean them).
#   The pragmatic approach: run ngsRelate as-is, label the output honestly,
#   use it as a diagnostic not a filter.
#
# REQUIREMENTS:
#   - ngsRelate binary (conda: ngsRelate or compile from source)
#   - Per-chr BEAGLE GL files in BEAGLE_DIR/<chr>.beagle.gz
#   - Sample list
#
# Each chromosome: ~2-5 min with 226 samples and thinned SNPs.
# All 28 in parallel: ~5 min wall time.
#
# =============================================================================
set -euo pipefail

BASE="${BASE:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"
BEAGLE_DIR="${BASE}/popstruct_thin/04_beagle_byRF_majmin"
SAMPLE_LIST="${BASE}/popstruct_thin/list_of_samples_one_per_line_same_bamfile_list.tsv"
OUTDIR="${BASE}/inversion_codebase_v8.5/MODULE_5A2_Discovery_Core/inv_detect_out_v9.3/ngsrelate_perchr"
NGSRELATE_BIN="${NGSRELATE_BIN:-ngsRelate}"
KIN_THRESHOLD="${KIN_THRESHOLD:-0.05}"

source activate assembly
mkdir -p "${OUTDIR}" logs

CHR=$(printf "C_gar_LG%02d" ${SLURM_ARRAY_TASK_ID})
echo "=== ngsRelate: ${CHR} === $(date)"

# ---- Find BEAGLE file ----
BEAGLE_FILE=""
for pattern in "${BEAGLE_DIR}/${CHR}.beagle.gz" \
               "${BEAGLE_DIR}/${CHR}_thin*.beagle.gz" \
               "${BEAGLE_DIR}/*${CHR}*.beagle.gz"; do
  f=$(ls ${pattern} 2>/dev/null | head -1)
  if [[ -n "$f" && -f "$f" ]]; then
    BEAGLE_FILE="$f"
    break
  fi
done

if [[ -z "${BEAGLE_FILE}" ]]; then
  echo "No BEAGLE file found for ${CHR}. Skipping."
  # Write empty result
  echo -e "a\tb\trab\tFratij\tR0\tR1\tKING" > "${OUTDIR}/${CHR}.ngsrelate.res"
  exit 0
fi

echo "  BEAGLE: ${BEAGLE_FILE}"

# ---- Count samples ----
N_SAMPLES=$(zcat "${BEAGLE_FILE}" | head -1 | tr '\t' '\n' | grep -c "Ind" | head -1 || echo 226)
# Alternative: from sample list
if [[ -f "${SAMPLE_LIST}" ]]; then
  N_SAMPLES=$(wc -l < "${SAMPLE_LIST}")
fi
echo "  Samples: ${N_SAMPLES}"

# ---- Run ngsRelate ----
echo "  Running ngsRelate..."
${NGSRELATE_BIN} \
  -G "${BEAGLE_FILE}" \
  -n "${N_SAMPLES}" \
  -O "${OUTDIR}/${CHR}.ngsrelate.res" \
  -p 4 \
  2> "${OUTDIR}/${CHR}.ngsrelate.log" || {
    echo "  ngsRelate failed for ${CHR}"
    exit 1
  }

N_PAIRS=$(wc -l < "${OUTDIR}/${CHR}.ngsrelate.res")
echo "  Pairs: ${N_PAIRS}"
echo "Done: ${CHR} $(date)"
