#!/usr/bin/env bash
# =============================================================================
# run_breakpoint_validation.sh — v2: dual-caller (DELLY2 + Manta)
# =============================================================================
set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${SCRIPT_DIR}/00_breakpoint_validation_config.sh"

bpv_init_dirs

echo "=== MODULE_5A2: Breakpoint Validation Pipeline (v2: dual-caller) ==="
echo "Script directory: ${SCRIPT_DIR}"
echo "DELLY INV VCF:    ${DELLY_INV_VCF}"
echo "Manta INV VCF:    ${MANTA_INV_VCF}"
echo "Markdup BAMs:     ${DELLY_MARKDUP_DIR}"
echo "Output root:      ${BPV_ROOT}"
echo ""

# STEP 1: Extract from BOTH callers + merge + match Snake
JID1=$(sbatch --parsable \
  --account=lt200308 \
  -o "${BPV_LOGS}/01_extract.%j.out" -e "${BPV_LOGS}/01_extract.%j.err" \
  "${SCRIPT_DIR}/01_extract_inv_candidates.sh")
echo "  [1/6] Extract candidates (DELLY+Manta):  Job ${JID1}"

# STEP 2: Per-sample BAM evidence (caller-agnostic pysam)
JID2=$(sbatch --parsable --dependency=afterok:${JID1} \
  --account=lt200308 \
  -p compute -N 1 -n 32 --mem=128G -t 1-00:00:00 \
  -J bpv_s02 \
  -o "${BPV_LOGS}/02_evidence.%j.out" -e "${BPV_LOGS}/02_evidence.%j.err" \
  --wrap="
    set -euo pipefail
    source ~/.bashrc
    mamba activate assembly
    source '${SCRIPT_DIR}/00_breakpoint_validation_config.sh'

    SNAKE2_ARG=''
    [[ -f '${SNAKE2_BANDS}' ]] && SNAKE2_ARG='--snake2_bands ${SNAKE2_BANDS}'
    DOSAGE_ARG=''
    [[ -d '${DOSAGE_DIR:-/dev/null}' ]] && DOSAGE_ARG='--dosage_dir ${DOSAGE_DIR}'

    python3 '${SCRIPT_DIR}/02_extract_breakpoint_evidence.py' \
      --candidates '${BPV_EVIDENCE}/matched_delly_snake_candidates.tsv' \
      --bam_dir '${DELLY_MARKDUP_DIR}' \
      --samples '${SAMPLES_ALL}' \
      \${SNAKE2_ARG} \${DOSAGE_ARG} \
      --outdir '${BPV_EVIDENCE}' \
      --bp_window ${BP_WINDOW} --min_mapq ${MIN_MAPQ}
  ")
echo "  [2/6] BAM evidence:    Job ${JID2} (after ${JID1})"

# STEP 3: Statistical tests + seed generation
JID3=$(sbatch --parsable --dependency=afterok:${JID2} \
  --account=lt200308 \
  -p compute -N 1 -n 8 --mem=32G -t 0-02:00:00 \
  -J bpv_s03 \
  -o "${BPV_LOGS}/03_tests.%j.out" -e "${BPV_LOGS}/03_tests.%j.err" \
  --wrap="
    set -euo pipefail
    source ~/.bashrc
    mamba activate assembly

    python3 '${SCRIPT_DIR}/03_statistical_tests_and_seeds.py' \
      --evidence_dir '${BPV_EVIDENCE}' \
      --stats_dir '${BPV_STATS}' \
      --seeds_dir '${BPV_SEEDS}' \
      --seed_min_concordance ${SEED_MIN_CONCORDANCE} \
      --seed_min_support_frac ${SEED_MIN_SUPPORT_FRAC} \
      --seed_min_samples ${SEED_MIN_SAMPLES}
  ")
echo "  [3/6] Statistical tests:   Job ${JID3} (after ${JID2})"

# STEP 4: Plots
JID4=$(sbatch --parsable --dependency=afterok:${JID3} \
  --account=lt200308 \
  -p compute -N 1 -n 8 --mem=32G -t 0-02:00:00 \
  -J bpv_s04 \
  -o "${BPV_LOGS}/04_plots.%j.out" -e "${BPV_LOGS}/04_plots.%j.err" \
  --wrap="
    set -euo pipefail
    source ~/.bashrc
    mamba activate assembly

    python3 '${SCRIPT_DIR}/04_validation_plots.py' \
      --evidence_dir '${BPV_EVIDENCE}' \
      --stats_dir '${BPV_STATS}' \
      --plot_dir '${BPV_PLOTS}'
  ")
echo "  [4/6] Plots:               Job ${JID4} (after ${JID3})"

# STEP 5: Cross-caller concordance report (reads unified table from STEP01)
JID5=$(sbatch --parsable --dependency=afterok:${JID3} \
  --account=lt200308 \
  -p compute -N 1 -n 4 --mem=16G -t 0-01:00:00 \
  -J bpv_s05 \
  -o "${BPV_LOGS}/05_concordance.%j.out" -e "${BPV_LOGS}/05_concordance.%j.err" \
  --wrap="
    set -euo pipefail
    source ~/.bashrc
    mamba activate assembly
    source '${SCRIPT_DIR}/00_breakpoint_validation_config.sh'

    python3 '${SCRIPT_DIR}/05_delly_manta_concordance.py' \
      --candidates '${BPV_EVIDENCE}/matched_delly_snake_candidates.tsv' \
      --tests '${BPV_STATS}/all_candidates_tests.tsv' \
      --outdir '${BPV_CONCORDANCE}'
  ")
echo "  [5/6] Concordance report:  Job ${JID5} (after ${JID3})"

# STEP 6: BND inversion signal (DELLY CT + Manta RAW bracket patterns)
BND_ARGS=""
[[ -f "${DELLY_BND_VCF}" ]] && BND_ARGS="${BND_ARGS} --delly_bnd_vcf ${DELLY_BND_VCF}"
[[ -f "${MANTA_RAW_MERGED_VCF}" ]] && BND_ARGS="${BND_ARGS} --manta_raw_vcf ${MANTA_RAW_MERGED_VCF}"
[[ -f "${DELLY_INV_VCF}" ]] && BND_ARGS="${BND_ARGS} --delly_inv_vcf ${DELLY_INV_VCF}"
[[ -f "${MANTA_INV_VCF}" ]] && BND_ARGS="${BND_ARGS} --manta_inv_vcf ${MANTA_INV_VCF}"

if [[ -n "${BND_ARGS}" ]]; then
  JID6=$(sbatch --parsable --dependency=afterok:${JID3} \
    --account=lt200308 \
    -p compute -N 1 -n 4 --mem=16G -t 0-01:00:00 \
    -J bpv_s06 \
    -o "${BPV_LOGS}/06_bnd_signal.%j.out" -e "${BPV_LOGS}/06_bnd_signal.%j.err" \
    --wrap="
      set -euo pipefail
      source ~/.bashrc
      mamba activate assembly

      python3 '${SCRIPT_DIR}/06_bnd_inversion_signal.py' \
        ${BND_ARGS} \
        --outdir '${BPV_ROOT}/06_bnd_signal'
    ")
  echo "  [6/6] BND inversion signal:  Job ${JID6} (after ${JID3})"
else
  echo "  [6/6] Skipped: No BND catalogs found"
fi

echo ""
echo "=== Pipeline submitted ==="
echo "Chain: ${JID1} → ${JID2} → ${JID3} → ${JID4}"
echo "                                   ├→ ${JID5:-skip} (concordance)"
echo "                                   └→ ${JID6:-skip} (BND signal)"
echo ""
echo "Monitor:  squeue -u \$(whoami)"
echo "Logs:     tail -f ${BPV_LOGS}/*.out"
