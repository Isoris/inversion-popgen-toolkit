#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -A lt200308
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -t 3:00:00
#SBATCH -J inv_4c_C01f
#SBATCH -o logs/inv_4c_C01f_%j.out
#SBATCH -e logs/inv_4c_C01f_%j.err

# =============================================================================
# LAUNCH_C01f_hypothesis.sh
#
# Phase 4c — hypothesis tests. Validates the groups sealed by 4b.4 using
# T8/T9/T10 + Layer D evidence. Emits verdicts (group_validation_before /
# _after, quality_flags, family_linkage columns) + per-candidate plots.
#
# Called by: orchestrator/run_phase4.sh (depends on 4b.4 seal).
#
# Upstream:    4a C01d pass-1 (candidate_scores.tsv.gz)
#              4b.4 seal (sealed groups in sample_registry)
#              2c triangles (triangle dir)
#              2c precomp
#              Module 2A relatedness pairs / samples / pruned list
# Downstream:  4d group cheats, 4e final classification
# =============================================================================

set -euo pipefail
source ~/.bashrc
mamba activate assembly

CONFIG="${CONFIG:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_modules/00_inversion_config.sh}"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
set -a
source "${CONFIG}"
set +a

BRIDGE_SH="${BRIDGE_SH:-${BASE}/utils/pipeline_bridge.sh}"
[[ -f "${BRIDGE_SH}" ]] || { echo "Missing bridge" >&2; exit 1; }
source "${BRIDGE_SH}"

inv_init_dirs
mkdir -p logs

# ── Parse orchestrator-style args (--chroms is informational) ────────────
CHROMS_CSV=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --chroms) CHROMS_CSV="$2"; shift 2 ;;
    --array=*) shift ;;
    --array)   shift 2 ;;
    *) echo "[C01f] unknown arg (ignored): $1" >&2; shift ;;
  esac
done

# ── Args ─────────────────────────────────────────────────────────────────
SCORES="${SCORES:-${SCORING_DIR}/pass1/candidate_scores.tsv.gz}"
TRIANGLES="${TRIANGLES:-${TRIANGLE_DIR}}"
PRECOMP_DIR_RUN="${PRECOMP_DIR_RUN:-${MDS_DIR}/snake_regions_multiscale/precomp}"
RELATEDNESS="${RELATEDNESS:-${BASE}/popstruct_thin/06_relatedness/ngsRelate.pairs.tsv}"
SAMPLES="${SAMPLES:-${SAMPLE_LIST}}"
PRUNED="${PRUNED:-${PRUNED_LIST}}"
OUTDIR_RUN="${OUTDIR_RUN:-${MDS_DIR}/snake_regions_multiscale/hypothesis_tests}"
TIER_MAX="${TIER_MAX:-2}"

mkdir -p "${OUTDIR_RUN}"

for f in "${SCORES}" "${SAMPLES}"; do
  [[ -f "$f" ]] || { echo "[C01f] ERROR: missing file — $f" >&2; exit 2; }
done
for d in "${TRIANGLES}" "${PRECOMP_DIR_RUN}"; do
  [[ -d "$d" ]] || { echo "[C01f] ERROR: missing dir — $d" >&2; exit 2; }
done

RELATE_OPTS=""
if [[ -f "${RELATEDNESS}" ]]; then
  RELATE_OPTS="--relatedness ${RELATEDNESS}"
else
  echo "[C01f] WARN: relatedness file not found — ${RELATEDNESS}" >&2
fi

PRUNED_OPTS=""
if [[ -f "${PRUNED}" ]]; then
  PRUNED_OPTS="--pruned_samples ${PRUNED}"
fi

echo "================================================================"
echo "  4c / C01f — hypothesis tests"
echo "  Scores:       ${SCORES}"
echo "  Triangles:    ${TRIANGLES}"
echo "  Precomp:      ${PRECOMP_DIR_RUN}"
echo "  Relatedness:  ${RELATEDNESS:-(skipped)}"
echo "  Samples:      ${SAMPLES}"
echo "  Pruned:       ${PRUNED:-(skipped)}"
echo "  Out dir:      ${OUTDIR_RUN}"
echo "  Tier max:     ${TIER_MAX}"
echo "  Chroms:       ${CHROMS_CSV:-all}"
echo "  SAMPLE_GROUP: ${SAMPLE_GROUP:-all_226}"
echo "  JobID:        ${SLURM_JOB_ID}"
echo "  Started:      $(date)"
echo "================================================================"

"${RSCRIPT_BIN}" \
  "${BASE}/inversion_modules/phase_4_postprocessing/4c_group_validation/STEP_C01f_hypothesis_tests.R" \
  --scores    "${SCORES}" \
  --triangles "${TRIANGLES}" \
  --precomp   "${PRECOMP_DIR_RUN}" \
  --samples   "${SAMPLES}" \
  --outdir    "${OUTDIR_RUN}" \
  --tier_max  "${TIER_MAX}" \
  ${RELATE_OPTS} ${PRUNED_OPTS}

echo "================================================================"
echo "  4c complete — $(date)"
echo "  Verdicts:       ${OUTDIR_RUN}/hypothesis_verdicts.tsv"
echo "  Decision tree:  ${OUTDIR_RUN}/hypothesis_decision_tree.tsv"
echo "  Test results:   ${OUTDIR_RUN}/hypothesis_test_results.tsv.gz"
echo "  Plots:          ${OUTDIR_RUN}/plots/"
echo "  Next (4d):      orchestrator/LAUNCH_group_cheats.sh"
echo "================================================================"
