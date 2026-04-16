#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=237G
#SBATCH -t 6:00:00
#SBATCH -J ghsl_v53
#SBATCH -o logs/ghsl_v53_%j.out
#SBATCH -e logs/ghsl_v53_%j.err
#SBATCH -A lt200308

# =============================================================================
# LAUNCH_C04_snake3_ghsl_v5.slurm
#
# Snake 3 v5.3: Within-sample haplotype divergence + z-score scoring
#
# Fixes from v5.2:
#   - Bug 1: trailing comma in annot .() block → "Item 16 missing"
#   - Bug 2: rank_mat missing dimnames → karyotype calls have no sample_id
#   - Bug 3: scoring formula inflated by founder structure → 100% PASS
#             Now uses z-score against chromosome-wide baseline
#
# Usage:
#   sbatch LAUNCH_C04_snake3_ghsl_v5.slurm                     # all chroms
#   sbatch LAUNCH_C04_snake3_ghsl_v5.slurm --chrom C_gar_LG01  # single chr
#   sbatch LAUNCH_C04_snake3_ghsl_v5.slurm --test 500          # first 500 windows
# =============================================================================

set -euo pipefail
source ~/.bashrc
mamba activate assembly

CONFIG="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_codebase_v8.5/00_inversion_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }

set -a
source "${CONFIG}"
set +a

inv_init_dirs

PRECOMP_DIR="${MDS_DIR}/snake_regions_multiscale/precomp"
GHSL_PREP_DIR="${MDS_DIR}/snake_regions_multiscale/ghsl_prep"
OUTDIR="${MDS_DIR}/snake3_ghsl_v5"
SCRIPT="${DISCOVERY2DIR}/snakes/STEP_C04_snake3_ghsl_v5.R"

CHROM_OPTS=""
TEST_OPTS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --chrom) CHROM_OPTS="--chrom $2"; shift 2 ;;
    --test)  TEST_OPTS="--test_windows $2"; shift 2 ;;
    *)       shift ;;
  esac
done

mkdir -p "${OUTDIR}"

echo "================================================================"
echo "  Snake 3 v5.3: Within-Sample Haplotype Divergence"
echo "  Precomp:  ${PRECOMP_DIR}"
echo "  GHSL:     ${GHSL_PREP_DIR}"
echo "  Output:   ${OUTDIR}"
echo "  Chrom:    ${CHROM_OPTS:-all}"
echo "  Test:     ${TEST_OPTS:-full}"
echo "  JobID:    ${SLURM_JOB_ID}"
echo "  Started:  $(date)"
echo "================================================================"

"${RSCRIPT_BIN}" "${SCRIPT}" \
  "${PRECOMP_DIR}" "${GHSL_PREP_DIR}" "${OUTDIR}" \
  ${CHROM_OPTS} ${TEST_OPTS} --ncores ${SLURM_CPUS_PER_TASK}

echo ""
echo "================================================================"
echo "  GHSL v5.3 complete — $(date)"
echo "  Output: ${OUTDIR}/"
echo "================================================================"
