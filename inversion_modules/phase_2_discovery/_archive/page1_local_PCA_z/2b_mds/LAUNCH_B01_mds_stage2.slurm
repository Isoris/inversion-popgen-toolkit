#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32GB
#SBATCH -t 01:00:00
#SBATCH -J inv_10v2_s2
#SBATCH -o logs/inv_10v2_s2.%j.out
#SBATCH -e logs/inv_10v2_s2.%j.err

# =============================================================================
# LAUNCH_B01_mds_stage2.slurm  (v8.5)
#
# Stage 2: merge per-chr MDS results into final .mds.rds + candidates
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

MDS_MODE_ARG="${1:-chunked_4x}"
A2_STEPS="${DISCOVERY2DIR}/steps"
mkdir -p "${MDS_DIR}" "${INVDIR}/logs"

echo "================================================================"
echo "  STEP10v2 Stage 2 — MDS merge (mode=${MDS_MODE_ARG})"
echo "  JobID: ${SLURM_JOB_ID}"
echo "  Started: $(date)"
echo "================================================================"

"${RSCRIPT_BIN}" "${A2_STEPS}/STEP_B01_mds_multimode_stage2.R" \
  --tmpdir "${MDS_DIR}/tmp" \
  --outprefix "${MDS_DIR}/inversion_localpca" \
  --mds_mode "${MDS_MODE_ARG}" \
  --gap_bp 500000 \
  --min_windows 3

echo ""
echo "[DONE] STEP10v2 Stage 2 — MDS merge complete — $(date)"
echo "  .mds.rds:    ${MDS_DIR}/inversion_localpca.mds.rds"
echo "  Candidates:  ${MDS_DIR}/inversion_localpca.candidate_regions.tsv.gz"
