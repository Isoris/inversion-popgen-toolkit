#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -A lt200308
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH -t 1:00:00
#SBATCH -J inv_C00_sv_prior
#SBATCH -o logs/inv_C00_sv_prior_%j.out
#SBATCH -e logs/inv_C00_sv_prior_%j.err

# =============================================================================
# LAUNCH_STEP_C00_sv_prior.slurm
#
# Phase 2 / 2c — build the SV-derived evidence prior per chromosome from
# DELLY2 + Manta catalogs. Run ONCE, before STEP_C01a_precompute.
#
# Upstream:    DELLY2 INV/DEL/BND VCFs  (MODULE_4D / MODULE_4B / MODULE_4F)
#              Manta  INV/DEL    VCFs  (MODULE_4H)
#              annotate_population_confidence output (POP_CONF_DIR)
# Downstream:  STEP_C01a_precompute.R  (via --sv_prior_dir)
#
# Usage:
#   sbatch LAUNCH_STEP_C00_sv_prior.slurm                   # all chromosomes
#   sbatch LAUNCH_STEP_C00_sv_prior.slurm C_gar_LG12        # single chrom
#
# Sizing: single-node R job, not parallel. 64G is comfortable headroom for
# a 226-sample DELLY/Manta INV catalog. Walltime ~10–20 min per chrom.
# Bump --cpus-per-task if you want data.table to use more threads.
# =============================================================================

set -euo pipefail
source ~/.bashrc
mamba activate assembly

# ── Load configs ─────────────────────────────────────────────────────────
CONFIG="${CONFIG:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_modules/00_inversion_config.sh}"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
set -a
source "${CONFIG}"
set +a

# Pipeline bridge provides SAMPLE_GROUP / RESULTS_REGISTRY_DIR / SAMPLE_REGISTRY
# (chat-16 additions). Not strictly required by STEP_C00 today — the evidence
# registry write-path doesn't FK-check — but sourcing now makes future-proof.
BRIDGE_SH="${BRIDGE_SH:-${BASE}/utils/pipeline_bridge.sh}"
[[ -f "${BRIDGE_SH}" ]] && source "${BRIDGE_SH}"

inv_init_dirs

mkdir -p logs

# ── Heads-up: SV_PRIOR_DIR mismatch ──────────────────────────────────────
# 00_inversion_config.sh exports SV_PRIOR_DIR=${INVDIR}/03_sv_prior, but
# STEP_C00_build_sv_prior.R writes to
# ${INVDIR}/06_mds_candidates/snake_regions_multiscale/sv_prior (hardcoded
# in the script, ~line 125). Until the script respects SV_PRIOR_DIR, the
# script's path is the source of truth. Resolve it here so the downstream
# STEP_C01a can be pointed at the same location.
SV_PRIOR_OUT="${MDS_DIR}/snake_regions_multiscale/sv_prior"

# ── Optional single-chromosome arg ───────────────────────────────────────
CHR_ARG="${1:-}"
CHR_OPTS=""
if [[ -n "${CHR_ARG}" ]]; then
  CHR_OPTS="--chr ${CHR_ARG}"
fi

echo "================================================================"
echo "  C00: Build SV prior (DELLY2 + Manta INV/DEL/BND)"
echo "  Chrom:       ${CHR_ARG:-all}"
echo "  DELLY INV:   ${DELLY_INV_VCF:-(default)}"
echo "  DELLY DEL:   ${DELLY_DEL_VCF:-(default)}"
echo "  DELLY BND:   ${DELLY_BND_VCF:-(default)}"
echo "  Manta INV:   ${MANTA_INV_VCF:-(default)}"
echo "  Manta DEL:   ${MANTA_DEL_VCF:-(default)}"
echo "  POP_CONF:    ${POP_CONF_DIR:-(default)}"
echo "  Samples:     ${SAMPLES_IND}"
echo "  Out dir:     ${SV_PRIOR_OUT}"
echo "  SAMPLE_GROUP=${SAMPLE_GROUP:-all_226}"
echo "  JobID:       ${SLURM_JOB_ID}"
echo "  Started:     $(date)"
echo "================================================================"

"${RSCRIPT_BIN}" \
  "${BASE}/inversion_modules/phase_2_discovery/2c_precomp/STEP_C00_build_sv_prior.R" \
  ${CHR_OPTS}

echo "================================================================"
echo "  C00 complete — $(date)"
echo "  SV prior RDS files: ${SV_PRIOR_OUT}/sv_prior_<chr>.rds"
echo "  Summary:            ${SV_PRIOR_OUT}/sv_prior_summary.tsv"
echo "  Next: sbatch LAUNCH_STEP_C01a_precompute.slurm"
echo "        (pass --sv_prior_dir ${SV_PRIOR_OUT})"
echo "================================================================"
