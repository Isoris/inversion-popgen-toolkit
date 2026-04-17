#!/usr/bin/env bash
# =============================================================================
# run_phase4.sh — Orchestrator for Phase 4 (catalog + groups + validation +
#                 group-dependent cheats + final re-score)
# =============================================================================
# Usage: bash run_phase4.sh [--dry-run] [--chroms LG01,LG02,...] [--resume]
#
# Submits SLURM jobs with afterok dependencies to enforce the v10 phasing:
#
#   Phase 4a — parallel, no ordering:
#       4a_scoring_pass1  (C01d without --boundary_dir / --hyp_dir)
#       4a_boundary       (C01g)
#
#   Phase 4b — depends on 4a_scoring_pass1 (needs the candidate list):
#       4b_decomposition  (C01i — registers 4 groups per candidate)
#
#   Phase 4c — depends on 4b_decomposition:
#       4c_hypothesis     (C01f — validates groups via T8/T9/T10 + Layer D)
#
#   Phase 4d — depends on 4c_hypothesis:
#       4d_group_cheats   (cheat6, cheat30, burden — gated by validation level)
#
#   Phase 4e — depends on 4a_boundary + 4c_hypothesis + 4d_group_cheats:
#       4e_final_score    (C01d pass-2 + characterize + classify)
#
# Requires:
#   - BASE exported or set below
#   - 00_inversion_config.sh in the module dir, sourced for paths
#   - LAUNCH_*.sh wrappers per module (existing convention)
# =============================================================================

set -euo pipefail

# ── Config ───────────────────────────────────────────────────────────────────
BASE="${BASE:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"
MODULE_DIR="${BASE}/inversion_modules"
CONFIG="${MODULE_DIR}/00_inversion_config.sh"
REGISTRY_DIR="${REGISTRY_DIR:-${BASE}/sample_registry}"

if [[ -f "${CONFIG}" ]]; then
  # shellcheck disable=SC1090
  source "${CONFIG}"
fi

SLURM_ACCT="${SLURM_ACCT:-lt200308}"
LOGDIR="${LOGDIR:-${MODULE_DIR}/logs/phase4_$(date +%Y%m%d_%H%M%S)}"
mkdir -p "${LOGDIR}"

# ── Args ─────────────────────────────────────────────────────────────────────
DRY_RUN=0
RESUME=0
CHROMS_CSV=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run) DRY_RUN=1; shift ;;
    --resume)  RESUME=1; shift ;;
    --chroms)  CHROMS_CSV="$2"; shift 2 ;;
    *) echo "[run_phase4] unknown arg: $1" >&2; exit 2 ;;
  esac
done

if [[ -z "${CHROMS_CSV}" ]]; then
  # Default: all 28 LG chromosomes
  CHROMS_CSV=$(seq -f "LG%02g" 1 28 | paste -sd, -)
fi
IFS=',' read -r -a CHROMS <<< "${CHROMS_CSV}"
N_CHR="${#CHROMS[@]}"

echo "[run_phase4] base=${BASE}"
echo "[run_phase4] registry=${REGISTRY_DIR}"
echo "[run_phase4] chroms=${CHROMS_CSV} (n=${N_CHR})"
echo "[run_phase4] logs=${LOGDIR}"
echo "[run_phase4] dry_run=${DRY_RUN}"

# ── Helper: sbatch with logging, returns job id ──────────────────────────────
submit_job () {
  local name="$1"; shift
  local dep="$1"; shift   # may be empty string
  local script="$1"; shift
  local extra_args=("$@")

  local dep_flag=""
  if [[ -n "${dep}" ]]; then
    dep_flag="--dependency=afterok:${dep}"
  fi

  if [[ "${DRY_RUN}" -eq 1 ]]; then
    echo "  [dry-run] sbatch -A ${SLURM_ACCT} -J ${name} ${dep_flag} ${script} ${extra_args[*]}"
    echo "99999999"  # fake job id
    return
  fi

  local jid
  jid=$(sbatch --parsable \
    -A "${SLURM_ACCT}" \
    -J "${name}" \
    -o "${LOGDIR}/${name}_%A_%a.out" \
    -e "${LOGDIR}/${name}_%A_%a.err" \
    ${dep_flag} \
    "${script}" "${extra_args[@]}")
  echo "[run_phase4] submitted ${name}: jobid=${jid}" >&2
  echo "${jid}"
}

# =============================================================================
# Phase 4a — parallel
# =============================================================================
# 4a_scoring_pass1 — one job per chromosome (array), C01d pass-1
# The launcher reads CHROM from $SLURM_ARRAY_TASK_ID
LAUNCH_C01D_PASS1="${MODULE_DIR}/phase_4_postprocessing/4a_existence_layers/LAUNCH_C01d_scoring_pass1.sh"
LAUNCH_C01G="${MODULE_DIR}/phase_4_postprocessing/4a_existence_layers/LAUNCH_C01g_boundary.sh"
# FIX 49 (chat 10): LAUNCH_C01I removed — replaced by run_phase4b.sh
# sub-orchestrator. 4b is now a 4-script DAG, not a single launcher.
LAUNCH_C01F="${MODULE_DIR}/phase_4_postprocessing/4c_group_validation/LAUNCH_C01f_hypothesis.sh"
LAUNCH_C01D_PASS2="${MODULE_DIR}/phase_4_postprocessing/4e_final_classification/LAUNCH_C01d_scoring_pass2.sh"
LAUNCH_GROUP_CHEATS="${MODULE_DIR}/phase_4_postprocessing/orchestrator/LAUNCH_group_cheats.sh"
LAUNCH_CHAR_CLASSIFY="${MODULE_DIR}/phase_4_postprocessing/4e_final_classification/LAUNCH_characterize_classify.sh"

# Use array notation: the launchers should accept $SLURM_ARRAY_TASK_ID
# and map it to a chromosome. Adjust below if your launchers take --chrom.

JID_SCORE1=$(submit_job "4a_score1"  ""  "${LAUNCH_C01D_PASS1}" --array="1-${N_CHR}" --chroms "${CHROMS_CSV}")
JID_BOUND=$(submit_job  "4a_bound"   ""  "${LAUNCH_C01G}"       --array="1-${N_CHR}" --chroms "${CHROMS_CSV}")

# =============================================================================
# Phase 4b — depends on 4a_scoring_pass1 (candidate list must exist)
# =============================================================================
# FIX 49 (chat 10): 4b is now a 4-script sub-DAG, not a single launcher.
# Invoke run_phase4b.sh which submits the four jobs with their internal
# dependency graph (decompose -> multi_recomb, nested_comp parallel, seal
# gates on all three). Grep the final seal jobid from the sub-orchestrator's
# stdout (emitted as "PHASE4B_SEAL_JID=<id>" on the last line — see
# run_phase4b.sh FIX 48).
#
# Earlier versions of this script called LAUNCH_C01i_decomposition.sh as
# a single job. That launcher is obsolete post-4b-rewrite. The sub-DAG
# approach preserves the seal-gate contract: 4c still depends on 4b
# completion, now expressed as dependency on the seal job specifically.
#
# NOTE: run_phase4b.sh itself submits jobs; it does not BE a SLURM job.
# We run it on the login node (or wherever run_phase4.sh is invoked).
# It respects ${DRY_RUN} via the --dry-run flag.
DRY_FLAG=""
[[ "${DRY_RUN}" -eq 1 ]] && DRY_FLAG="--dry-run"

echo "[run_phase4] invoking run_phase4b.sh sub-orchestrator..."
PHASE4B_OUT=$(bash "${MODULE_DIR}/orchestrator/run_phase4b.sh" \
  ${DRY_FLAG} \
  --chroms "${CHROMS_CSV}" \
  2>&1 | tee /dev/stderr)
JID_DECOMP=$(echo "${PHASE4B_OUT}" | grep '^PHASE4B_SEAL_JID=' | tail -n1 | cut -d= -f2)
if [[ -z "${JID_DECOMP}" ]]; then
  echo "[run_phase4] ERROR: could not parse PHASE4B_SEAL_JID from run_phase4b.sh output" >&2
  exit 3
fi
echo "[run_phase4] 4b seal job id = ${JID_DECOMP}"

# =============================================================================
# Phase 4c — depends on 4b_decomposition (groups must be registered)
# =============================================================================
JID_HYP=$(submit_job    "4c_hyp"     "${JID_DECOMP}"  "${LAUNCH_C01F}"  --chroms "${CHROMS_CSV}")

# =============================================================================
# Phase 4d — depends on 4c_hypothesis (validation levels must be set)
# =============================================================================
JID_CHEATS=$(submit_job "4d_cheats"  "${JID_HYP}"     "${LAUNCH_GROUP_CHEATS}" --chroms "${CHROMS_CSV}")

# =============================================================================
# Phase 4e — depends on 4a_boundary + 4c_hypothesis + 4d_group_cheats
# =============================================================================
# Multiple dependencies: afterok:jid1:jid2:jid3
DEP_4E="${JID_BOUND}:${JID_HYP}:${JID_CHEATS}"
JID_FINAL=$(submit_job  "4e_final"   "${DEP_4E}"      "${LAUNCH_C01D_PASS2}" --chroms "${CHROMS_CSV}")
JID_CHAR=$(submit_job   "4e_char"    "${JID_FINAL}"   "${LAUNCH_CHAR_CLASSIFY}" --chroms "${CHROMS_CSV}")

echo ""
echo "[run_phase4] submitted jobs:"
echo "  4a_score1  = ${JID_SCORE1}"
echo "  4a_bound   = ${JID_BOUND}"
echo "  4b_decomp  = ${JID_DECOMP}"
echo "  4c_hyp     = ${JID_HYP}"
echo "  4d_cheats  = ${JID_CHEATS}"
echo "  4e_final   = ${JID_FINAL}"
echo "  4e_char    = ${JID_CHAR}"
echo ""
echo "[run_phase4] monitor with:"
echo "  squeue -u \$USER -j ${JID_SCORE1},${JID_BOUND},${JID_DECOMP},${JID_HYP},${JID_CHEATS},${JID_FINAL},${JID_CHAR}"
echo ""
echo "[run_phase4] logs at ${LOGDIR}"
