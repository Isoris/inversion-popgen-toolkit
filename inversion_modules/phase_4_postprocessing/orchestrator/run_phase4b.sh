#!/usr/bin/env bash
# =============================================================================
# run_phase4b.sh — Phase 4b rewrite orchestrator
# =============================================================================
# Submits 4 SLURM jobs for phase 4b:
#
#   4b.1  decompose     (parallel with 4b.2 and 4b.3 in the ideal case,
#                       but 4b.2 reads per_window_class.rds from 4b.1, so
#                       4b.2 depends on 4b.1)
#   4b.2  multi_recomb  (depends on 4b.1)
#   4b.3  nested_comp   (parallel with 4b.1 — reads Engine B only, no
#                       dependency on decompose)
#   4b.4  seal          (depends on 4b.1, 4b.2, 4b.3 all completing)
#
# So the actual DAG is:
#
#                    ┌─→ 4b.1 ──→ 4b.2 ──┐
#                    │                   ├─→ 4b.4
#                    └─→ 4b.3 ───────────┘
#
# 4b.3 runs in parallel with 4b.1+4b.2; 4b.4 gates on all three.
#
# Usage: bash run_phase4b.sh [--dry-run] [--chroms LG01,LG02,...]
#                            [--candidates path] [--tier_max N]
# =============================================================================

set -euo pipefail

BASE="${BASE:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"
MODULE_DIR="${MODULE_DIR:-${BASE}/inversion_modules}"
SLURM_ACCT="${SLURM_ACCT:-lt200308}"
LOGDIR="${LOGDIR:-${MODULE_DIR}/logs/phase4b_$(date +%Y%m%d_%H%M%S)}"
mkdir -p "${LOGDIR}"

# ── Defaults ─────────────────────────────────────────────────────────────────
DRY_RUN=0
CANDIDATES="${MODULE_DIR}/catalog/candidate_scores.tsv.gz"
TIER_MAX=3
CHROMS_CSV=""

DECOMP_OUT="${MODULE_DIR}/phase4b_out/decompose"
RECOMB_OUT="${MODULE_DIR}/phase4b_out/multi_recomb"
NESTED_OUT="${MODULE_DIR}/phase4b_out/nested_comp"
SEAL_OUT="${MODULE_DIR}/phase4b_out/seal"

Q_CACHE_DIR="${Q_CACHE_DIR:-${BASE}/unified_ancestry/local_Q}"

# ── Args ─────────────────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run)    DRY_RUN=1; shift ;;
    --chroms)     CHROMS_CSV="$2"; shift 2 ;;
    --candidates) CANDIDATES="$2"; shift 2 ;;
    --tier_max)   TIER_MAX="$2"; shift 2 ;;
    *) echo "[run_phase4b] unknown arg: $1" >&2; exit 2 ;;
  esac
done

echo "[run_phase4b] candidates=${CANDIDATES}"
echo "[run_phase4b] tier_max=${TIER_MAX}"
echo "[run_phase4b] logs=${LOGDIR}"
echo "[run_phase4b] q_cache=${Q_CACHE_DIR}"
[[ -d "${Q_CACHE_DIR}" ]] || echo "[run_phase4b] NOTE: Q cache absent; nested_comp will write unknown_no_engine_b stubs"

# ── sbatch wrapper ───────────────────────────────────────────────────────────
submit_job() {
  local name="$1"; shift
  local dep="$1"; shift
  local time="$1"; shift
  local cpus="$1"; shift
  local mem="$1"; shift
  local cmd="$*"

  local dep_flag=""
  if [[ -n "${dep}" ]]; then
    dep_flag="--dependency=afterok:${dep}"
  fi

  if [[ "${DRY_RUN}" -eq 1 ]]; then
    echo "  [dry-run] sbatch -J ${name} ${dep_flag} -t ${time} -c ${cpus} --mem=${mem} -- ${cmd}"
    echo "99999999"
    return
  fi

  # Build a tiny heredoc to submit (avoids needing an external .sh file per step)
  local jid
  jid=$(sbatch --parsable \
    -A "${SLURM_ACCT}" \
    -J "${name}" \
    -t "${time}" -c "${cpus}" --mem="${mem}" \
    -o "${LOGDIR}/${name}_%j.out" \
    -e "${LOGDIR}/${name}_%j.err" \
    ${dep_flag} \
    --wrap="${cmd}")
  echo "[run_phase4b] submitted ${name}: jobid=${jid}" >&2
  echo "${jid}"
}

# ── Phase 4b.1: decompose ────────────────────────────────────────────────────
CMD_DECOMP="Rscript ${MODULE_DIR}/phase_4_postprocessing/4b_group_proposal/STEP_C01i_decompose.R \
  --candidates ${CANDIDATES} \
  --outdir ${DECOMP_OUT} \
  --tier_max ${TIER_MAX}"
JID_DECOMP=$(submit_job "4b1_decomp" "" "04:00:00" 4 16G "${CMD_DECOMP}")

# ── Phase 4b.2: multi_recomb (depends on 4b.1) ───────────────────────────────
CMD_RECOMB="Rscript ${MODULE_DIR}/phase_4_postprocessing/4b_group_proposal/STEP_C01i_b_multi_recomb.R \
  --candidates ${CANDIDATES} \
  --decomp_dir ${DECOMP_OUT} \
  --outdir ${RECOMB_OUT} \
  --tier_max ${TIER_MAX}"
JID_RECOMB=$(submit_job "4b2_recomb" "${JID_DECOMP}" "03:00:00" 2 8G "${CMD_RECOMB}")

# ── Phase 4b.3: nested_composition (parallel with 4b.1) ──────────────────────
CMD_NESTED="python3 ${MODULE_DIR}/phase_4_postprocessing/4b_group_proposal/STEP_C01i_c_nested_composition.py \
  --candidates ${CANDIDATES} \
  --q_cache_dir ${Q_CACHE_DIR} \
  --outdir ${NESTED_OUT} \
  --tier_max ${TIER_MAX}"
JID_NESTED=$(submit_job "4b3_nested" "" "02:00:00" 2 8G "${CMD_NESTED}")

# ── Phase 4b.4: seal (depends on 4b.1, 4b.2, 4b.3) ───────────────────────────
CMD_SEAL="Rscript ${MODULE_DIR}/phase_4_postprocessing/4b_group_proposal/STEP_C01i_d_seal.R \
  --candidates ${CANDIDATES} \
  --decomp_dir ${DECOMP_OUT} \
  --recomb_dir ${RECOMB_OUT} \
  --nested_dir ${NESTED_OUT} \
  --outdir ${SEAL_OUT} \
  --tier_max ${TIER_MAX}"

# Multi-dep: afterok:jid1:jid2:jid3
DEP_SEAL="${JID_DECOMP}:${JID_RECOMB}:${JID_NESTED}"
JID_SEAL=$(submit_job "4b4_seal" "${DEP_SEAL}" "01:00:00" 2 8G "${CMD_SEAL}")

echo ""
echo "[run_phase4b] dependency chain:"
echo "  4b1_decomp  = ${JID_DECOMP}"
echo "  4b2_recomb  = ${JID_RECOMB}  (depends on 4b1)"
echo "  4b3_nested  = ${JID_NESTED}  (parallel with 4b1)"
echo "  4b4_seal    = ${JID_SEAL}  (depends on 4b1, 4b2, 4b3)"
echo ""
echo "[run_phase4b] monitor:"
echo "  squeue -u \$USER -j ${JID_DECOMP},${JID_RECOMB},${JID_NESTED},${JID_SEAL}"
echo ""
echo "[run_phase4b] logs at ${LOGDIR}"
echo "[run_phase4b] outputs at ${MODULE_DIR}/phase4b_out/"
# FIX 48 (chat 10): machine-readable "final job id" line for parent
# orchestrators (run_phase4.sh). The seal job is the exit gate of 4b —
# 4c must gate on it. Emitted as the LAST line of stdout with a fixed
# prefix, so callers can `grep '^PHASE4B_SEAL_JID=' | cut -d= -f2`
# deterministically without parsing anything else.
echo "PHASE4B_SEAL_JID=${JID_SEAL}"
