#!/usr/bin/env bash
# =============================================================================
# LAUNCH_q_residual_dosage.slurm — Export Q-corrected residual dosage
#
# Per-chromosome array (1-28). For each chromosome, subtracts
# E[dosage|Q,F] from observed dosage, producing a residual BEAGLE file
# that downstream tools (ngsLD, lostruct, plot_ld_panels.R) can consume.
#
# Usage:
#   sbatch Modules/MODULE_5B_inversion_followup/launchers/LAUNCH_q_residual_dosage.slurm
# =============================================================================
#SBATCH --job-name=q_residual
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-28
#SBATCH --output=logs/q_residual_%a.out
#SBATCH --error=logs/q_residual_%a.err

set -euo pipefail

# v12.2: Launcher relocated to Modules/MODULE_5B_inversion_followup/launchers/
# The ancestry config lives at the canonical unified_ancestry/ location.
# Allow override via $CONFIG env var; otherwise use ${BASE} from environment.
: "${BASE:=/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"
CONFIG="${CONFIG:-${BASE}/unified_ancestry/00_ancestry_config.sh}"
[[ -f "$CONFIG" ]] || { echo "Config not found: $CONFIG" >&2; exit 1; }
source "$CONFIG"

eval "$(conda shell.bash hook)"
conda activate "${CONDA_ENV:-assembly}"
mkdir -p logs

# ── Build binary if needed ──
# QRES_BIN and MODULE_5B_ENGINES_DIR now come from 00_ancestry_config.sh (v12.2)
if [[ ! -x "$QRES_BIN" ]]; then
  anc_log "Compiling export_q_residual_dosage..."
  cd "${MODULE_5B_ENGINES_DIR}" && make export_q_residual_dosage && cd -
fi

# ── Output directory ──
OUTDIR="${BASE}/unified_ancestry/q_residual_dosage"
mkdir -p "$OUTDIR"

# ── Map array index to chromosome ──
mapfile -t CHROMS < <(awk '{print $1}' "$REF_FAI" | sort -V)
N_CHROMS=${#CHROMS[@]}
IDX=$((SLURM_ARRAY_TASK_ID - 1))
if [[ $IDX -ge $N_CHROMS ]]; then exit 0; fi
CHR="${CHROMS[$IDX]}"

# ── Check inputs ──
BGL=$(ls "${BEAGLE_DIR}"/*${CHR}*.beagle.gz 2>/dev/null | head -1)
[[ -n "$BGL" && -f "$BGL" ]] || { anc_log "No BEAGLE for ${CHR}"; exit 0; }
[[ -f "$BEST_FOPT" ]] || { anc_log "No F matrix: $BEST_FOPT"; exit 1; }
[[ -f "$BEST_QOPT" ]] || { anc_log "No Q matrix: $BEST_QOPT"; exit 1; }

# ── Skip if exists ──
OUT_FILE="${OUTDIR}/${CHR}.residual.dosage.tsv.gz"
if [[ -f "$OUT_FILE" && "${FORCE:-0}" -ne 1 ]]; then
  anc_log "${CHR} exists — skip"
  exit 0
fi

anc_log "Q-residual dosage: ${CHR}"
"$QRES_BIN" \
  --beagle "$BGL" \
  --fopt "$BEST_FOPT" \
  --local_q "$BEST_QOPT" \
  --sample_list "$SAMPLE_LIST" \
  --chr "$CHR" \
  --out_prefix "${OUTDIR}/${CHR}"

anc_log "Done: ${CHR}"
