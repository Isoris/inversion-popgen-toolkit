#!/usr/bin/env bash
# =============================================================================
# LAUNCH_instant_q_precompute.slurm — SLURM array for Engine B (v12.1)
#
# v12.1: Sources pipeline_bridge.sh for cross-module env vars
# =============================================================================
#SBATCH --job-name=iq_precomp
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-28
#SBATCH --output=logs/iq_precomp_%a.out
#SBATCH --error=logs/iq_precomp_%a.err

set -euo pipefail

# ── Load config ──
CONFIG="${CONFIG:-$(dirname "${BASH_SOURCE[0]}")/../00_ancestry_config.sh}"
[[ -f "$CONFIG" ]] || { echo "Config not found: $CONFIG" >&2; exit 1; }
source "$CONFIG"

# v12.1: Source pipeline bridge
BRIDGE_SH="${BASE}/inversion_codebase_v8.5/utils/pipeline_bridge.sh"
[[ -f "$BRIDGE_SH" ]] && source "$BRIDGE_SH"

eval "$(conda shell.bash hook)"
conda activate "${CONDA_ENV:-assembly}"

# ── Map array index to chromosome ──
mapfile -t CHROMS < <(awk '{print $1}' "$REF_FAI" | sort -V)
N_CHROMS=${#CHROMS[@]}

IDX=$((SLURM_ARRAY_TASK_ID - 1))
if [[ $IDX -ge $N_CHROMS ]]; then
  echo "[LAUNCH] Array index ${SLURM_ARRAY_TASK_ID} exceeds $N_CHROMS — skip"
  exit 0
fi

CHR="${CHROMS[$IDX]}"
anc_log "Task ${SLURM_ARRAY_TASK_ID}/${N_CHROMS}: chromosome ${CHR}"

# ── Resolve BEAGLE ──
BGL=$(ls "${BEAGLE_DIR}"/*${CHR}*.beagle.gz 2>/dev/null | head -1)
if [[ -z "$BGL" || ! -f "$BGL" ]]; then
  anc_log "No BEAGLE for ${CHR} — skip"
  exit 0
fi

# ── Check binary ──
if [[ ! -x "$INSTANT_Q_BIN" ]]; then
  anc_log "Compiling instant_q..."
  cd "$(dirname "$INSTANT_Q_BIN")"
  make -j "${SLURM_CPUS_PER_TASK:-4}"
  cd -
fi

mkdir -p "${LOCAL_Q_DIR}" logs

# ── Check cache ──
SUMMARY="${LOCAL_Q_DIR}/${CHR}.local_Q_summary.tsv.gz"
if [[ -f "$SUMMARY" && "${FORCE:-0}" -ne 1 ]]; then
  anc_log "Cache exists: ${SUMMARY} — skip (set FORCE=1 to recompute)"
  exit 0
fi

# ── Run precompute ──
TMPDIR_IQ=$(mktemp -d /tmp/iq_${CHR}_XXXXXX)

anc_log "Running: ${CHR} (win=${IQ_WINDOW_SIZE}, step=${IQ_WINDOW_STEP}, iter=${IQ_EM_ITER})"
T_START=$(date +%s)

"$INSTANT_Q_BIN" \
  --beagle "$BGL" \
  --fopt "$BEST_FOPT" \
  --qinit "$BEST_QOPT" \
  --precompute \
  --outdir "$TMPDIR_IQ" \
  --chr "$CHR" \
  --window_size "$IQ_WINDOW_SIZE" \
  --window_step "$IQ_WINDOW_STEP" \
  --em_iter "$IQ_EM_ITER" \
  --ncores "${SLURM_CPUS_PER_TASK:-4}" \
  --sample_output

T_END=$(date +%s)
ELAPSED=$((T_END - T_START))

# ── Move to cache ──
for f in "${TMPDIR_IQ}/${CHR}".*; do
  [[ -f "$f" ]] || continue
  bn=$(basename "$f")
  if [[ "$bn" == *summary.tsv || "$bn" == *samples.tsv ]]; then
    gzip -c "$f" > "${LOCAL_Q_DIR}/${bn}.gz"
  else
    cp "$f" "${LOCAL_Q_DIR}/"
  fi
done

rm -rf "$TMPDIR_IQ"
anc_log "Done: ${CHR} in ${ELAPSED}s"

# ── Update manifest ──
MANIFEST="${LOCAL_Q_DIR}/manifest.tsv"
if [[ ! -f "$MANIFEST" ]]; then
  echo -e "chrom\tsummary_file\tsamples_file\ttimestamp" > "$MANIFEST"
fi
echo -e "${CHR}\t${LOCAL_Q_DIR}/${CHR}.local_Q_summary.tsv.gz\t${LOCAL_Q_DIR}/${CHR}.local_Q_samples.tsv.gz\t$(date -Iseconds)" >> "$MANIFEST"

anc_log "Manifest updated"
