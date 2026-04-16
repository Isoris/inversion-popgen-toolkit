#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=16G
#SBATCH -t 0-01:00:00
#SBATCH -J tra_prep
#SBATCH -o logs/01_prep.%j.out
#SBATCH -e logs/01_prep.%j.err
set -euo pipefail
source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${SCRIPT_DIR}/../00_module4f_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
set -a
source "${CONFIG}"
set +a

dv_init_dirs
dv_log "=== STEP 1: Prepare inputs for DELLY TRA pipeline ==="

# ── Validate reused inputs from DEL pipeline ───────────────────────────────
dv_check_file "$REF" "Reference FASTA"
dv_check_file "$REF_FAI" "Reference index"
dv_check_file "$EXCL_BED" "Exclusion BED (from DEL pipeline)"
dv_check_file "$SAMPLES_ALL" "Sample list (all 226)"
dv_check_file "$SAMPLES_UNRELATED" "Sample list (unrelated 81)"
[[ -x "${DELLY_BIN}" ]] || dv_die "DELLY binary not executable: ${DELLY_BIN}"

# ── Verify markdup BAMs exist ──────────────────────────────────────────────
dv_log "Verifying markdup BAMs in ${DIR_MARKDUP}..."
MISSING=0
while IFS= read -r sid; do
  bam="${DIR_MARKDUP}/${sid}.markdup.bam"
  [[ -f "${bam}" && -f "${bam}.bai" ]] || { dv_err "Missing: ${bam}"; ((MISSING++)) || true; }
done < "$SAMPLES_ALL"
[[ ${MISSING} -eq 0 ]] || dv_die "${MISSING} markdup BAMs missing. Run DEL pipeline step 01 first."

N_ALL=$(wc -l < "$SAMPLES_ALL")
N_UNREL=$(wc -l < "$SAMPLES_UNRELATED")
dv_log "  Markdup BAMs: ${N_ALL} OK"
dv_log "  Samples all: ${N_ALL}"
dv_log "  Samples unrelated: ${N_UNREL}"

# ── Sort annotation BEDs ──────────────────────────────────────────────────
dv_log "Preparing annotation BEDs..."
ANNOT_OUT="${DIR_ANNOT}/beds"
mkdir -p "${ANNOT_OUT}"

for bed_label in GENE EXON CDS; do
  src="${ANNOT_DIR}/features/${bed_label}.bed"
  if [[ -f "$src" ]]; then
    bedtools sort -i "$src" > "${ANNOT_OUT}/${bed_label}.sorted.bed"
    dv_log "  ${bed_label}: $(wc -l < "${ANNOT_OUT}/${bed_label}.sorted.bed") features"
  else
    dv_log "  WARNING: ${bed_label}.bed not found at ${src}"
  fi
done

if [[ -n "${REPEAT_BED}" && -f "${REPEAT_BED}" ]]; then
  bedtools sort -i "${REPEAT_BED}" > "${ANNOT_OUT}/REPEATS.sorted.bed"
  dv_log "  REPEATS: $(wc -l < "${ANNOT_OUT}/REPEATS.sorted.bed") intervals"
fi

DELLY_VER=$("${DELLY_BIN}" 2>&1 | grep -oP 'Version: \K[0-9.]+' || echo "unknown")
dv_log "  DELLY ${DELLY_VER} at ${DELLY_BIN}"
dv_log "=== STEP 1 COMPLETE ==="
