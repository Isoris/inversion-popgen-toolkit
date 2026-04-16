#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 127
#SBATCH --mem=237G
#SBATCH -t 1-00:00:00
#SBATCH -J manta_prep
#SBATCH -o logs/01_prep.%j.out
#SBATCH -e logs/01_prep.%j.err
# =============================================================================
# 01_prep_inputs.sh — Prepare all inputs for Manta SV pipeline
# =============================================================================
# This script does:
#   1. Loads config (SLURM-safe)
#   2. Reads BAM manifest
#   3. Verifies markdup BAMs exist (reuses DELLY markdup)
#      — If missing, runs markdup (same pipeline as DELLY 01_prep)
#   4. Builds custom Manta config ini (minCandidateVariantSize = 50)
#   5. Builds call_regions BED (inverse of exclude BED, bgzipped + tabixed)
#   6. Builds sample lists (if not already present from DELLY)
#   7. Sorts annotation BEDs
# =============================================================================
set -euo pipefail
source ~/.bashrc
mamba activate manta_py2
export PATH="/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin:$PATH"

# ── Config loading (SLURM-safe) ────────────────────────────────────────────
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${SCRIPT_DIR}/../00_module4g_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
set -a
source "${CONFIG}"
set +a

mv_init_dirs
mv_log "=== STEP 1: Prepare inputs for Manta SV pipeline ==="

# ─────────────────────────────────────────────────────────────────────────────
# 1A. Read BAM manifest
# ─────────────────────────────────────────────────────────────────────────────
mv_log "--- 1A: Reading BAM manifest ---"
mv_read_manifest

mv_log "  First 5 samples:"
for i in 0 1 2 3 4; do
  [[ $i -lt ${MANIFEST_N_SAMPLES} ]] || break
  mv_log "    ${MANIFEST_SIDS[$i]}  →  $(basename "${MANIFEST_RAW_BAMS[$i]}")"
done

# ─────────────────────────────────────────────────────────────────────────────
# 1B. Verify markdup BAMs (reuse from DELLY)
# ─────────────────────────────────────────────────────────────────────────────
mv_log "--- 1B: Verifying markdup BAMs ---"
mv_log "  Markdup directory: ${DIR_MARKDUP}"

if [[ -d "${DIR_MARKDUP}" ]]; then
  MISSING_MARKDUP=0
  for i in $(seq 0 $((MANIFEST_N_SAMPLES - 1))); do
    sid="${MANIFEST_SIDS[$i]}"
    mbam="${DIR_MARKDUP}/${sid}.markdup.bam"
    if [[ ! -f "${mbam}" || ! -f "${mbam}.bai" ]]; then
      mv_err "Missing markdup BAM or index for: ${sid}"
      ((MISSING_MARKDUP++)) || true
    fi
  done

  if [[ ${MISSING_MARKDUP} -gt 0 ]]; then
    mv_die "${MISSING_MARKDUP} markdup BAMs missing. Run DELLY 01_prep first or set DIR_MARKDUP to create fresh."
  fi
  mv_log "  All ${MANIFEST_N_SAMPLES} markdup BAMs verified from DELLY."
else
  mv_die "Markdup directory not found: ${DIR_MARKDUP}. Run DELLY pipeline first to create markdup BAMs, or update DIR_MARKDUP in config."
fi

# ─────────────────────────────────────────────────────────────────────────────
# 1C. Build custom Manta config ini
# ─────────────────────────────────────────────────────────────────────────────
mv_log "--- 1C: Building custom Manta config ini ---"
mv_log "  minCandidateVariantSize = 50 (default is 8)"

cat > "${MANTA_CUSTOM_INI}" << 'INIEOF'
#
# Custom Manta configuration for catfish lcWGS cohort
#
# Key change: minCandidateVariantSize raised from 8 (default) to 50.
# Rationale: Clair3 handles SNPs + indels ≤50 bp. Manta covers ≥50 bp SVs.
# This avoids redundant small indel calls and speeds up runtime.
#

[manta]

# Minimum size for SV/indel candidates (default: 8)
minCandidateVariantSize = 50

# RNA min candidate size (not relevant for DNA, kept at default)
rnaMinCandidateVariantSize = 1000

# Minimum graph edge observations (default: 3)
# Keep default — appropriate for ~9X lcWGS
minEdgeObservations = 3

# Graph node max edge count (default: 10)
graphNodeMaxEdgeCount = 10

# Minimum spanning support for candidates (default: 3)
minCandidateSpanningCount = 3

# Only score/report SVs at or above this size (default: 50)
minScoredVariantSize = 50

# Minimum QUAL for diploid VCF output (default: 10)
minDiploidVariantScore = 10

# Minimum QUAL for somatic VCF output (not used, kept default)
minSomaticScore = 10
INIEOF

mv_log "  Written: ${MANTA_CUSTOM_INI}"

# ─────────────────────────────────────────────────────────────────────────────
# 1D. Build call regions BED (inverse of exclude BED)
# ─────────────────────────────────────────────────────────────────────────────
# Manta uses --callRegions (bgzip+tabix BED) to restrict calling.
# We invert the DELLY exclude BED to create callable regions.
# ─────────────────────────────────────────────────────────────────────────────
mv_log "--- 1D: Building call regions BED ---"
mv_check_file "${REF_FAI}" "Reference FASTA index"
mv_check_file "${EXCL_BED}" "Exclude BED (from DELLY)"

CALL_BED_RAW="${OUTDIR}/call_regions.bed"

# Use bedtools complement to invert exclude → callable
# The genome file for bedtools is just chrom\tsize from the .fai
GENOME_FILE="${OUTDIR}/genome.bed"
awk -F'\t' '{print $1"\t"$2}' "${REF_FAI}" > "${GENOME_FILE}"

# Sort the exclude BED first (bedtools requires sorted input)
EXCL_SORTED="${OUTDIR}/exclude.sorted.bed"
bedtools sort -i "${EXCL_BED}" -g "${GENOME_FILE}" > "${EXCL_SORTED}"

bedtools complement -i "${EXCL_SORTED}" -g "${GENOME_FILE}" > "${CALL_BED_RAW}"

# bgzip + tabix
bgzip -c "${CALL_BED_RAW}" > "${CALL_REGIONS_BED_GZ}"
tabix -p bed "${CALL_REGIONS_BED_GZ}"

N_CALL=$(wc -l < "${CALL_BED_RAW}")
CALL_BP=$(awk '{s+=$3-$2}END{print s}' "${CALL_BED_RAW}")
GENOME_BP=$(awk '{s+=$2}END{print s}' "${GENOME_FILE}")
mv_log "  Call regions: ${N_CALL} intervals, ${CALL_BP} bp callable ($(echo "scale=1; ${CALL_BP}*100/${GENOME_BP}" | bc)% of genome)"
mv_log "  Written: ${CALL_REGIONS_BED_GZ}"

# ─────────────────────────────────────────────────────────────────────────────
# 1E. Verify / build sample lists
# ─────────────────────────────────────────────────────────────────────────────
mv_log "--- 1E: Verifying sample lists ---"
if [[ -f "${SAMPLES_ALL}" ]]; then
  N_ALL=$(wc -l < "${SAMPLES_ALL}")
  mv_log "  Reusing samples_all: ${N_ALL} from ${SAMPLES_ALL}"
else
  printf '%s\n' "${MANIFEST_SIDS[@]}" | sort > "${SAMPLES_ALL}"
  N_ALL=$(wc -l < "${SAMPLES_ALL}")
  mv_log "  Built samples_all: ${N_ALL}"
fi

if [[ -f "${SAMPLES_UNRELATED}" ]]; then
  N_UNREL=$(wc -l < "${SAMPLES_UNRELATED}")
  mv_log "  Reusing samples_unrelated: ${N_UNREL} from ${SAMPLES_UNRELATED}"
else
  mv_check_file "${NATORA_KEEP}" "NAToRA keep list"
  cp "${NATORA_KEEP}" "${SAMPLES_UNRELATED}"
  N_UNREL=$(wc -l < "${SAMPLES_UNRELATED}")
  mv_log "  Built samples_unrelated: ${N_UNREL}"
fi

# ─────────────────────────────────────────────────────────────────────────────
# 1F. Sort annotation BEDs
# ─────────────────────────────────────────────────────────────────────────────
mv_log "--- 1F: Preparing annotation BEDs ---"
ANNOT_OUT="${DIR_ANNOT}/beds"
mkdir -p "${ANNOT_OUT}"

for bed_label in GENE EXON CDS; do
  src="${ANNOT_DIR}/features/${bed_label}.bed"
  if [[ -f "$src" ]]; then
    bedtools sort -i "$src" > "${ANNOT_OUT}/${bed_label}.sorted.bed"
    mv_log "  ${bed_label}: $(wc -l < "${ANNOT_OUT}/${bed_label}.sorted.bed") features"
  else
    mv_log "  WARNING: ${bed_label}.bed not found at ${src}"
  fi
done

if [[ -n "${REPEAT_BED}" && -f "${REPEAT_BED}" ]]; then
  bedtools sort -i "${REPEAT_BED}" > "${ANNOT_OUT}/REPEATS.sorted.bed"
  mv_log "  REPEATS: $(wc -l < "${ANNOT_OUT}/REPEATS.sorted.bed") intervals"
fi

# ─────────────────────────────────────────────────────────────────────────────
# 1G. Verify Manta installation
# ─────────────────────────────────────────────────────────────────────────────
mv_log "--- 1G: Checking Manta installation ---"
[[ -f "${MANTA_CONFIG}" ]] || mv_die "configManta.py not found: ${MANTA_CONFIG}"
[[ -f "${MANTA_CONVERT_INV}" ]] || mv_die "convertInversion.py not found: ${MANTA_CONVERT_INV}"

# Check Python 2 availability (Manta requires it)
if command -v python2 &>/dev/null; then
  PY2=$(command -v python2)
elif python --version 2>&1 | grep -q 'Python 2'; then
  PY2=$(command -v python)
else
  mv_log "  WARNING: Python 2 not found. Manta may fail. Ensure python2 is in PATH."
  PY2="python"
fi
mv_log "  Manta: ${MANTA_INSTALL}"
mv_log "  Python for Manta: ${PY2}"
mv_log "  convertInversion.py: ${MANTA_CONVERT_INV}"

# ─────────────────────────────────────────────────────────────────────────────
# Done
# ─────────────────────────────────────────────────────────────────────────────
mv_log "=== STEP 1 COMPLETE ==="
mv_log "  Markdup BAMs:    ${DIR_MARKDUP}/ (${MANIFEST_N_SAMPLES} samples)"
mv_log "  Custom ini:      ${MANTA_CUSTOM_INI}"
mv_log "  Call regions:    ${CALL_REGIONS_BED_GZ} (${N_CALL} intervals)"
mv_log "  Samples all:     ${SAMPLES_ALL} (${N_ALL})"
mv_log "  Samples unrel:   ${SAMPLES_UNRELATED} (${N_UNREL})"
mv_log "Next: submit 02_manta_discovery.sh"
