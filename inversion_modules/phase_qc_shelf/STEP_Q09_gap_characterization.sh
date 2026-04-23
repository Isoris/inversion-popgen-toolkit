#!/bin/bash
# =============================================================================
# STEP_Q09_gap_characterization.sh
# =============================================================================
# For each chromosome, detect zero-SNP / low-density "dropout" windows from
# the Engine B precompute and annotate each region with:
#   - percent N bases (reference assembly gap)
#   - percent soft-masked bases (repeat masker output)
#   - percent uppercase (confident non-repeat) bases
#   - mean mosdepth coverage if available
#   - n samples with local coverage < threshold
#
# Two outputs per chromosome:
#   - ${QC_TRACKS}/gap_features.<CHR>.tsv      per-gap annotated features
#   - ${QC_TRACKS}/dropouts.<CHR>.bed          BED of all gap regions
#
# Two aggregate outputs:
#   - ${QC_TRACKS}/gap_features.ALL.tsv         concatenated features
#   - ${QC_TRACKS}/dropouts.ALL.bed             concatenated BED
#
# Usage:
#   bash STEP_Q09_gap_characterization.sh <CHR>
#   bash STEP_Q09_gap_characterization.sh ALL        # stitch concatenated
# =============================================================================
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/00_config.sh"

CHR="${1:-}"
[[ -z "${CHR}" ]] && qc_die "Usage: $0 <CHR|ALL>"

# Global stitch mode
if [[ "${CHR}" == "ALL" ]]; then
  qc_log "Q09 ALL: stitching per-chromosome gap tables"
  out_tsv="${QC_TRACKS}/gap_features.ALL.tsv"
  out_bed="${QC_TRACKS}/dropouts.ALL.bed"

  first=1
  : > "${out_tsv}.tmp"
  : > "${out_bed}.tmp"
  for f in "${QC_TRACKS}"/gap_features.*.tsv; do
    [[ -f "${f}" ]] || continue
    [[ "${f}" == *"ALL"* ]] && continue
    if [[ "${first}" == "1" ]]; then
      cat "${f}" >> "${out_tsv}.tmp"
      first=0
    else
      tail -n +2 "${f}" >> "${out_tsv}.tmp"
    fi
  done
  for f in "${QC_TRACKS}"/dropouts.*.bed; do
    [[ -f "${f}" ]] || continue
    [[ "${f}" == *"ALL"* ]] && continue
    cat "${f}" >> "${out_bed}.tmp"
  done
  mv "${out_tsv}.tmp" "${out_tsv}"
  sort -k1,1 -k2,2n "${out_bed}.tmp" > "${out_bed}"
  rm -f "${out_bed}.tmp"

  n_tsv=$(wc -l < "${out_tsv}" 2>/dev/null || echo 0)
  n_bed=$(wc -l < "${out_bed}" 2>/dev/null || echo 0)
  qc_log "  -> ${out_tsv}  (${n_tsv} lines)"
  qc_log "  -> ${out_bed}  (${n_bed} regions)"
  exit 0
fi

# Per-chromosome characterization
qc_log "Q09 ${CHR}: gap characterization"

precomp="${PRECOMP_DIR}/${CHR}.precomp.rds"
[[ -f "${precomp}" ]] || precomp="${PRECOMP_DIR}/main_qcpass.${CHR}.precomp.rds"
[[ -f "${precomp}" ]] || { qc_log "SKIP ${CHR}: no precomp"; exit 0; }

out_tsv="${QC_TRACKS}/gap_features.${CHR}.tsv"
out_bed="${QC_TRACKS}/dropouts.${CHR}.bed"

args=(
  --precomp "${precomp}"
  --chrom   "${CHR}"
  --out_tsv "${out_tsv}"
  --out_bed "${out_bed}"
)

# Optional inputs — added only if they exist so users without repeats/mosdepth
# still get N-content annotation.
if [[ -n "${REFERENCE_FASTA:-}" && -f "${REFERENCE_FASTA}" ]]; then
  args+=( --reference "${REFERENCE_FASTA}" )
  qc_log "  using reference: ${REFERENCE_FASTA}"
else
  qc_log "  WARN: REFERENCE_FASTA not set or missing; no N-content or softmask stats"
fi

if [[ -n "${REPEAT_BED:-}" && -f "${REPEAT_BED}" ]]; then
  args+=( --repeat_bed "${REPEAT_BED}" )
  qc_log "  using repeat BED: ${REPEAT_BED}"
fi

if [[ -d "${HET_DIR:-}/mosdepth_output" ]]; then
  args+=( --mosdepth_dir "${HET_DIR}/mosdepth_output" )
  qc_log "  using mosdepth dir: ${HET_DIR}/mosdepth_output"
fi

${RSCRIPT_BIN} --vanilla "${here}/R/q09_gap_characterization.R" "${args[@]}" \
  2> "${QC_LOGS}/q09_${CHR}.log"

if [[ -f "${out_tsv}" ]]; then
  n=$(tail -n +2 "${out_tsv}" | wc -l)
  qc_log "  -> ${out_tsv}  (${n} gap regions)"
fi
[[ -f "${out_bed}" ]] && qc_log "  -> ${out_bed}"
qc_log "Q09 ${CHR} DONE."
