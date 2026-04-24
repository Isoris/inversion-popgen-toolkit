#!/bin/bash
# =============================================================================
# STEP_Q07c_hobs_windower.sh
# =============================================================================
# Run the hobs_windower C binary on each of the three per-group .hwe.gz
# outputs from Q07b, then merge the windowed results into a single TSV
# for the chromosome. Computes Hobs, Hexp, and the HoverE ratio at seven
# multi-scale sliding windows.
#
# The hobs_windower binary implements Claire Mérot's sliding-window Hobs
# logic (per-site Hobs = Hexp * (1 - F) from ANGSD .hwe.gz, aggregated
# into windows at multiple scales with mean / median / SD / outlier burden).
#
# At an inversion locus the signature is:
#   HoverE_HET  → 2.0   (every heterokaryotype sample is het at every
#                         differentiated SNP between the two arrangements)
#   HoverE_HOM1 → 0     (homozygous-for-arrangement-1 samples have no
#                         het genotypes at arrangement-differentiating SNPs)
#   HoverE_HOM2 → 0     (same, opposite arrangement)
#
# Outside the inversion: HoverE → 1.0 for all three groups (HWE).
#
# Outputs:
#   ${QC_TRACKS}/hobs_sites.<CHR>.<GROUP>.tsv.gz   (per-site Hobs/Hexp/F)
#   ${QC_TRACKS}/hobs_win.<CHR>.<GROUP>.<SCALE>.tsv.gz   (per-scale windows)
#   ${QC_TRACKS}/hobs_merged.<CHR>.tsv.gz          (wide merged TSV for Q04)
#
# Usage:
#   bash STEP_Q07c_hobs_windower.sh <CHR>
#   bash STEP_Q07c_hobs_windower.sh ALL
#
# Env:
#   HOBS_WINDOWER_BIN  path to compiled hobs_windower
#                      default: ${BASE}/unified_ancestry/engines/hobs_hwe/scripts/hobs_windower
#   HOBS_SCALES        scales to run, default:
#                      "5kb:5000:1000,10kb:10000:2000,50kb:50000:10000,100kb:100000:20000,250kb:250000:50000,500kb:500000:100000,1Mb:1000000:200000"
# =============================================================================
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/00_config.sh"

CHR="${1:-}"
[[ -z "${CHR}" ]] && qc_die "Usage: $0 <CHR|ALL>"

: "${HOBS_WINDOWER_BIN:=${BASE}/unified_ancestry/engines/hobs_hwe/scripts/hobs_windower}"
if [[ ! -x "${HOBS_WINDOWER_BIN}" ]]; then
  for cand in \
    "${BASE}/unified_ancestry/engines/hobs_hwe/hobs_windower" \
    "${BASE}/unified_ancestry/src/hobs_windower" \
    "${BASE}/unified_ancestry/bin/hobs_windower"; do
    if [[ -x "${cand}" ]]; then HOBS_WINDOWER_BIN="${cand}"; break; fi
  done
fi
[[ -x "${HOBS_WINDOWER_BIN}" ]] || qc_die "hobs_windower binary not found (try compiling from unified_ancestry/engines/hobs_hwe/scripts/hobs_windower.c)"

: "${HOBS_SCALES:=5kb:5000:1000,10kb:10000:2000,50kb:50000:10000,100kb:100000:20000,250kb:250000:50000,500kb:500000:100000,1Mb:1000000:200000}"

# ALL mode
if [[ "${CHR}" == "ALL" ]]; then
  CHR_LIST="${BEAGLE_DIR}/chr.list"
  [[ -f "${CHR_LIST}" ]] || qc_die "chr.list not found at ${CHR_LIST}"
  while IFS= read -r c; do
    [[ -z "${c}" ]] && continue
    bash "${BASH_SOURCE[0]}" "${c}"
  done < "${CHR_LIST}"
  exit 0
fi

qc_log "Q07c ${CHR}: windowing per-group Hobs"
qc_log "  binary: ${HOBS_WINDOWER_BIN}"

# Determine chromosome size from reference .fai
: "${REF:=${BASE}/00-samples/fClaHyb_Gar_LG.fa}"
: "${REF_FAI:=${REF}.fai}"
[[ -f "${REF_FAI}" ]] || qc_die "reference .fai not found at ${REF_FAI}"
chrom_size=$(awk -v c="${CHR}" '$1==c {print $2; exit}' "${REF_FAI}")
[[ -n "${chrom_size}" ]] || qc_die "chromosome ${CHR} not in ${REF_FAI}"
qc_log "  chrom_size: ${chrom_size} bp"

hwe_dir="${QC_TRACKS}/hwe_per_group"
win_dir="${QC_TRACKS}/hobs_win"
mkdir -p "${win_dir}"

# Run hobs_windower on each group
run_windower() {
  local group="$1"
  local hwe_file="${hwe_dir}/${CHR}.${group}.hwe.gz"
  local out_prefix="${win_dir}/${CHR}.${group}"
  if [[ ! -f "${hwe_file}" ]]; then
    qc_log "  ${group}: ${hwe_file} missing, skipping"
    return 0
  fi
  qc_log "  ${group}: windowing ${hwe_file}"
  "${HOBS_WINDOWER_BIN}" \
    "${hwe_file}" \
    "${out_prefix}" \
    "${chrom_size}" \
    --scales "${HOBS_SCALES}" \
    2> "${QC_LOGS}/q07c_${CHR}_${group}.log"
}

for grp in HOM1 HET HOM2; do
  run_windower "${grp}"
done

# Merge windowed outputs into a single wide TSV for Q04/Q10.
# hobs_windower writes <prefix>.win<label>.tsv per scale.
# We pick one primary scale (10kb default) for the main Q04 panel, and keep
# others available for cross-scale checks.
: "${HOBS_PRIMARY_SCALE:=10kb}"
out_merged="${QC_TRACKS}/hobs_merged.${CHR}.tsv.gz"
qc_log "  merging scale=${HOBS_PRIMARY_SCALE} -> ${out_merged}"

${RSCRIPT_BIN} --vanilla "${here}/R/q07c_merge_hobs.R" \
  --chrom "${CHR}" \
  --win_dir "${win_dir}" \
  --scale "${HOBS_PRIMARY_SCALE}" \
  --out "${out_merged}" \
  2> "${QC_LOGS}/q07c_merge_${CHR}.log"

[[ -f "${out_merged}" ]] && qc_log "  -> ${out_merged} ($(zcat "${out_merged}" | wc -l) rows)"
qc_log "Q07c ${CHR} DONE"
