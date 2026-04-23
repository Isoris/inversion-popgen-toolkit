#!/bin/bash
# =============================================================================
# STEP_Q07b_hobs_per_group.sh
# =============================================================================
# Run patched ANGSD (Isoris/angsd_fixed_HWE) with -doHWE 1 separately on each
# karyotype group (Hom1 / Het / Hom2) from Q07's invgt_assignments. Emits
# three .hwe.gz files per chromosome, one per group.
#
# Full Mérot approach: F is estimated independently per group, so Hexp at
# an inversion SNP is computed from within-group allele frequencies. This is
# the scientifically correct denominator for inversion-informative Hobs.
#
# Outputs:
#   ${QC_TRACKS}/hwe_per_group/<CHR>.<GROUP>.hwe.gz     (.hwe.gz per group)
#   ${QC_TRACKS}/hwe_per_group/<CHR>.<GROUP>.arg        (ANGSD log for audit)
#
# Usage:
#   bash STEP_Q07b_hobs_per_group.sh <CHR>
#   bash STEP_Q07b_hobs_per_group.sh ALL       # iterates chr.list
#
# Env:
#   ANGSD_PATCHED_BIN  path to Isoris/angsd_fixed_HWE (default: $BASE/angsd_fixed_HWE/angsd)
#   BAM_DIR            where per-sample BAMs live (from 00_config.sh)
#   HWE_THREADS        cpus per ANGSD call (default 2 — we run 3 in parallel)
#   HWE_MAJORMINOR     -doMajorMinor mode (default 1, from allele frequencies)
#   HWE_MINMAF         min MAF filter (default 0.05)
#   HWE_MINMAPQ        min mapping quality (default 20)
#   HWE_MINQ           min base quality (default 20)
# =============================================================================
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/00_config.sh"

CHR="${1:-}"
[[ -z "${CHR}" ]] && qc_die "Usage: $0 <CHR|ALL>"

# Patched ANGSD binary
: "${ANGSD_PATCHED_BIN:=${BASE}/angsd_fixed_HWE/angsd}"
if [[ ! -x "${ANGSD_PATCHED_BIN}" ]]; then
  # Fall back to a couple common locations before giving up
  for cand in \
    "${BASE}/angsd_fixed_HWE/angsd" \
    "${BASE}/programs/angsd_fixed_HWE/angsd" \
    "/lustrefs/disk/project/lt200308-agbsci/13-programs/angsd_fixed_HWE/angsd"; do
    if [[ -x "${cand}" ]]; then ANGSD_PATCHED_BIN="${cand}"; break; fi
  done
fi
[[ -x "${ANGSD_PATCHED_BIN}" ]] || qc_die "patched ANGSD not found (looked at ${ANGSD_PATCHED_BIN})"

# ALL mode — iterate
if [[ "${CHR}" == "ALL" ]]; then
  CHR_LIST="${BEAGLE_DIR}/chr.list"
  [[ -f "${CHR_LIST}" ]] || qc_die "chr.list not found at ${CHR_LIST}"
  while IFS= read -r c; do
    [[ -z "${c}" ]] && continue
    bash "${BASH_SOURCE[0]}" "${c}"
  done < "${CHR_LIST}"
  exit 0
fi

# Per-chromosome
qc_log "Q07b ${CHR}: per-group ANGSD -doHWE using patched binary"
qc_log "  ANGSD: ${ANGSD_PATCHED_BIN}"

invgt_file="${QC_TRACKS}/invgt_assignments.${CHR}.tsv"
[[ -f "${invgt_file}" ]] || { qc_log "SKIP ${CHR}: no invgt_assignments (run Q07 first)"; exit 0; }

out_dir="${QC_TRACKS}/hwe_per_group"
mkdir -p "${out_dir}"

# Resolve sample→BAM mapping: each line of SAMPLE_LIST_POPSTATS is a sample_id.
# Per-sample BAMs live at ${BAM_DIR}/<sample>.${BAM_SUFFIX:-markdup.bam}.
: "${BAM_SUFFIX:=markdup.bam}"
build_bamlist() {
  local group="$1" out="$2"
  : > "${out}"
  awk -F'\t' -v g="${group}" 'NR>1 && $2==g {print $1}' "${invgt_file}" | \
    while IFS= read -r sid; do
      # Normalize group casing (Hom1 / HOM1 / hom1 all accepted)
      :
      for bam in \
        "${BAM_DIR}/${sid}.${BAM_SUFFIX}" \
        "${BAM_DIR}/${sid}.bam" \
        "${BAM_DIR}/${sid}.markdup.bam"; do
        if [[ -f "${bam}" ]]; then
          echo "${bam}" >> "${out}"
          break
        fi
      done
    done
}

# Normalize invgt labels so we accept Hom1/HOM1/hom_ref
# (Q07 writes whatever the downstream R script emits; normalize for safety)
tmp_invgt="${out_dir}/${CHR}.invgt_normalized.tsv"
awk -F'\t' 'BEGIN {OFS="\t"}
  NR==1 {print; next}
  {
    g = toupper($2)
    gsub(/^HOM_REF$|^HOMOZYGOUS_REF$|^REF$/, "HOM1", g)
    gsub(/^HETEROZYGOUS$/, "HET", g)
    gsub(/^HOM_ALT$|^HOMOZYGOUS_ALT$|^INV$/, "HOM2", g)
    $2 = g
    print
  }' "${invgt_file}" > "${tmp_invgt}"

# Spawn per-group ANGSD calls in parallel (3 groups, 2 cpus each = 6 cpus total)
: "${HWE_THREADS:=2}"
: "${HWE_MAJORMINOR:=1}"
: "${HWE_MINMAF:=0.05}"
: "${HWE_MINMAPQ:=20}"
: "${HWE_MINQ:=20}"

run_group() {
  local group="$1"
  local bamlist="${out_dir}/${CHR}.${group}.bamlist"
  local out_prefix="${out_dir}/${CHR}.${group}"

  build_bamlist "${group}" "${bamlist}"
  local n_bams
  n_bams=$(wc -l < "${bamlist}")
  if (( n_bams == 0 )); then
    qc_log "  ${group}: no BAMs resolved, skipping"
    return 0
  fi
  qc_log "  ${group}: ${n_bams} samples  ->  ${out_prefix}.hwe.gz"

  # Skip if output already exists and is non-empty (idempotent rerun)
  if [[ -f "${out_prefix}.hwe.gz" ]] && \
     [[ $(zcat "${out_prefix}.hwe.gz" 2>/dev/null | wc -l) -gt 1 ]] && \
     [[ "${FORCE_Q07B:-0}" != "1" ]]; then
    qc_log "  ${group}: cached output exists, skipping (FORCE_Q07B=1 to override)"
    return 0
  fi

  "${ANGSD_PATCHED_BIN}" \
    -bam "${bamlist}" \
    -ref "${REF:-${BASE}/00-samples/fClaHyb_Gar_LG.fa}" \
    -out "${out_prefix}" \
    -GL 1 \
    -doMajorMinor "${HWE_MAJORMINOR}" \
    -doMaf 1 \
    -SNP_pval 1e-6 \
    -doHWE 1 \
    -maxHetFreq 1.0 \
    -minMaf "${HWE_MINMAF}" \
    -minMapQ "${HWE_MINMAPQ}" \
    -minQ   "${HWE_MINQ}" \
    -r "${CHR}:" \
    -nThreads "${HWE_THREADS}" \
    2> "${QC_LOGS}/q07b_${CHR}_${group}.log"

  if [[ -f "${out_prefix}.hwe.gz" ]]; then
    local n_sites
    n_sites=$(zcat "${out_prefix}.hwe.gz" | tail -n +2 | wc -l)
    qc_log "  ${group}: done, ${n_sites} sites"
  else
    qc_log "  ${group}: FAILED — no .hwe.gz produced. See ${QC_LOGS}/q07b_${CHR}_${group}.log"
  fi
}

# Parallel launch
run_group HOM1 &
pid_hom1=$!
run_group HET &
pid_het=$!
run_group HOM2 &
pid_hom2=$!

wait "${pid_hom1}" "${pid_het}" "${pid_hom2}"

qc_log "Q07b ${CHR} DONE"
