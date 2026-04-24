#!/bin/bash
# =============================================================================
# STEP_Q10_register.sh
# =============================================================================
# Ingest phase_qc_shelf outputs for one chromosome into the four registries:
#
#   1. interval_registry  — adds an inversion candidate row
#      (chrom, start_bp, end_bp, method="phase_qc_shelf", confidence)
#
#   2. sample_registry    — adds three karyotype groups (Hom1/Het/Hom2)
#      from Q07 invgt_assignments, named inv_<cid>_HOM1 / HET / HOM2
#      with dimension="karyotype" and parent_group="<cohort_group>"
#
#   3. results_registry   — registers the pairwise Fst_Hom1_Hom2 profile
#      (kind=pairwise, stat=fst, who=hom1_vs_hom2, where=<cid>)
#      plus summary stats (Fst mean inside vs outside shelf, n_sites, etc.)
#
#   4. evidence_registry  — writes an evidence bundle:
#      - diagnostic.CHR.pdf + diagnostic.CHR.clean.pdf
#      - shelf_heatmap.CHR.pdf
#      - gap_features.CHR.tsv
#      - popstats_invgt.CHR.tsv
#      - invgt_assignments.CHR.tsv
#      - summary.json (computed stats)
#
# All four writes are idempotent: if the cid already exists, the bridge
# will UPDATE in place (append-only audit log kept in the registries).
#
# Usage:
#   bash STEP_Q10_register.sh <CHR>
#   bash STEP_Q10_register.sh ALL
#
# Env vars:
#   SHELF_START_MB / SHELF_END_MB / BP1_MB / BP2_MB   inversion coords
#   SAMPLE_GROUP     parent cohort group (default: all_226)
#   METHOD_TAG       method name in registry (default: phase_qc_shelf)
#   CID_PREFIX       candidate id prefix (default: inv)
#   FORCE_OVERWRITE  =1 to overwrite existing cid entries (default: 0)
#
# Bridge script calls the R API from inversion-popgen-toolkit/registries/api/R/
# which provides load_sample_registry, load_evidence_registry,
# load_interval_registry, load_results_registry.
# =============================================================================
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/00_config.sh"

CHR="${1:-}"
[[ -z "${CHR}" ]] && qc_die "Usage: $0 <CHR|ALL>"

# Registry paths — prefer env, else derive from BASE
: "${REGISTRIES:=${BASE}/inversion-popgen-toolkit/registries}"
: "${REGISTRY_API_R:=${REGISTRIES}/api/R}"
: "${SAMPLE_GROUP:=all_226}"
: "${METHOD_TAG:=phase_qc_shelf}"
: "${CID_PREFIX:=inv}"
: "${FORCE_OVERWRITE:=0}"

# ALL mode: iterate over every chromosome that has a gap_features.CHR.tsv
if [[ "${CHR}" == "ALL" ]]; then
  qc_log "Q10 ALL: registering every chromosome with phase_qc_shelf outputs"
  registered=0
  skipped=0
  for f in "${QC_TRACKS}"/invgt_assignments.*.tsv; do
    [[ -f "${f}" ]] || continue
    chr=$(basename "${f}" | sed 's/^invgt_assignments\.\(.*\)\.tsv$/\1/')
    if [[ -z "${chr}" || "${chr}" == "ALL" ]]; then continue; fi
    qc_log "--- Q10 ${chr} ---"
    if bash "${BASH_SOURCE[0]}" "${chr}"; then
      registered=$((registered + 1))
    else
      skipped=$((skipped + 1))
      qc_log "  Q10 ${chr}: FAILED (continuing with next chrom)"
    fi
  done
  qc_log "Q10 ALL: done. registered=${registered}, skipped=${skipped}"
  exit 0
fi

# Per-chromosome registration
qc_log "Q10 ${CHR}: registering into registries at ${REGISTRIES}"

precomp="${PRECOMP_DIR}/${CHR}.precomp.rds"
[[ -f "${precomp}" ]] || precomp="${PRECOMP_DIR}/main_qcpass.${CHR}.precomp.rds"
[[ -f "${precomp}" ]] || { qc_log "SKIP ${CHR}: no precomp"; exit 0; }

invgt_file="${QC_TRACKS}/invgt_assignments.${CHR}.tsv"
popstats_file="${QC_TRACKS}/popstats_invgt.${CHR}.tsv"
gap_file="${QC_TRACKS}/gap_features.${CHR}.tsv"
hobs_file="${QC_TRACKS}/hobs_merged.${CHR}.tsv.gz"

[[ -f "${invgt_file}" ]]   || { qc_log "SKIP ${CHR}: no invgt_assignments"; exit 0; }
[[ -f "${popstats_file}" ]] || { qc_log "SKIP ${CHR}: no popstats_invgt"; exit 0; }

# Collect evidence files (exist-or-skip gracefully in the R worker)
evidence_files=()
for f in \
  "${QC_FIGS}/diagnostic.${CHR}.pdf" \
  "${QC_FIGS}/diagnostic.${CHR}.clean.pdf" \
  "${QC_FIGS}/shelf_heatmap.${CHR}.pdf" \
  "${popstats_file}" \
  "${invgt_file}" \
  "${gap_file}" \
  "${hobs_file}"; do
  [[ -f "${f}" ]] && evidence_files+=( "${f}" )
done

qc_log "  ${#evidence_files[@]} evidence files collected"

args=(
  --chrom             "${CHR}"
  --registries_root   "${REGISTRIES}"
  --registry_api_r    "${REGISTRY_API_R}"
  --precomp           "${precomp}"
  --invgt_file        "${invgt_file}"
  --popstats_file     "${popstats_file}"
  --sample_group      "${SAMPLE_GROUP}"
  --method_tag        "${METHOD_TAG}"
  --cid_prefix        "${CID_PREFIX}"
  --force_overwrite   "${FORCE_OVERWRITE}"
)
[[ -n "${SHELF_START_MB:-}" ]] && args+=( --shelf_start_mb "${SHELF_START_MB}" )
[[ -n "${SHELF_END_MB:-}"   ]] && args+=( --shelf_end_mb   "${SHELF_END_MB}"   )
[[ -n "${BP1_MB:-}"         ]] && args+=( --breakpoint1_mb "${BP1_MB}"         )
[[ -n "${BP2_MB:-}"         ]] && args+=( --breakpoint2_mb "${BP2_MB}"         )
[[ -f "${gap_file}" ]]         && args+=( --gap_file       "${gap_file}"       )
[[ -f "${hobs_file}" ]]        && args+=( --hobs_file      "${hobs_file}"      )

# Pass evidence files as a single colon-separated list
if (( ${#evidence_files[@]} > 0 )); then
  IFS=':' ev="${evidence_files[*]}"
  args+=( --evidence_files "${ev}" )
fi

${RSCRIPT_BIN} --vanilla "${here}/R/q10_register.R" "${args[@]}" \
  2> "${QC_LOGS}/q10_${CHR}.log"

qc_log "Q10 ${CHR}: DONE"
qc_log "  log: ${QC_LOGS}/q10_${CHR}.log"
