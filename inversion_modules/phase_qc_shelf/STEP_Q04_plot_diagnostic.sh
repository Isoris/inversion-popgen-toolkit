#!/bin/bash
# =============================================================================
# STEP_Q04_plot_diagnostic.sh
# =============================================================================
# Read the three tracks (Q01 snp_density, Q02 uncertainty, Q03 coverage) plus
# the precomp Z profile and produce a single stacked 5-track PDF per chromosome
# with shelf regions highlighted.
#
# Usage:
#   SHELF_START_MB=15 SHELF_END_MB=18 bash STEP_Q04_plot_diagnostic.sh C_gar_LG28
#   bash STEP_Q04_plot_diagnostic.sh all     # one PDF per chrom, no shelf shading
# =============================================================================
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/00_config.sh"

CHR="${1:-}"
: "${SHELF_START_MB:=}"
: "${SHELF_END_MB:=}"
[[ -z "${CHR}" ]] && qc_die "Usage: $0 <CHR|all>"

run_one() {
  local chr="$1"
  local out_pdf="${QC_FIGS}/diagnostic.${chr}${OUT_SUFFIX:-}.pdf"

  local precomp="${PRECOMP_DIR}/${chr}.precomp.rds"
  [[ -f "${precomp}" ]] || precomp="${PRECOMP_DIR}/main_qcpass.${chr}.precomp.rds"
  [[ -f "${precomp}" ]] || { qc_log "SKIP ${chr}: no precomp"; return 0; }

  qc_log "Q04 ${chr}: composing diagnostic figure${OUT_SUFFIX:+ ($OUT_SUFFIX)}"

  # Look up theta track (prefer finer scale)
  local theta_file=""
  for scale in win50000.step10000 win10000.step2000 win5000.step1000 win500000.step500000; do
    local cand="${QC_TRACKS}/theta.${chr}.${scale}.tsv.gz"
    if [[ -f "${cand}" ]]; then theta_file="${cand}"; break; fi
  done

  local args=(
    --precomp "${precomp}"
    --chrom "${chr}"
    --snp_track "${QC_TRACKS}/snp_density.${chr}.tsv"
    --unc_track "${QC_TRACKS}/uncertainty.${chr}.tsv"
    --cov_track "${QC_TRACKS}/coverage.${chr}.tsv"
    --ancestry_track "${QC_TRACKS}/ancestry_window.${chr}.tsv"
    --out "${out_pdf}"
  )
  # Multi-scale ancestry
  local anc_multi="${QC_TRACKS}/ancestry_window_multiscale.${chr}.tsv"
  [[ -f "${anc_multi}" ]] && args+=( --ancestry_multiscale "${anc_multi}" )
  # Popstats invgt (Hom1/Het/Hom2 θπ + Fst) from Q07
  local popstats_invgt="${QC_TRACKS}/popstats_invgt.${chr}.tsv"
  [[ -f "${popstats_invgt}" ]] && args+=( --popstats_invgt "${popstats_invgt}" )
  # Invgt assignments (for sample class strip / coloring)
  local invgt_assign="${QC_TRACKS}/invgt_assignments.${chr}.tsv"
  [[ -f "${invgt_assign}" ]] && args+=( --invgt_assignments "${invgt_assign}" )
  # Per-group Hobs/HoverE track from Q07c (Engine H, patched ANGSD)
  local hobs_merged="${QC_TRACKS}/hobs_merged.${chr}.tsv.gz"
  [[ -f "${hobs_merged}" ]] && args+=( --hobs_track "${hobs_merged}" )
  # sim_mat (separate RDS if present; otherwise embedded in precomp)
  local sim_rds_v1="${PRECOMP_DIR}/${chr}.sim_mat.rds"
  local sim_rds_v2="${PRECOMP_DIR}/main_qcpass.${chr}.sim_mat.rds"
  if   [[ -f "${sim_rds_v1}" ]]; then args+=( --sim_mat "${sim_rds_v1}" )
  elif [[ -f "${sim_rds_v2}" ]]; then args+=( --sim_mat "${sim_rds_v2}" )
  fi
  # Reference-N BED (for dark-stipple assembly-gap layer on the ideogram)
  # Prefer a per-chromosome BED in QC_TRACKS, then a genome-wide one.
  local ref_n_bed=""
  if [[ -f "${QC_TRACKS}/ref_n.${chr}.bed" ]]; then
    ref_n_bed="${QC_TRACKS}/ref_n.${chr}.bed"
  elif [[ -n "${REF_N_BED:-}" && -f "${REF_N_BED}" ]]; then
    ref_n_bed="${REF_N_BED}"
  fi
  [[ -n "${ref_n_bed}" ]] && args+=( --ref_n_bed "${ref_n_bed}" )

  [[ -n "${theta_file}" ]] && args+=( --theta_track "${theta_file}" )
  [[ -n "${SHELF_START_MB}" ]] && args+=( --shelf_start_mb "${SHELF_START_MB}" )
  [[ -n "${SHELF_END_MB}"   ]] && args+=( --shelf_end_mb   "${SHELF_END_MB}"   )
  [[ -n "${BREAKPOINT1_MB:-}" ]] && args+=( --breakpoint1_mb "${BREAKPOINT1_MB}" )
  [[ -n "${BREAKPOINT2_MB:-}" ]] && args+=( --breakpoint2_mb "${BREAKPOINT2_MB}" )
  [[ -n "${SMOOTH_WIN:-}"   ]] && args+=( --smooth_win     "${SMOOTH_WIN}"     )
  [[ -n "${SNP_DENSITY_SCALE_KB:-}" ]] && args+=( --snp_density_scale_kb "${SNP_DENSITY_SCALE_KB}" )
  [[ "${Q04_NO_STRIPS:-}" == "1" ]] && args+=( --no_nodata_strips yes )

  ${RSCRIPT_BIN} --vanilla "${here}/R/q04_compose_plot.R" "${args[@]}" \
    2> "${QC_LOGS}/q04_${chr}${OUT_SUFFIX:-}.log"

  qc_log "Q04 ${chr}: wrote ${out_pdf}"
}

if [[ "${CHR}" == "all" ]]; then
  while IFS= read -r c; do
    [[ -z "${c}" ]] && continue
    run_one "${c}"
  done < "${BEAGLE_DIR}/chr.list"
else
  run_one "${CHR}"
fi

qc_log "Q04 DONE. Figures in ${QC_FIGS}/"
