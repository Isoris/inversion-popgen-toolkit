#!/bin/bash
# =============================================================================
# STEP_Q07_popstats.sh
# =============================================================================
# Run Engine F (region_popstats C binary) to compute cohort-wide and
# between-group population statistics at pestPG scale (default 50kb/10kb).
#
# Two grouping strategies produced in one call (region_popstats supports
# multiple groups natively):
#
#   Strategy A — ancestry-group (global): assigns each sample to its
#                dominant genome-wide ancestry from Engine B's local_Q.
#                Used chromosome-wide for demographic/substructure contrast.
#
#   Strategy B — inversion-genotype (region-local): k-means(3) clustering on
#                local PC1 loadings inside the shelf region (or any user
#                --group_region). Produces Hom1/Het/Hom2 group labels.
#                Used to compute the classical inversion signature
#                (Hom1-vs-Hom2 Fst spike at shelf, reduced theta_pi in
#                homozygote groups).
#
# Outputs (long-format, one row per window):
#   ${QC_TRACKS}/popstats_ancestry.<CHR>.tsv    (Strategy A)
#   ${QC_TRACKS}/popstats_invgt.<CHR>.tsv       (Strategy B)
#   ${QC_TRACKS}/invgt_assignments.<CHR>.tsv    (sample -> Hom1/Het/Hom2)
#
# Usage:
#   bash STEP_Q07_popstats.sh <CHR> [--shelf_start_mb N --shelf_end_mb M]
#   SHELF_START_MB=15 SHELF_END_MB=18 bash STEP_Q07_popstats.sh C_gar_LG28
# =============================================================================
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/00_config.sh"

CHR="${1:-}"
[[ -z "${CHR}" ]] && qc_die "Usage: $0 <CHR|all>"
: "${SHELF_START_MB:=}"
: "${SHELF_END_MB:=}"
: "${POPSTATS_WIN_BP:=50000}"
: "${POPSTATS_STEP_BP:=10000}"
: "${MIN_GROUP_N:=10}"
: "${INV_K:=3}"

# ---- Check Engine F binary --------------------------------------------------
if [[ ! -x "${REGION_POPSTATS_BIN}" ]]; then
  qc_log "ERROR: Engine F binary not found or not executable: ${REGION_POPSTATS_BIN}"
  qc_log "  Compile with: cd ${UNIFIED_ANCESTRY_DIR}/engines/fst_dxy && make"
  qc_log "  Or run: bash install.sh"
  exit 1
fi

run_one() {
  local chr="$1"
  qc_log "Q07 ${chr}: popstats via Engine F"

  local beagle="${BEAGLE_DIR}/main_qcpass.${chr}.beagle.gz"
  [[ -f "${beagle}" ]] || { qc_log "SKIP ${chr}: no beagle"; return 0; }
  local precomp="${PRECOMP_DIR}/${chr}.precomp.rds"
  [[ -f "${precomp}" ]] || precomp="${PRECOMP_DIR}/main_qcpass.${chr}.precomp.rds"
  [[ -f "${precomp}" ]] || { qc_log "SKIP ${chr}: no precomp"; return 0; }
  # Resolve local_Q_samples across layouts: flat, scale_1x, scale_dense,
  # scale_thin, and the v3 multi-K layout. scale_dense is preferred (higher
  # SNP resolution) for invgt group construction.
  local localq_samples=""
  for cand in \
    "${LOCAL_Q_DIR}/scale_dense/${chr}.local_Q_samples.tsv.gz" \
    "${LOCAL_Q_DIR}/scale_thin/${chr}.local_Q_samples.tsv.gz" \
    "${LOCAL_Q_DIR}/scale_1x/${chr}.local_Q_samples.tsv.gz" \
    "${LOCAL_Q_DIR}/${chr}.local_Q_samples.tsv.gz" \
    "${LOCAL_Q_DIR}/scale_dense/K08/${chr}.local_Q_samples.tsv.gz" \
    "${LOCAL_Q_DIR}/scale_thin/K08/${chr}.local_Q_samples.tsv.gz"; do
    if [[ -f "${cand}" ]]; then
      localq_samples="${cand}"
      break
    fi
  done
  [[ -n "${localq_samples}" ]] || qc_log "WARN ${chr}: no local_Q_samples found, Q07 will use MDS1-only grouping"

  # ---- Build group sample lists via R helper ------------------------------
  local groups_dir="${QC_OUT}/popstats_groups/${chr}"
  mkdir -p "${groups_dir}"
  local inv_assign_out="${QC_TRACKS}/invgt_assignments.${chr}.tsv"

  local r_args=(
    --precomp  "${precomp}"
    --chrom    "${chr}"
    --sample_list "${SAMPLE_LIST_POPSTATS}"
    --groups_dir  "${groups_dir}"
    --invgt_out   "${inv_assign_out}"
    --min_group_n "${MIN_GROUP_N}"
    --inv_k       "${INV_K}"
  )
  [[ -f "${localq_samples}" ]] && r_args+=( --localq_samples "${localq_samples}" )
  [[ -n "${SHELF_START_MB}" ]] && r_args+=( --shelf_start_mb "${SHELF_START_MB}" )
  [[ -n "${SHELF_END_MB}"   ]] && r_args+=( --shelf_end_mb   "${SHELF_END_MB}"   )

  ${RSCRIPT_BIN} --vanilla "${here}/R/q07_build_groups.R" "${r_args[@]}" \
    2> "${QC_LOGS}/q07_groups_${chr}.log"

  # ---- Strategy A: ancestry groups ----------------------------------------
  local anc_out="${QC_TRACKS}/popstats_ancestry.${chr}.tsv"
  local anc_spec=""
  # Build "name1:path1,name2:path2,..." spec from files in groups_dir/ancestry/
  if [[ -d "${groups_dir}/ancestry" ]]; then
    for f in "${groups_dir}/ancestry"/*.txt; do
      [[ -f "${f}" ]] || continue
      local n
      n=$(wc -l < "${f}")
      [[ "${n}" -lt "${MIN_GROUP_N}" ]] && { qc_log "  ancestry: skip $(basename "${f}" .txt) (n=${n} < ${MIN_GROUP_N})"; continue; }
      local name
      name=$(basename "${f}" .txt)
      if [[ -z "${anc_spec}" ]]; then
        anc_spec="${name}:${f}"
      else
        anc_spec="${anc_spec},${name}:${f}"
      fi
    done
  fi

  if [[ -n "${anc_spec}" ]]; then
    qc_log "Q07 ${chr}: ancestry groups=${anc_spec//,/; }"
    ${REGION_POPSTATS_BIN} \
      --beagle "${beagle}" \
      --sample_list "${SAMPLE_LIST_POPSTATS}" \
      --groups "${anc_spec}" \
      --fixed_win "${POPSTATS_WIN_BP}:${POPSTATS_STEP_BP}" \
      --chr "${chr}" \
      --out "${anc_out}" \
      2> "${QC_LOGS}/q07_popstats_ancestry_${chr}.log" || \
      qc_log "  ancestry popstats returned non-zero (see log)"
    if [[ -f "${anc_out}" ]]; then
      local nrows
      nrows=$(( $(wc -l < "${anc_out}") - 1 ))
      qc_log "  -> ${anc_out} (${nrows} rows)"
    fi
  else
    qc_log "Q07 ${chr}: no ancestry groups with n >= ${MIN_GROUP_N}, skipping Strategy A"
  fi

  # ---- Strategy B: inversion-genotype groups (from PCA) -------------------
  local invgt_out="${QC_TRACKS}/popstats_invgt.${chr}.tsv"
  local invgt_spec=""
  if [[ -d "${groups_dir}/invgt" ]]; then
    for f in "${groups_dir}/invgt"/*.txt; do
      [[ -f "${f}" ]] || continue
      local n
      n=$(wc -l < "${f}")
      [[ "${n}" -lt 3 ]] && { qc_log "  invgt: skip $(basename "${f}" .txt) (n=${n} < 3)"; continue; }
      local name
      name=$(basename "${f}" .txt)
      if [[ -z "${invgt_spec}" ]]; then
        invgt_spec="${name}:${f}"
      else
        invgt_spec="${invgt_spec},${name}:${f}"
      fi
    done
  fi

  if [[ -n "${invgt_spec}" ]]; then
    qc_log "Q07 ${chr}: invgt groups=${invgt_spec//,/; }"
    ${REGION_POPSTATS_BIN} \
      --beagle "${beagle}" \
      --sample_list "${SAMPLE_LIST_POPSTATS}" \
      --groups "${invgt_spec}" \
      --fixed_win "${POPSTATS_WIN_BP}:${POPSTATS_STEP_BP}" \
      --chr "${chr}" \
      --out "${invgt_out}" \
      2> "${QC_LOGS}/q07_popstats_invgt_${chr}.log" || \
      qc_log "  invgt popstats returned non-zero (see log)"
    if [[ -f "${invgt_out}" ]]; then
      local nrows
      nrows=$(( $(wc -l < "${invgt_out}") - 1 ))
      qc_log "  -> ${invgt_out} (${nrows} rows)"
    fi
  else
    qc_log "Q07 ${chr}: no invgt groups built (need --shelf_start_mb and --shelf_end_mb)"
  fi
}

if [[ "${CHR}" == "all" ]]; then
  while IFS= read -r c; do
    [[ -z "${c}" ]] && continue
    run_one "${c}"
  done < "${BEAGLE_DIR}/chr.list"
else
  run_one "${CHR}"
fi

qc_log "Q07 DONE."
