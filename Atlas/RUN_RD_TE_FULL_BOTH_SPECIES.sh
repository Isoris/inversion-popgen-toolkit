#!/usr/bin/env bash
# =============================================================================
# RUN_RD_TE_FULL_BOTH_SPECIES.sh — comprehensive TE density JSONs (login node)
# =============================================================================
#
# Runs STEP_RD_TE_density_full.R for both catfish species in sequence:
#   - Clarias gariepinus    (haplotype A of the F1 hybrid: fClaHyb_Gar_LG)
#   - Clarias macrocephalus (haplotype B of the F1 hybrid: fClaHyb_Mac_LG)
#
# Output structure:
#   /project/lt200308-agbsci/02-TE_catfish/EDTA2.2_results/
#     Cgar_TE_density/      ← 28 JSONs from C. gariepinus (~6-8 MB each)
#     Cmac_TE_density/      ← 27 JSONs from C. macrocephalus
#     RUN_REPORT_TE_FULL.txt
#
# Each JSON contains ~85 density arrays per chromosome:
#   - 26 TE classes × 3 stratifications (all/young/old) = 78 per-class arrays
#   - 5 aggregates (all_TE, young_TE_all, old_TE_all, insertion_count,
#     intact_element_count)
#   - 2 TSD layers (intact-only + intact+TIR-attribute-derived)
#
# Runtime: ~3-5 min per species on login node (TEanno.gff3 is 271-289 MB).
# Memory: ~2-3 GB peak. No SLURM allocation needed — login node is fine.
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Configuration — adjust if your layout shifts
# -----------------------------------------------------------------------------
EDTA_DIR="/project/lt200308-agbsci/02-TE_catfish/EDTA2.2_results"
OUT_BASE="${EDTA_DIR}"

# File stems — same convention as the TSD wrapper. The F1-hybrid haplotype
# names; deliberately NOT C_macro.fa.mod (the older independent macrocephalus
# assembly), which is not coordinate-aligned with the F1 hybrid haplotype.
GAR_STEM="fClaHyb_Gar_LG.fa.mod"
MAC_STEM="fClaHyb_Mac_LG.fa.mod"

# Window parameters (must match the Inversion Atlas's scrubber grid)
WINDOW_BP=5000
STEP_BP=5000

# Young-TE identity threshold (≥ this → young, recently inserted/active)
IDENTITY_YOUNG_THRESHOLD=0.95

RSCRIPT="/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="${SCRIPT_DIR}/STEP_RD_TE_density_full.R"

REPORT="${OUT_BASE}/RUN_REPORT_TE_FULL.txt"

# -----------------------------------------------------------------------------
# Pre-flight
# -----------------------------------------------------------------------------
[[ -d "${EDTA_DIR}" ]]  || { echo "ERROR: EDTA_DIR not found: ${EDTA_DIR}"; exit 1; }
[[ -f "${SCRIPT}" ]]    || { echo "ERROR: R script not found: ${SCRIPT}"; exit 1; }
[[ -x "${RSCRIPT}" ]]   || { echo "ERROR: Rscript not at ${RSCRIPT}"; exit 1; }

# -----------------------------------------------------------------------------
# Helper: run for one species
# Args: species_label  stem  out_subdir
# -----------------------------------------------------------------------------
run_one_species() {
  local species_label="$1"
  local stem="$2"
  local out_subdir="$3"

  local teanno_gff3="${EDTA_DIR}/${stem}.EDTA.TEanno.gff3"
  local intact_gff3="${EDTA_DIR}/${stem}.EDTA.intact.gff3"

  # .fai discovery — same fallback chain as the TSD wrapper
  local fai="${EDTA_DIR}/${stem}.fai"
  if [[ ! -f "${fai}" ]]; then
    for alt in "${EDTA_DIR}/${stem%.mod}.fai" \
               "${EDTA_DIR}/${stem}.fa.fai"   \
               "${EDTA_DIR}/${stem%.mod}.fa.fai"; do
      if [[ -f "${alt}" ]]; then fai="${alt}"; break; fi
    done
  fi

  echo
  echo "═══════════════════════════════════════════════════════════════"
  echo "  ${species_label}"
  echo "═══════════════════════════════════════════════════════════════"
  echo "  teanno_gff3 : ${teanno_gff3}  ($(du -h "${teanno_gff3}" 2>/dev/null | awk '{print $1}'))"
  echo "  intact_gff3 : ${intact_gff3}  ($(du -h "${intact_gff3}" 2>/dev/null | awk '{print $1}'))"
  echo "  fai         : ${fai}"
  echo "  out_dir     : ${out_subdir}"

  if [[ ! -f "${teanno_gff3}" ]]; then
    echo "  ✗ TEanno.gff3 not found — SKIPPING ${species_label}"
    return 1
  fi
  if [[ ! -f "${intact_gff3}" ]]; then
    echo "  ✗ intact.gff3 not found — SKIPPING ${species_label}"
    return 1
  fi
  if [[ ! -f "${fai}" ]]; then
    echo "  ✗ .fai not found — SKIPPING ${species_label}"
    echo "    (searched: ${stem}.fai, ${stem%.mod}.fai, ${stem}.fa.fai, ${stem%.mod}.fa.fai)"
    return 1
  fi

  mkdir -p "${out_subdir}"

  local t0=$(date +%s)
  "${RSCRIPT}" "${SCRIPT}" \
    --teanno_gff3 "${teanno_gff3}" \
    --intact_gff3 "${intact_gff3}" \
    --fai         "${fai}" \
    --window_bp   "${WINDOW_BP}" \
    --step_bp     "${STEP_BP}" \
    --out_dir     "${out_subdir}" \
    --species     "${species_label}" \
    --identity_young_threshold "${IDENTITY_YOUNG_THRESHOLD}"
  local t1=$(date +%s)

  local n_json=$(ls "${out_subdir}"/*_repeat_density_TEfull.json 2>/dev/null | wc -l)
  local total_size=$(du -sh "${out_subdir}" 2>/dev/null | awk '{print $1}')
  echo "  ✓ ${species_label}: ${n_json} JSONs (${total_size}) in $((t1 - t0))s"

  # Append to report
  {
    echo
    echo "── ${species_label} ──"
    echo "teanno_gff3 : ${teanno_gff3}"
    echo "intact_gff3 : ${intact_gff3}"
    echo "fai         : ${fai}"
    echo "out_dir     : ${out_subdir}"
    echo "n_jsons     : ${n_json}"
    echo "total_size  : ${total_size}"
    echo "wall_time_s : $((t1 - t0))"
    echo "  per-chrom n_classes (from JSON):"
    for f in "${out_subdir}"/*_repeat_density_TEfull.json; do
      [[ -f "$f" ]] || continue
      bn=$(basename "$f" _repeat_density_TEfull.json)
      n_cls=$(grep -oE '"n_classes":[0-9]+' "$f" | head -1 | grep -oE '[0-9]+')
      n_win=$(grep -oE '"n_windows":[0-9]+' "$f" | head -1 | grep -oE '[0-9]+')
      sz_kb=$(du -k "$f" | awk '{print $1}')
      printf "    %-15s n_windows=%-7s n_classes=%-4s size=%sKB\n" \
             "$bn" "${n_win:-?}" "${n_cls:-?}" "${sz_kb:-?}"
    done
  } >> "${REPORT}"

  return 0
}

# -----------------------------------------------------------------------------
# Initialize report
# -----------------------------------------------------------------------------
{
  echo "================================================================"
  echo "  Comprehensive TE density layer — RUN REPORT"
  echo "  $(date)"
  echo "  $(hostname)"
  echo "================================================================"
  echo "EDTA_DIR    : ${EDTA_DIR}"
  echo "OUT_BASE    : ${OUT_BASE}"
  echo "Gar stem    : ${GAR_STEM}"
  echo "Mac stem    : ${MAC_STEM}"
  echo "window_bp   : ${WINDOW_BP}"
  echo "step_bp     : ${STEP_BP}"
  echo "young_thr   : ${IDENTITY_YOUNG_THRESHOLD}"
  echo "R script    : ${SCRIPT}"
  echo "Rscript bin : ${RSCRIPT}"
} > "${REPORT}"

# -----------------------------------------------------------------------------
# Run both species
# -----------------------------------------------------------------------------
GAR_OK=0; MAC_OK=0

run_one_species "Clarias gariepinus"    "${GAR_STEM}" "${OUT_BASE}/Cgar_TE_density" && GAR_OK=1 || true
run_one_species "Clarias macrocephalus" "${MAC_STEM}" "${OUT_BASE}/Cmac_TE_density" && MAC_OK=1 || true

# -----------------------------------------------------------------------------
# Final summary
# -----------------------------------------------------------------------------
echo
echo "═══════════════════════════════════════════════════════════════"
echo "  Summary"
echo "═══════════════════════════════════════════════════════════════"
echo "  C. gariepinus    : $([[ ${GAR_OK} -eq 1 ]] && echo '✓ ok' || echo '✗ failed')"
echo "  C. macrocephalus : $([[ ${MAC_OK} -eq 1 ]] && echo '✓ ok' || echo '✗ failed')"
echo "  Report          : ${REPORT}"
echo

# Tail of report for quick look
echo "──── tail of report ────"
tail -50 "${REPORT}"
