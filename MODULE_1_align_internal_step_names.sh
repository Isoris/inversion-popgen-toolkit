#!/usr/bin/env bash
# =============================================================================
# MODULE_1_align_internal_step_names.sh
#
# PROBLEM: MODULE_1 scripts were renamed on disk to the A/B-track convention
# (STEP_A02_..., SLURM_B01_..., etc.) but their INTERNAL $STEP= variables
# still carry the old S01..S16 names. As a result:
#   - Sidecar files (${STEP}.arg, ${STEP}.results) are written under the
#     OLD names, so users can't easily link sidecar -> script.
#   - Header comments reference legacy script names.
#   - Any grep for "STEP_A03" in sidecar contents fails.
#
# FIX: For each script, rewrite its $STEP assignment to match its filename
# (minus the extension). Also rewrite any header comment referencing the
# legacy S0X_ name.
#
# Idempotent: safe to re-run. Prints a summary of all edits.
#
# Usage:
#   bash MODULE_1_align_internal_step_names.sh [--dry-run] [--module-dir PATH]
#
# Defaults to the inversion-popgen-toolkit clone in cwd.
# =============================================================================
set -euo pipefail

DRY_RUN=0
MODULE_DIR="./inversion-popgen-toolkit/Modules/MODULE_1_read_prep"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run)    DRY_RUN=1; shift ;;
    --module-dir) MODULE_DIR="$2"; shift 2 ;;
    --help|-h)    sed -n '3,25p' "$0"; exit 0 ;;
    *) echo "[WARN] unknown arg: $1"; shift ;;
  esac
done

[[ -d "$MODULE_DIR" ]] || { echo "[ERR] not a directory: $MODULE_DIR"; exit 1; }

# The correspondence is: <new filename stem> → <old S0X_ name stem>.
# Pairs derived from the actual $STEP="S0X_..." assignments currently in each file.
declare -A MAPPING=(
  # steps/
  [STEP_A02_fastp_eof_redo_table]=S02_fastp_eof_redo_table
  [STEP_A03_extract_fastp_qc]=S03_extract_fastp_qc
  [STEP_A05_species_assign_kmers]=S05_species_assign_kmers
  [STEP_A06_plot_species_assign]=S06_plot_species_assign
  [STEP_A07_inventory_fastp_vs_bam]=S07_inventory_fastp_vs_bam
  [STEP_A08_check_bam_bai_pairs]=S08_check_bam_bai_pairs
  [STEP_B03_get_tlen_percentiles]=S11_get_tlen_percentiles
  [STEP_B04_get_insert_size_stats]=S12_get_insert_size_stats
  [STEP_B07_summarize_bam_qc]=S15_summarize_bam_qc
  [STEP_B08_make_bam_provenance]=S16_make_bam_provenance
  # slurm/
  [SLURM_A01_fastp_trim_reads]=S01_fastp_trim_reads
  [SLURM_A04_species_assign_mash]=S04_species_assign_mash
  [SLURM_B01_map_minimap2]=S09_map_minimap2
  [SLURM_B02_merge_markdup_clip]=S10_merge_markdup_clip
  [SLURM_B05_filter_bam_popgen]=S13_filter_bam_popgen
  [SLURM_B06_qc_depth_table]=S14_qc_depth_table
)

N_EDITS=0
N_FILES=0

apply_one() {
  local new_stem="$1"           # e.g. STEP_A02_fastp_eof_redo_table
  local old_stem="$2"           # e.g. S02_fastp_eof_redo_table
  local dir ext path
  case "$new_stem" in
    STEP_*) dir="${MODULE_DIR}/steps" ;;
    SLURM_*) dir="${MODULE_DIR}/slurm" ;;
  esac
  # Figure out extension on disk
  local found=""
  for ext in sh slurm py; do
    if [[ -f "${dir}/${new_stem}.${ext}" ]]; then
      found="${dir}/${new_stem}.${ext}"; break
    fi
  done
  if [[ -z "$found" ]]; then
    echo "  [SKIP] no file on disk for $new_stem"
    return
  fi
  path="$found"
  N_FILES=$((N_FILES+1))

  # Current internal STEP value, if any
  local cur
  cur=$(grep -oE 'STEP="[^"]+"' "$path" | head -1 | sed 's/STEP="\(.*\)"/\1/') || true

  if [[ -z "$cur" ]]; then
    echo "  [INFO] $path has no STEP= (ok for python); adding comment only"
  elif [[ "$cur" == "$new_stem" ]]; then
    echo "  [SKIP] $(basename "$path") already aligned"
    return
  fi

  # Dry-run report or apply
  if (( DRY_RUN )); then
    echo "  [DRY]  $(basename "$path") :  STEP=\"$cur\"  ->  STEP=\"$new_stem\""
  else
    # 1. Replace STEP="<old>" -> STEP="<new>"
    sed -i.bak -E "s|^(STEP=)\"${cur}\"|\\1\"${new_stem}\"|" "$path"
    # 2. Replace any standalone mention of the old stem in header comments
    #    (e.g. "# S02_fastp_eof_redo_table.sh"). Limit to lines starting with #.
    sed -i -E "\|^#|{ s|\\b${old_stem}\\b|${new_stem}|g }" "$path"
    # Remove .bak after confirming no corruption
    rm -f "${path}.bak"
    echo "  [OK]   $(basename "$path") :  STEP -> \"$new_stem\""
    N_EDITS=$((N_EDITS+1))
  fi
}

MODE_LABEL=$( (( DRY_RUN )) && echo "(dry run)" || echo "(apply)" )
echo "══════════════════════════════════════════════════════════════════════"
echo "  MODULE_1 internal STEP alignment $MODE_LABEL"
echo "  module_dir: $MODULE_DIR"
echo "══════════════════════════════════════════════════════════════════════"

for new_stem in "${!MAPPING[@]}"; do
  apply_one "$new_stem" "${MAPPING[$new_stem]}"
done

echo ""
echo "──────────────────────────────────────────────────────────────────────"
echo "  files inspected: $N_FILES"
echo "  files edited:    $N_EDITS"
if (( DRY_RUN )); then
  echo "  dry-run only: no files written. Drop --dry-run to apply."
fi
