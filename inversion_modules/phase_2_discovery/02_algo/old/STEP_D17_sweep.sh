#!/usr/bin/env bash
# =============================================================================
# STEP_D17_sweep.sh
#
# Runs STEP_D17_multipass.R + STEP_D17c_overlay_plot_v2.R across a small
# parameter grid for one chromosome. One PDF per setting, named with the
# parameter values so pages line up tile-by-tile across PDFs.
#
# Sweep dimensions (8 settings total):
#   post_trim_nms     : OFF, 0.20, 0.30, 0.50
#   trim_drop_z       : 0.5, 1.0
#
# Usage:
#   bash STEP_D17_sweep.sh \
#     --precomp <precomp.rds> --sim_mat <sim_mat.rds> \
#     --chr <label> --outroot <root_dir> \
#     [--rscript <path_to_Rscript>] \
#     [--multipass <path>] [--overlay <path>]
#
# Output layout:
#   <outroot>/<chr>/sweep_<setting_id>/
#       <chr>_d17_catalogue.tsv
#       <chr>_d17c_overlay.pdf
#   <outroot>/<chr>/_summary.tsv
# =============================================================================

set -euo pipefail

PRECOMP=""
SIM_MAT=""
CHR=""
OUTROOT=""
RSCRIPT="Rscript"
MULTIPASS_R="$(dirname "$0")/STEP_D17_multipass.R"
OVERLAY_R="$(dirname "$0")/STEP_D17c_overlay_plot_v2.R"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --precomp)   PRECOMP="$2"; shift 2 ;;
    --sim_mat)   SIM_MAT="$2"; shift 2 ;;
    --chr)       CHR="$2"; shift 2 ;;
    --outroot)   OUTROOT="$2"; shift 2 ;;
    --rscript)   RSCRIPT="$2"; shift 2 ;;
    --multipass) MULTIPASS_R="$2"; shift 2 ;;
    --overlay)   OVERLAY_R="$2"; shift 2 ;;
    *) echo "[sweep] unknown arg: $1" >&2; exit 1 ;;
  esac
done

if [[ -z "$PRECOMP" || -z "$SIM_MAT" || -z "$CHR" || -z "$OUTROOT" ]]; then
  echo "[sweep] usage: --precomp <f> --sim_mat <f> --chr <label> --outroot <dir>" >&2
  exit 2
fi
for f in "$PRECOMP" "$SIM_MAT" "$MULTIPASS_R" "$OVERLAY_R"; do
  [[ -f "$f" ]] || { echo "[sweep] missing: $f" >&2; exit 3; }
done

CHR_ROOT="${OUTROOT}/${CHR}"
mkdir -p "$CHR_ROOT"
SUMMARY="${CHR_ROOT}/_summary.tsv"
printf "setting_id\tl2_method\tcross_term\tscore_min\tcatalogue_rows\tphase_qc_rows\n" \
  > "$SUMMARY"

# Settings grid: compares scan vs refiner (z mean) vs refiner (zq, lower-tail).
# All use nested architecture + fragmented_emit ON, dz0.5 + post_nms 0.30.
# Format: l2_method  cross_term  score_min
SETTINGS=(
  "scan      -    -"
  "refiner   z    0.05"
  "refiner   zq   1.0"
  "refiner   zq   1.5"
)

for setting in "${SETTINGS[@]}"; do
  read -r L2M CROSS SMIN <<< "$setting"
  if [[ "$L2M" == "scan" ]]; then
    SID="l2scan"
  else
    SID="l2refiner_${CROSS}_smin${SMIN}"
  fi
  SETDIR="${CHR_ROOT}/sweep_${SID}"
  mkdir -p "$SETDIR"

  echo "================================================================="
  echo "[sweep] setting: $SID"
  echo "[sweep]   l2_method        = $L2M"
  if [[ "$L2M" == "refiner" ]]; then
    echo "[sweep]   cross_term       = $CROSS"
    echo "[sweep]   score_min        = $SMIN"
  fi
  echo "================================================================="

  # ---- run multipass ------------------------------------------------------
  MP_ARGS=(
    --precomp "$PRECOMP"
    --sim_mat "$SIM_MAT"
    --chr     "$CHR"
    --outdir  "$SETDIR"
    --architecture nested
    --l2_method    "$L2M"
    --trim_core TRUE
    --trim_drop_z 0.5
    --post_trim_nms TRUE
    --post_trim_nms_overlap 0.30
    --fragmented_emit TRUE
  )
  if [[ "$L2M" == "refiner" ]]; then
    MP_ARGS+=( --refiner_cross_term "$CROSS"
               --refiner_score_min  "$SMIN" )
  fi

  "$RSCRIPT" "$MULTIPASS_R" "${MP_ARGS[@]}" \
    2>&1 | tee "${SETDIR}/multipass.log"

  CAT_FILE="${SETDIR}/${CHR}_d17_catalogue.tsv"
  PQC_FILE="${SETDIR}/${CHR}_phase_qc_candidates.tsv"
  if [[ ! -f "$CAT_FILE" ]]; then
    echo "[sweep]   no catalogue produced; skipping overlay for $SID"
    printf "%s\t%s\t%s\t%s\tNA\tNA\n" \
      "$SID" "$L2M" "$CROSS" "$SMIN" >> "$SUMMARY"
    continue
  fi

  # rows = (lines - 1) for header
  CAT_ROWS=$(($(wc -l < "$CAT_FILE") - 1))
  PQC_ROWS="NA"
  if [[ -f "$PQC_FILE" ]]; then
    PQC_ROWS=$(($(wc -l < "$PQC_FILE") - 1))
  fi
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$SID" "$L2M" "$CROSS" "$SMIN" "$CAT_ROWS" "$PQC_ROWS" >> "$SUMMARY"

  # ---- run overlay --------------------------------------------------------
  "$RSCRIPT" "$OVERLAY_R" \
    --precomp   "$PRECOMP" \
    --sim_mat   "$SIM_MAT" \
    --catalogue "$CAT_FILE" \
    --chr       "$CHR" \
    --outdir    "$SETDIR" \
    --n_tiles   5 \
    --z_mode    diagonal \
    2>&1 | tee "${SETDIR}/overlay.log"

  # Rename the PDF so the filename carries the setting id (helps when
  # opening many at once / stacking pages).
  SRC_PDF="${SETDIR}/${CHR}_d17c_overlay.pdf"
  DST_PDF="${SETDIR}/${CHR}_d17c_overlay__${SID}.pdf"
  if [[ -f "$SRC_PDF" ]]; then
    cp "$SRC_PDF" "$DST_PDF"
  fi
done

echo ""
echo "[sweep] all done."
echo "[sweep] summary: $SUMMARY"
column -t -s $'\t' "$SUMMARY"
echo ""
echo "[sweep] PDFs:"
find "$CHR_ROOT" -name "*_d17c_overlay__*.pdf" | sort
