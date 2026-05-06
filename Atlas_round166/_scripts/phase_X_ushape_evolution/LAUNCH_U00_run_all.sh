#!/usr/bin/env bash
# =============================================================================
# LAUNCH_U00_run_all.sh — offline fallback runner for U-shape evolution
# =============================================================================
# Primary entry point is the popstats server (POST /api/ushape/candidate).
# This script is the *fallback* — it runs the same logic as a one-shot batch
# against per-chromosome dosage TSVs.
#
# Usage:
#   LAUNCH_U00_run_all.sh \
#       --candidates  config/candidates.tsv \
#       --groups      config/candidate_groups.tsv \
#       --dosage_dir  /path/to/dosage/ \
#       --out_dir     output/ \
#       [--config     config/config_ushape.yaml] \
#       [--gw_windows /path/to/genomewide_window_stats]
#
# Required input shapes:
#   candidates.tsv         candidate_id, chrom, start_bp, end_bp,
#                          left_breakpoint, right_breakpoint
#   candidate_groups.tsv   candidate_id, sample_id, karyotype_group
#                          (karyotype_group ∈ {HOMO_1, HET, HOMO_2})
#   <chrom>.dosage.tsv.gz  chrom, pos, sample1, sample2, ... (0/1/2/NA)
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"

# defaults
CFG="config/config_ushape.yaml"
CAND=""; GROUPS=""; DOSAGE=""; OUT=""

while (( "$#" )); do
  case "$1" in
    --candidates) CAND="$2"; shift 2;;
    --groups)     GROUPS="$2"; shift 2;;
    --dosage_dir) DOSAGE="$2"; shift 2;;
    --out_dir)    OUT="$2"; shift 2;;
    --config)     CFG="$2"; shift 2;;
    *) echo "unknown arg: $1"; exit 2;;
  esac
done

for v in CAND GROUPS DOSAGE OUT; do
  [[ -n "${!v}" ]] || { echo "missing --${v,,}"; exit 2; }
done

mkdir -p "$OUT/window_stats" "$OUT/plots" "$OUT/logs"

# minimal yq-free YAML reader: grep numeric scalars
get() { awk -v k="$1:" '$1==k{print $2}' "$CFG"; }

N_WIN_INSIDE=$(get n_windows_inside)
MAX_FLANK=$(get max_flank_bp)
MIN_FLANK=$(get min_flank_bp)
MIN_WINDOW=$(get min_window_bp)
EDGE_FRAC=$(get edge_fraction)
CENTER_FRAC=$(get center_fraction)
MIN_GROUP_N=$(get min_group_n)

USC=$(get u_score_high)
IFL=$(get inside_flank_high)
IPK=$(get internal_peak_high)
ASY=$(get asymmetry_log2_high)
FEN=$(get fst_enrichment_high)
OLD=$(get oldness_min)
DXY=$(get dxy_min_inside)
FLT=$(get flatness_max_for_flat)
MWS=$(get min_inside_windows_for_shape)

NRAND=$(get n_random_background)

echo "[U00] step U01 — window stats"
Rscript cli/STEP_U01_compute_window_stats.R \
  --candidate_table "$CAND" --group_table "$GROUPS" \
  --dosage_dir "$DOSAGE" --out_dir "$OUT" \
  --n_windows_inside "$N_WIN_INSIDE" \
  --max_flank_bp "$MAX_FLANK" --min_flank_bp "$MIN_FLANK" \
  --min_window_bp "$MIN_WINDOW" \
  --edge_fraction "$EDGE_FRAC" --center_fraction "$CENTER_FRAC" \
  --min_group_n "$MIN_GROUP_N" \
  2>&1 | tee "$OUT/logs/U01.log"

echo "[U00] step U02 — score + classify"
Rscript cli/STEP_U02_score_and_classify_shapes.R \
  --raw_summary      "$OUT/candidate_raw_summary.tsv" \
  --window_stats_dir "$OUT/window_stats" --out_dir "$OUT" \
  --u_score_high "$USC" --inside_flank_high "$IFL" \
  --internal_peak_high "$IPK" --asymmetry_log2_high "$ASY" \
  --fst_enrichment_high "$FEN" --oldness_min "$OLD" \
  --dxy_min_inside "$DXY" --min_inside_windows_for_shape "$MWS" \
  --flatness_max_for_flat "$FLT" \
  2>&1 | tee "$OUT/logs/U02.log"

echo "[U00] step U03 — flank-resample matched background"
Rscript cli/STEP_U03_matched_background_scores.R \
  --candidate_table "$CAND" \
  --window_stats_dir "$OUT/window_stats" \
  --out_dir "$OUT" --n_random "$NRAND" \
  2>&1 | tee "$OUT/logs/U03.log"

echo "[U00] step U04 — clustering"
Rscript cli/STEP_U04_cluster_shape_profiles.R \
  --shape_scores "$OUT/candidate_shape_scores.tsv" \
  --classes      "$OUT/candidate_shape_classes.tsv" \
  --out_dir      "$OUT" \
  2>&1 | tee "$OUT/logs/U04.log"

echo "[U00] step U05 — PERMANOVA / ANOSIM validation"
Rscript cli/STEP_U05_validate_shape_classes.R \
  --feature_matrix "$OUT/shape_feature_matrix.tsv" \
  --class_table    "$OUT/candidate_shape_classes.tsv" \
  --cluster_table  "$OUT/shape_cluster_assignments.tsv" \
  --out_dir        "$OUT" \
  2>&1 | tee "$OUT/logs/U05.log"

echo "[U00] step U06 — JSON export"
Rscript cli/STEP_U06_export_ushape_json.R \
  --candidate_table     "$CAND" \
  --window_stats_dir    "$OUT/window_stats" \
  --raw_summary         "$OUT/candidate_raw_summary.tsv" \
  --shape_scores        "$OUT/candidate_shape_scores.tsv" \
  --shape_classes       "$OUT/candidate_shape_classes.tsv" \
  --cluster_assignments "$OUT/shape_cluster_assignments.tsv" \
  --out_json            "$OUT/ushape_evolution_v1.json" \
  --n_windows_inside "$N_WIN_INSIDE" \
  --max_flank_bp "$MAX_FLANK" --min_flank_bp "$MIN_FLANK" \
  --min_window_bp "$MIN_WINDOW" \
  --edge_fraction "$EDGE_FRAC" --center_fraction "$CENTER_FRAC" \
  --u_score_high "$USC" --inside_flank_high "$IFL" \
  --internal_peak_high "$IPK" --asymmetry_log2_high "$ASY" \
  --fst_enrichment_high "$FEN" --oldness_min "$OLD" \
  --dxy_min_inside "$DXY" --flatness_max_for_flat "$FLT" \
  2>&1 | tee "$OUT/logs/U06.log"

echo "[U00] step U07 — plots"
Rscript cli/STEP_U07_plot_candidate_profiles.R \
  --ushape_json "$OUT/ushape_evolution_v1.json" \
  --out_dir     "$OUT/plots" \
  2>&1 | tee "$OUT/logs/U07.log"

echo "[U00] DONE"
echo "  -> $OUT/ushape_evolution_v1.json"
echo "  -> $OUT/plots/"
