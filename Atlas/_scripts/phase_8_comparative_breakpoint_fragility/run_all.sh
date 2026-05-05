#!/usr/bin/env bash
# =============================================================================
# run_all.sh — end-to-end runner for phase_8_comparative_breakpoint_fragility
# =============================================================================
#
# Runs steps 00 -> 07 in order against the configs in ./config/.
#
# Required inputs (must exist before running):
#   config/species_manifest.tsv
#   config/candidate_breakpoints.tsv
#   config/te_file_manifest.tsv         (or pass --te_root and let step 00 build it)
#
# Optional:
#   config/synteny_mapping.tsv          (without it, focal-only run)
#   config/species_map.tsv              (for step 00, manual species overrides)
#
# Usage:
#   bash run_all.sh                                    # use existing manifest
#   bash run_all.sh --te_root /path/to/messy/te_dir    # rebuild manifest first
#
# Failures abort. All logs land in output/logs/. Run manifests live as
# .run.json sidecars next to every output.
# =============================================================================
set -euo pipefail

cd "$(dirname "$0")"

PY=${PY:-python3}
RSCRIPT=${RSCRIPT:-Rscript}

TE_ROOT=""
while (( "$#" )); do
  case "$1" in
    --te_root) TE_ROOT="$2"; shift 2;;
    *) echo "unknown arg: $1"; exit 2;;
  esac
done

mkdir -p output/{normalized_te,breakpoint_windows,density/focal,density/comparative,comparative/map,classification,per_candidate_json,per_chromosome_json,summary_tables,plots,logs}

echo "[run_all] step 00: TE manifest"
if [[ -n "$TE_ROOT" ]]; then
  $PY 00_manifest/build_te_manifest.py \
    --te_root  "$TE_ROOT" \
    --out      config/te_file_manifest.tsv \
    --species_map config/species_map.tsv 2>/dev/null \
    --log_dir  output/logs || \
  $PY 00_manifest/build_te_manifest.py \
    --te_root  "$TE_ROOT" \
    --out      config/te_file_manifest.tsv \
    --log_dir  output/logs
fi
test -s config/te_file_manifest.tsv || { echo "missing config/te_file_manifest.tsv"; exit 1; }

echo "[run_all] step 01: normalize TE annotations"
$PY 01_normalize_te/normalize_te_annotations.py \
  --manifest config/te_file_manifest.tsv \
  --species  config/species_manifest.tsv \
  --out_dir  output/normalized_te \
  --log_dir  output/logs

echo "[run_all] step 02: build breakpoint windows"
SYNTENY_ARGS=()
if [[ -s config/synteny_mapping.tsv ]]; then
  SYNTENY_ARGS=(--synteny config/synteny_mapping.tsv)
fi
$PY 02_breakpoint_windows/make_breakpoint_windows.py \
  --candidates config/candidate_breakpoints.tsv \
  --species_manifest config/species_manifest.tsv \
  "${SYNTENY_ARGS[@]}" \
  --out_dir output/breakpoint_windows \
  --log_dir output/logs

echo "[run_all] step 03: compute TE density (focal)"
FOCAL_ID=$(awk -F '\t' 'NR==1{for(i=1;i<=NF;i++){h[$i]=i}} NR>1 && tolower($h["is_focal"])=="true"{print $h["species_id"]; exit}' config/species_manifest.tsv)
FOCAL_FAI=$(awk -F '\t' 'NR==1{for(i=1;i<=NF;i++){h[$i]=i}} NR>1 && tolower($h["is_focal"])=="true"{print $h["fai_path"]; exit}' config/species_manifest.tsv)
test -n "$FOCAL_ID" || { echo "no focal species in species_manifest.tsv"; exit 1; }

$PY 03_density/compute_te_density.py \
  --windows    output/breakpoint_windows/${FOCAL_ID}.windows.bed \
  --te_bed     output/normalized_te/${FOCAL_ID}.te.bed.gz \
  --fai        "$FOCAL_FAI" \
  --species_id "$FOCAL_ID" \
  --is_focal \
  --out_dir    output/density/focal \
  --log_dir    output/logs

if [[ -s config/synteny_mapping.tsv ]]; then
  echo "[run_all] step 04a: map homologous breakpoints"
  $PY 04_comparative/map_homologous_breakpoints.py \
    --candidates config/candidate_breakpoints.tsv \
    --species_manifest config/species_manifest.tsv \
    --synteny  config/synteny_mapping.tsv \
    --out_dir  output/comparative/map \
    --log_dir  output/logs

  echo "[run_all] step 04b: compute comparative TE density"
  $PY 04_comparative/compute_comparative_te_density.py \
    --species_manifest config/species_manifest.tsv \
    --windows_dir output/breakpoint_windows \
    --te_dir      output/normalized_te \
    --out_dir     output/density/comparative \
    --log_dir     output/logs
else
  echo "[run_all] no synteny_mapping.tsv — skipping comparative steps"
fi

echo "[run_all] step 05: classify"
$PY 05_classification/classify_breakpoint_fragility.py \
  --candidates config/candidate_breakpoints.tsv \
  --species_manifest config/species_manifest.tsv \
  --homologs_dir output/comparative/map \
  --focal_density_dir output/density/focal \
  --comparative_density_dir output/density/comparative \
  --out_dir output/classification \
  --log_dir output/logs

echo "[run_all] step 06: export final JSON + summary TSV"
$PY 06_json_export/export_breakpoint_fragility_json.py \
  --candidates config/candidate_breakpoints.tsv \
  --species_manifest config/species_manifest.tsv \
  --focal_density_dir output/density/focal \
  --comparative_density_dir output/density/comparative \
  --homologs_dir output/comparative/map \
  --classification_dir output/classification \
  --out_root output \
  --log_dir output/logs

if command -v "$RSCRIPT" >/dev/null 2>&1; then
  echo "[run_all] step 07: plots"
  $RSCRIPT 07_plots/plot_breakpoint_te_density.R \
    --summary  output/summary_tables/breakpoint_fragility_summary.tsv \
    --json_dir output/per_candidate_json \
    --out_dir  output/plots || echo "[run_all] plots failed (non-fatal)"
else
  echo "[run_all] $RSCRIPT not available; skipping plots"
fi

echo "[run_all] DONE"
echo "  -> output/per_candidate_json/"
echo "  -> output/per_chromosome_json/"
echo "  -> output/summary_tables/breakpoint_fragility_summary.tsv"
echo "  -> output/plots/"
