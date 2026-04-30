#!/bin/bash
# =============================================================================
# sync_jsons.sh — copy precomp JSONs from cluster paths into Atlas/json/
# =============================================================================
# The cluster keeps verbose filenames (C_gar_LG28_phase2_ghsl.json) because
# the production pipeline cares about phase numbers and species. The atlas
# folder uses short names (LG28_ghsl.json) because you read them every day.
# This script bridges the two.
#
# Edit the SOURCE_* paths below for your layout. Then:
#
#   ./sync_jsons.sh                # copy all newer files
#   ./sync_jsons.sh --dry-run      # show what would be copied
#   ./sync_jsons.sh LG28           # only this chromosome
#
# Idempotent: only copies files whose source mtime is newer than dest.
# =============================================================================

set -euo pipefail

# ---- Source paths on the cluster (or wherever the JSONs live) ---------------
# Edit these for your layout.
SOURCE_DOSAGE_DIR="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_codebase_v8.5/phase_2_discovery/2a_local_pca/precomp_json"
SOURCE_GHSL_DIR="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_codebase_v8.5/phase_2_discovery/2e_ghsl_discovery/05_json_out"
SOURCE_THETA_DIR="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/phase_2_outputs/2f_theta_discovery/json_out"

# ---- Destination — assumed to be Atlas/json/ relative to this script --------
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DEST_DIR="${SCRIPT_DIR}/../json"
mkdir -p "$DEST_DIR"

# ---- Argument parsing -------------------------------------------------------
DRY_RUN=0
ONLY_CHROM=""
for a in "$@"; do
  case "$a" in
    --dry-run|-n) DRY_RUN=1 ;;
    LG[0-9][0-9]) ONLY_CHROM="$a" ;;
    -h|--help)
      sed -n '2,20p' "$0"; exit 0 ;;
    *) echo "[sync] unknown arg: $a" >&2; exit 2 ;;
  esac
done

# ---- Helper: copy if source newer than dest --------------------------------
copy_one() {
  local src="$1"
  local dst="$2"
  if [[ ! -f "$src" ]]; then
    return  # source missing, silently skip (other streams may exist)
  fi
  if [[ -f "$dst" && "$src" -ot "$dst" ]]; then
    return  # dest is newer or equal, skip
  fi
  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "[would-copy] $src"
    echo "          → $dst"
  else
    cp -p "$src" "$dst"
    echo "[copied] $(basename "$dst")"
  fi
}

# ---- Iterate chromosomes ---------------------------------------------------
chrom_list=()
if [[ -n "$ONLY_CHROM" ]]; then
  chrom_list=("$ONLY_CHROM")
else
  for n in $(seq -w 1 28); do
    chrom_list+=("LG${n}")
  done
fi

for short in "${chrom_list[@]}"; do
  full="C_gar_${short}"

  # Dosage: <full>_phase4a_dosage.json or similar; adjust as needed
  copy_one "${SOURCE_DOSAGE_DIR}/${full}/${full}_precomp.json" \
           "${DEST_DIR}/${short}_dosage.json"

  # GHSL
  copy_one "${SOURCE_GHSL_DIR}/${full}/${full}_phase2_ghsl.json" \
           "${DEST_DIR}/${short}_ghsl.json"

  # Theta
  copy_one "${SOURCE_THETA_DIR}/${full}/${full}_phase2_theta.json" \
           "${DEST_DIR}/${short}_theta.json"
done

# ---- Summary ---------------------------------------------------------------
n_files=$(find "$DEST_DIR" -maxdepth 1 -name '*.json' | wc -l)
echo ""
echo "[sync] $n_files JSON files in $DEST_DIR"
if [[ "$DRY_RUN" -eq 1 ]]; then
  echo "[sync] dry-run only; nothing was copied"
fi
