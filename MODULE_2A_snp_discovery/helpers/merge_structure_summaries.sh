#!/usr/bin/env bash
###############################################################################
# helpers/merge_structure_summaries.sh
#
# Merges best-seed-by-K tables from multiple structure runs (e.g. all + pruned)
# into one combined summary with a sample_set column.
#
# Usage:
#   bash helpers/merge_structure_summaries.sh
#
# Looks for:
#   ${THIN_DIR}/05_structure_all/canonical_exports/best_seed_by_K.tsv
#   ${THIN_DIR}/05_structure_pruned/canonical_exports/best_seed_by_K.tsv
#   (and any other 05_structure_*/canonical_exports/best_seed_by_K.tsv)
#
# Produces:
#   ${THIN_DIR}/combined_best_seed_by_K.tsv
###############################################################################
set -euo pipefail
source "$(dirname "$0")/../config.sh"

OUTFILE="${THIN_DIR}/combined_best_seed_by_K.tsv"
TMPFILE="${THIN_DIR}/.tmp_combined_best.tsv"

echo "[$(timestamp)] Merging structure summaries..."

first=1
for dir in "${THIN_DIR}"/05_structure_*/canonical_exports; do
  F="${dir}/best_seed_by_K.tsv"
  [[ -s "$F" ]] || continue

  # Extract sample_set from directory name (05_structure_all -> all, 05_structure_pruned -> pruned)
  SAMPLE_SET="$(basename "$(dirname "$dir")" | sed 's/^05_structure_//')"

  if (( first )); then
    # Header + sample_set column
    head -1 "$F" | awk -v OFS="\t" '{print "sample_set", $0}' > "$TMPFILE"
    first=0
  fi

  # Data rows with sample_set prepended
  tail -n +2 "$F" | awk -v ss="$SAMPLE_SET" -v OFS="\t" '{print ss, $0}' >> "$TMPFILE"

  echo "[INFO] Added sample_set=${SAMPLE_SET} ($(( $(wc -l < "$F") - 1 )) rows)"
done

if [[ -s "$TMPFILE" ]]; then
  mv "$TMPFILE" "$OUTFILE"
  echo "[DONE] Combined: $OUTFILE ($(( $(wc -l < "$OUTFILE") - 1 )) rows)"
else
  echo "[WARN] No best_seed_by_K.tsv files found"
  rm -f "$TMPFILE"
fi
