#!/usr/bin/env bash
# =============================================================================
# smart_tree.sh — Compact directory listing with pattern collapsing
#
# Like tree but when 100 files share a pattern (e.g. CGA*.markdup.bam),
# shows one line: CGA{009..461}.markdup.bam (226 files, 1.6 TB)
#
# Usage: bash smart_tree.sh [directory] [--depth N] [--min-collapse 3]
# =============================================================================

set -euo pipefail

DIR="${1:-.}"
MAX_DEPTH=3
MIN_COLLAPSE=3

shift || true
while [[ $# -gt 0 ]]; do
  case "$1" in
    --depth) MAX_DEPTH="$2"; shift 2 ;;
    --min-collapse) MIN_COLLAPSE="$2"; shift 2 ;;
    *) shift ;;
  esac
done

echo "=== ${DIR} ==="
echo ""

find "$DIR" -maxdepth "$MAX_DEPTH" -type d | sort | while read -r d; do
  # Indent based on depth
  rel="${d#$DIR}"
  rel="${rel#/}"
  depth=$(echo "$rel" | tr -cd '/' | wc -c)
  indent=""
  for ((i=0; i<depth; i++)); do indent="  ${indent}"; done

  if [[ -z "$rel" ]]; then
    dirlab="."
  else
    dirlab="$rel/"
  fi

  # Count files and total size in this dir (non-recursive)
  n_files=$(find "$d" -maxdepth 1 -type f 2>/dev/null | wc -l)
  if [[ "$n_files" -eq 0 ]]; then
    n_subdirs=$(find "$d" -maxdepth 1 -type d 2>/dev/null | wc -l)
    n_subdirs=$((n_subdirs - 1))
    if [[ "$n_subdirs" -gt 0 ]]; then
      echo "${indent}📁 ${dirlab} (${n_subdirs} subdirs)"
    fi
    continue
  fi

  dir_size=$(du -sh "$d" --max-depth=0 2>/dev/null | cut -f1)
  echo "${indent}📁 ${dirlab} (${n_files} files, ${dir_size})"

  # Group files by extension + pattern
  find "$d" -maxdepth 1 -type f -printf '%f\t%s\n' 2>/dev/null | sort | \
  awk -F'\t' -v min_c="$MIN_COLLAPSE" -v indent="${indent}  " '
  {
    fname = $1; fsize = $2
    # Extract extension
    if (match(fname, /\.[^.]+$/)) {
      ext = substr(fname, RSTART)
    } else {
      ext = ""
    }
    # Extract prefix pattern: replace digits/sample-IDs with wildcard
    pattern = fname
    # Common patterns: CGA123 -> CGA*, LG01 -> LG*, chr1 -> chr*
    gsub(/[0-9]+/, "*", pattern)
    
    key = pattern
    counts[key]++
    sizes[key] += fsize
    if (!(key in first)) first[key] = fname
    last[key] = fname
    exts[key] = ext
  }
  END {
    # Sort keys
    n = asorti(counts, sorted)
    for (i = 1; i <= n; i++) {
      k = sorted[i]
      c = counts[k]
      sz = sizes[k]
      
      # Human-readable size
      if (sz > 1099511627776) s = sprintf("%.1f TB", sz/1099511627776)
      else if (sz > 1073741824) s = sprintf("%.1f GB", sz/1073741824)
      else if (sz > 1048576) s = sprintf("%.1f MB", sz/1048576)
      else if (sz > 1024) s = sprintf("%.1f KB", sz/1024)
      else s = sz " B"
      
      if (c >= min_c) {
        # Collapsed: show pattern
        printf "%s  %s (%d files, %s)\n", indent, k, c, s
      } else {
        # Show individual files
        printf "%s  %s (%s)\n", indent, first[k], s
        if (c > 1) printf "%s  %s (%s)\n", indent, last[k], "..."
      }
    }
  }'
  echo ""
done
