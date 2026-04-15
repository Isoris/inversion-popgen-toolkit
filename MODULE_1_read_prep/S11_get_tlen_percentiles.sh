#!/usr/bin/env bash
###############################################################################
# S11_get_tlen_percentiles.sh
#
# Computes approximate p95 and p99 absolute TLEN values from a properly
# paired, primary, high-quality BAM.
#
# Usage:
#   bash S11_get_tlen_percentiles.sh sample.markdup.clip.bam
#   bash S11_get_tlen_percentiles.sh sample.markdup.clip.bam > tlen_pcts.tsv
#
# Filters: -f 0x2 (proper pair) -F 0xF0C (exclude unmapped/secondary/dup/supp) -q 60
#
# Output (stdout): tab-separated p95 and p99 values
# Sidecar: <bam_dir>/S11_get_tlen_percentiles.<sample>.arg
###############################################################################
set -euo pipefail

STEP="S11_get_tlen_percentiles"
timestamp(){ date '+%F %T'; }

BAM="${1:?Usage: $0 sample.bam}"
[[ -s "$BAM" ]] || { echo "ERROR: BAM not found: $BAM" >&2; exit 1; }

SAMPLE="$(basename "$BAM" | sed 's/\.merged.*//; s/\.bam$//')"
SDIR="$(dirname "$BAM")"
ARGFILE="${SDIR}/${STEP}.${SAMPLE}.arg"

MINMAPQ=60
REQUIRE_FLAGS=0x2     # proper pair
EXCLUDE_FLAGS=0xF0C   # unmapped, mate-unmapped, secondary, QCfail, dup, supp

{
  echo -e "key\tvalue"
  echo -e "step\t${STEP}"
  echo -e "script_name\t$(basename "$0")"
  echo -e "datetime\t$(timestamp)"
  echo -e "host\t$(hostname)"
  echo -e "bam\t${BAM}"
  echo -e "sample\t${SAMPLE}"
  echo -e "min_mapq\t${MINMAPQ}"
  echo -e "require_flags\t${REQUIRE_FLAGS}"
  echo -e "exclude_flags\t${EXCLUDE_FLAGS}"
} > "$ARGFILE"

samtools view -f "$REQUIRE_FLAGS" -F "$EXCLUDE_FLAGS" -q "$MINMAPQ" "$BAM" \
| awk '$9 != 0 { t = $9; if (t < 0) t = -t; print t }' \
| sort -n \
| awk -v sample="$SAMPLE" '
  { a[NR] = $1 }
  END {
    n = NR
    if (n == 0) { print "ERROR: no TLEN values" > "/dev/stderr"; exit 1 }
    p95 = int(n * 0.95); if (p95 < 1) p95 = 1
    p99 = int(n * 0.99); if (p99 < 1) p99 = 1
    printf "sample\tp95\tp99\tn_reads\n"
    printf "%s\t%d\t%d\t%d\n", sample, a[p95], a[p99], n
  }'
