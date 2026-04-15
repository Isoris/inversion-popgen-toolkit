#!/usr/bin/env bash
###############################################################################
# S12_get_insert_size_stats.sh
#
# Computes detailed insert size statistics from a merged/markdup BAM:
# raw stats, robust MAD-trimmed stats, layout detection.
#
# Uses properly paired, primary, non-duplicate, non-supplementary reads
# with MAPQ >= 30 and same-reference mate pairs only.
#
# Usage:
#   bash S12_get_insert_size_stats.sh sample.markdup.bam > insert_stats.tsv
#
# Notes:
#   - PE151 means read_length=151, not insert size
#   - Standard Illumina PE layout is inward-facing (FR / ><)
#   - Trimming window: median ± 10*MAD (robust outlier removal)
#
# Sidecar: <bam_dir>/S12_get_insert_size_stats.<sample>.arg
###############################################################################
set -euo pipefail

STEP="S12_get_insert_size_stats"
timestamp(){ date '+%F %T'; }

[[ $# -eq 1 ]] || { echo "Usage: bash $0 sample.markdup.bam" >&2; exit 1; }
BAM="$1"
[[ -f "$BAM" ]] || { echo "ERROR: BAM not found: $BAM" >&2; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found" >&2; exit 1; }

SAMPLE="$(basename "$BAM" | sed 's/\.merged.*//; s/\.bam$//')"
SDIR="$(dirname "$BAM")"
ARGFILE="${SDIR}/${STEP}.${SAMPLE}.arg"

# Filters:
# -f 2      proper pair
# -F 3844   exclude unmapped(4), mate-unmapped(8), secondary(256), QCfail(512), dup(1024), supp(2048)
# -q 30     MAPQ >= 30
MINMAPQ=30
REQUIRE_FLAGS=2
EXCLUDE_FLAGS=3844

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
  echo -e "trimming_method\tmedian_pm_10MAD"
  echo -e "same_chr_only\tyes"
} > "$ARGFILE"

python3 - "$BAM" <<'PY'
import sys, subprocess, statistics

bam = sys.argv[1]

cmd = ["samtools", "view", "-f", "2", "-F", "3844", "-q", "30", bam]
vals = []

proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
for line in proc.stdout:
    if not line.strip():
        continue
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 9:
        continue
    # Keep only same-reference mate pairs
    if fields[6] != "=":
        continue
    try:
        tlen = abs(int(fields[8]))
    except ValueError:
        continue
    if tlen == 0:
        continue
    vals.append(tlen)

stderr = proc.stderr.read()
ret = proc.wait()
if ret != 0:
    sys.stderr.write(stderr)
    sys.exit(ret)

header = ("bam\tN_raw\tmean_raw\tsd_raw\tmedian_raw\tmad_raw\t"
          "trim_low\ttrim_high\tN_trimmed\tmean_trimmed\tsd_trimmed\t"
          "layout\tread_length")
print(header)

if len(vals) == 0:
    print(f"{bam}\t0\tNA\tNA\tNA\tNA\tNA\tNA\t0\tNA\tNA\t><\t151")
    sys.exit(0)

n_raw = len(vals)
mean_raw = sum(vals) / n_raw
sd_raw = (sum((x - mean_raw) ** 2 for x in vals) / (n_raw - 1)) ** 0.5 if n_raw > 1 else 0.0
median_raw = statistics.median(vals)
mad_raw = statistics.median([abs(x - median_raw) for x in vals])

# Robust trimming: median ± 10*MAD
trim_low = max(0, median_raw - 10 * mad_raw)
trim_high = median_raw + 10 * mad_raw
trimmed = [x for x in vals if trim_low <= x <= trim_high]

n_trim = len(trimmed)
mean_trim = sum(trimmed) / n_trim if n_trim > 0 else float("nan")
sd_trim = (
    (sum((x - mean_trim) ** 2 for x in trimmed) / (n_trim - 1)) ** 0.5
    if n_trim > 1 else 0.0
)

print(
    f"{bam}\t{n_raw}\t{mean_raw:.2f}\t{sd_raw:.2f}\t{median_raw:.2f}\t{mad_raw:.2f}\t"
    f"{trim_low:.2f}\t{trim_high:.2f}\t{n_trim}\t{mean_trim:.2f}\t{sd_trim:.2f}\t"
    f"><\t151"
)
PY
