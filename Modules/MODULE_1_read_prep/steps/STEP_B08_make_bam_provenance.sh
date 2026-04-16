#!/usr/bin/env bash
###############################################################################
# S16_make_bam_provenance.sh
#
# Builds a provenance table linking each sample to its raw minimap2 BAM
# and final popgen-filtered BAM path. This table is used by downstream
# modules (MODULE_2, MODULE_3) as the definitive sample manifest.
#
# Usage:
#   bash S16_make_bam_provenance.sh
#
# Outputs:
#   ${OUT}
#   S16_make_bam_provenance.arg
#   S16_make_bam_provenance.results
###############################################################################
set -euo pipefail

STEP="S16_make_bam_provenance"
timestamp(){ date '+%F %T'; }

ROOT="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
LIST="${ROOT}/popstruct_thin/list_of_samples_one_per_line_same_bamfile_list.tsv"
OUT="${ROOT}/pa_roary_results/00_manifests/sample_bam_minimap2_vs_P99TLENMAPQ30.tsv"

mkdir -p "$(dirname "$OUT")"

ARGFILE="$(dirname "$OUT")/${STEP}.arg"
RESULTSFILE="$(dirname "$OUT")/${STEP}.results"

# ---- Write .arg ----
{
  echo -e "key\tvalue"
  echo -e "step\t${STEP}"
  echo -e "script_name\t$(basename "$0")"
  echo -e "datetime\t$(timestamp)"
  echo -e "host\t$(hostname)"
  echo -e "root\t${ROOT}"
  echo -e "sample_list\t${LIST}"
  echo -e "output\t${OUT}"
  echo -e "raw_bam_search\t${ROOT}/01-minimap2-bams"
  echo -e "filtered_bam_search\t${ROOT}/02-merged_per_sample/<sample>/*.filtered.bam"
} > "$ARGFILE"

n_found=0; n_missing=0

{
  echo -e "Sample\tBAM_minimap2\tBAM_filtered_P99TLEN"

  while IFS= read -r sample; do
    [[ -z "$sample" ]] && continue

    raw_bam="$(find "${ROOT}/01-minimap2-bams" \
      -type f -name "${sample}.*.sr.bam" \
      ! -name "*.tmp.*.bam" \
      ! -path "*/qc_dedup/*" \
      | sort | head -n 1)"

    filt_bam="$(find "${ROOT}/02-merged_per_sample/${sample}" \
      -type f -name "*.filtered.bam" \
      | sort | head -n 1)"

    if [[ -n "$raw_bam" && -n "$filt_bam" ]]; then
      n_found=$((n_found+1))
    else
      n_missing=$((n_missing+1))
    fi

    echo -e "${sample}\t${raw_bam:-NA}\t${filt_bam:-NA}"
  done < "$LIST"
} > "$OUT"

n_total=$(( $(wc -l < "$OUT") - 1 ))

# ---- Write .results ----
{
  echo -e "key\tvalue\tdescription"
  echo -e "arg_file\t${ARGFILE}\tParameter record"
  echo -e "provenance_table\t${OUT}\tSample-to-BAM provenance mapping"
  echo -e "n_samples\t${n_total}\tTotal samples in table"
  echo -e "n_both_found\t${n_found}\tSamples with both raw and filtered BAM"
  echo -e "n_any_missing\t${n_missing}\tSamples missing at least one BAM"
} > "$RESULTSFILE"

echo "[DONE] Wrote: $OUT ($n_total samples, $n_found complete, $n_missing missing)"
