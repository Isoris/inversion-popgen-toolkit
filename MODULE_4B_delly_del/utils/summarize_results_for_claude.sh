#!/usr/bin/env bash
# summarize_delly_results_for_claude.sh
# Make a compact text summary of DEL pipeline result files for sharing with Claude

set -euo pipefail

OUT="${1:-claude_delly_results_summary.txt}"

files=(
  "07_final_catalogs/catalog_226.DEL.vcf.gz"
  "07_final_catalogs/catalog_226.DEL.bed"
  "07_final_catalogs/catalog_226.DEL.GT_matrix.tsv"
  "07_final_catalogs/catalog_81.DEL.germline.PASS.vcf.gz"
  "07_final_catalogs/catalog_81.DEL.germline.PASS.bed"
  "07_final_catalogs/catalog_81.DEL.germline.PASS.GT_matrix.tsv"
  "08_annotation/catalog_226.functional_class.tsv"
  "08_annotation/catalog_81.functional_class.tsv"
  "08_annotation/catalog_226.DELs_in_repeats.bed"
  "08_annotation/catalog_226.DELs_not_in_repeats.bed"
  "09_depth_support/depth_support_226.tsv"
  "10_mate_distance_qc/mate_distance_qc_226.tsv"
  "11_summary/per_sample_DEL_counts.tsv"
  "11_summary/per_chromosome_DEL_counts.tsv"
  "11_summary/private_vs_shared_DEL.tsv"
  "11_summary/DEL_svlen_distribution.tsv"
  "11_summary/DEL_window_counts_1Mb.tsv"
  "11_summary/pairwise_shared_DEL.tsv"
  "11_summary/DEL_binary_genotype_matrix.tsv"
  "11_summary/delly_DEL_counts.tsv"
  "11_summary/delly_DEL_summary_report.txt"
  "../00-samples/fClaHyb_Gar_LG.fa.fai"
  "samples_all_226.txt"
  "samples_unrelated_81.txt"
  "06_plot_delly_results.R"
  "05_summary_report.sh"
  "04_annotation_layers.sh"
)

human_size() {
  du -h "$1" 2>/dev/null | cut -f1
}

print_text_preview() {
  local f="$1"
  echo "Preview:"
  sed -n '1,8p' "$f" 2>/dev/null || true
}

print_gz_preview() {
  local f="$1"
  echo "Preview:"
  zcat "$f" 2>/dev/null | sed -n '1,12p' || true
}

print_vcf_preview() {
  local f="$1"
  echo "VCF header preview:"
  zcat "$f" 2>/dev/null | grep '^##' | sed -n '1,12p' || true
  echo
  echo "Column header + first 3 records:"
  zcat "$f" 2>/dev/null | awk 'BEGIN{c=0} /^#CHROM/{print; next} !/^#/{print; c++; if(c==3) exit}' || true
}

{
  echo "DEL PIPELINE RESULT SUMMARY FOR CLAUDE"
  echo "Generated: $(date '+%F %T')"
  echo "Working directory: $(pwd)"
  echo

  for f in "${files[@]}"; do
    echo "=================================================================="
    echo "FILE: $f"

    if [[ ! -e "$f" ]]; then
      echo "Status: MISSING"
      echo
      continue
    fi

    echo "Status: PRESENT"
    echo "Size: $(human_size "$f")"

    case "$f" in
      *.vcf.gz)
        echo "Type: gzipped VCF"
        print_vcf_preview "$f"
        ;;
      *.gz)
        echo "Type: gzipped text"
        print_gz_preview "$f"
        ;;
      *.tsv|*.bed|*.txt|*.sh|*.R|*.fai)
        echo "Type: text"
        echo "Lines: $(wc -l < "$f" 2>/dev/null || echo NA)"
        print_text_preview "$f"
        ;;
      *)
        echo "Type: other"
        ;;
    esac

    echo
  done

  echo "=================================================================="
  echo "INTERPRETATION NOTES"
  echo "- catalog_226.* = strict DEL catalog for all 226 samples"
  echo "- catalog_81.*  = strict DEL catalog for unrelated 81-sample panel"
  echo "- functional_class.tsv = DEL functional class (intergenic/intronic/exon/CDS)"
  echo "- DELs_in_repeats / DELs_not_in_repeats = repeat-overlap partitions"
  echo "- depth_support_226.tsv = depth-based support summary"
  echo "- mate_distance_qc_226.tsv = size/extreme-event QC summary"
  echo "- per_sample_DEL_counts.tsv = per-sample DEL burden"
  echo "- per_chromosome_DEL_counts.tsv = per-chromosome DEL counts"
  echo "- private_vs_shared_DEL.tsv = DEL frequency class summary"
  echo "- DEL_svlen_distribution.tsv = DEL size list for histogram/statistics"
  echo "- DEL_window_counts_1Mb.tsv = genome-window DEL counts"
  echo "- pairwise_shared_DEL.tsv = sample-by-sample shared DEL counts"
  echo "- DEL_binary_genotype_matrix.tsv = binary DEL presence/absence matrix for PCA/heatmaps"
  echo "- fai = chromosome sizes for genome plots"
  echo
  echo "CURRENT STRICT FILTER LOGIC"
  echo "- SVTYPE=DEL"
  echo "- PRECISE=1"
  echo "- QUAL >= 500"
  echo "- PE >= 5"
  echo "- PASS where required"
  echo
} > "$OUT"

echo "Wrote: $OUT"
