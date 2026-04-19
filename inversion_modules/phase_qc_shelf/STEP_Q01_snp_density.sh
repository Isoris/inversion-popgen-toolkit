#!/bin/bash
# =============================================================================
# STEP_Q01_snp_density.sh
# =============================================================================
# Build per-chromosome SNP density tracks from ANGSD outputs, binned at BIN_MB.
#
# Signal: low SNP density in a region indicates either (a) low-polymorphism
# segment (heterochromatin, conservation) or (b) mapping/calling failure.
# Both matter for diagnosing a Z plateau.
#
# Input:  ${BEAGLE_DIR}/main_qcpass.<CHR>.pos.fixed (or .pos)
#         These are 2-col files: chr<TAB>position (one line per called site)
#
# Output: ${QC_TRACKS}/snp_density.<CHR>.tsv
#         Columns: chrom  bin_start_bp  bin_end_bp  bin_mid_mb  n_snps  density_per_kb
#
# Usage:
#   bash STEP_Q01_snp_density.sh <CHR>       # single chrom
#   bash STEP_Q01_snp_density.sh all         # loop chr.list
# =============================================================================
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/00_config.sh"

CHR="${1:-}"
[[ -z "${CHR}" ]] && qc_die "Usage: $0 <CHR|all>"

run_one() {
  local chr="$1"
  local pos_file="${BEAGLE_DIR}/main_qcpass.${chr}.pos.fixed"
  [[ -f "${pos_file}" ]] || pos_file="${BEAGLE_DIR}/main_qcpass.${chr}.pos"
  [[ -f "${pos_file}" ]] || { qc_log "SKIP ${chr}: no .pos file"; return 0; }

  local out="${QC_TRACKS}/snp_density.${chr}.tsv"
  qc_log "Q01 ${chr}: ${pos_file} -> ${out}"

  # Auto-detect header, accept gzipped or plain
  local reader="cat"
  [[ "${pos_file}" == *.gz ]] && reader="${ZCAT_BIN}"

  # Binning with awk. BIN_BP = BIN_MB * 1e6
  awk -v OFS='\t' \
      -v chr="${chr}" \
      -v bin_bp=$(awk "BEGIN{printf \"%d\", ${BIN_MB}*1e6}") '
    NR==1 {
      # Detect header: if second field is non-numeric, skip
      if ($2 !~ /^[0-9]+$/) next
    }
    {
      # pos files: col 1 = chrom name (match), col 2 = position (1-based)
      if ($1 != chr) next
      p = $2 + 0
      b = int((p-1) / bin_bp)
      counts[b]++
      if (b > maxb) maxb = b
    }
    END {
      print "chrom", "bin_start_bp", "bin_end_bp", "bin_mid_mb", "n_snps", "density_per_kb"
      for (b = 0; b <= maxb; b++) {
        n = (b in counts) ? counts[b] : 0
        s = b * bin_bp + 1
        e = (b + 1) * bin_bp
        mid = (s + e) / 2 / 1e6
        density = n / (bin_bp / 1000)
        printf "%s\t%d\t%d\t%.4f\t%d\t%.4f\n", chr, s, e, mid, n, density
      }
    }
  ' <(${reader} "${pos_file}") > "${out}"

  local n_bins n_snps
  n_bins=$(( $(wc -l < "${out}") - 1 ))
  n_snps=$(awk 'NR>1 {s+=$5} END {print s+0}' "${out}")
  qc_log "Q01 ${chr}: ${n_bins} bins, ${n_snps} SNPs"
}

if [[ "${CHR}" == "all" ]]; then
  chr_list="${BEAGLE_DIR}/chr.list"
  [[ -f "${chr_list}" ]] || qc_die "chr.list not found: ${chr_list}"
  while IFS= read -r c; do
    [[ -z "${c}" ]] && continue
    run_one "${c}"
  done < "${chr_list}"
else
  run_one "${CHR}"
fi

qc_log "Q01 DONE. Tracks in ${QC_TRACKS}/"
