#!/bin/bash
# =============================================================================
#  STEP_Q01_snp_density.sh
#  -------------------------------------------------------------------
#  Per-chromosome SNP density tracks from ANGSD called-site positions.
# =============================================================================
#
#  PURPOSE   Bin per-chromosome ANGSD called-site positions into fixed-width
#            windows; emit a SNPs-per-window track for the Q04 diagnostic plotter.
#  INPUTS    ${BEAGLE_DIR}/main_qcpass.<CHR>.pos.fixed   (2-col: chr, 1-based pos)
#            ${BEAGLE_DIR}/chr.list                       (needed when CHR=all)
#  OUTPUTS   ${QC_TRACKS}/snp_density.<CHR>.tsv
#            cols: chrom  bin_start_bp  bin_end_bp  bin_mid_mb  n_snps  density_per_kb
#  CONFIG    BIN_MB  (from 00_config.sh, typical 0.05 = 50 kb bins)
#  USAGE     bash STEP_Q01_snp_density.sh <CHR>      # single chromosome
#            bash STEP_Q01_snp_density.sh all        # loop chr.list
#  STATUS    STABLE              LAST UPDATED   2026-04-20
#  CALLED BY run_chrom.sh · run_all.sh · slurm/array_28chrom.sh
#
# -------------------------------------------------------------------
#  What this step is doing
# -------------------------------------------------------------------
#  This is Q01 — the first QC step in phase_qc_shelf. It reduces the per-site
#  position list from ANGSD into a binned density track, one row per bin. The
#  density signal later overlays the chromosome ideogram in Q04 so a reader
#  can instantly see where SNPs are sparse vs dense.
#
#  Low-density bins have two biological readings, and the paper methods
#  must discriminate them: (a) genuine low-polymorphism biology — highly
#  conserved exons, heterochromatin, recent selective sweeps — versus (b)
#  mapping/calling failure — repeat-masked regions, reference-N blocks,
#  low-mappability sequence. Interpreting a local-PCA Z-plateau or a Fst
#  outlier requires deciding which of the two is in play at that locus.
# =============================================================================
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/00_config.sh"

CHR="${1:-}"
[[ -z "${CHR}" ]] && qc_die "Usage: $0 <CHR|all>"

# Shared helpers available from 00_config.sh:
#   qc_log / qc_warn / qc_err / qc_die
#   qc_banner_open / qc_banner_close
#   qc_config_snapshot <STEP> VAR1 VAR2 ...
#   qc_preview_file <label> <path> [n_lines]
#   qc_log_grep_tips

qc_config_snapshot "Q01" BIN_MB ZCAT_BIN

run_one() {
  local chr="$1"
  local t0 t1 elapsed
  t0=$(date +%s)

  local pos_file="${BEAGLE_DIR}/main_qcpass.${chr}.pos.fixed"
  [[ -f "${pos_file}" ]] || pos_file="${BEAGLE_DIR}/main_qcpass.${chr}.pos"
  [[ -f "${pos_file}" ]] || { qc_warn "SKIP ${chr}: no .pos file in ${BEAGLE_DIR}"; return 0; }

  local out="${QC_TRACKS}/snp_density.${chr}.tsv"

  qc_banner_open "Q01" "${chr}"
  qc_preview_file "input" "${pos_file}" 3

  # Auto-detect gzipped or plain
  local reader="cat"
  [[ "${pos_file}" == *.gz ]] && reader="${ZCAT_BIN}"

  # Binning with awk. BIN_BP = BIN_MB * 1e6
  awk -v OFS='\t' \
      -v chr="${chr}" \
      -v bin_bp=$(awk "BEGIN{printf \"%d\", ${BIN_MB}*1e6}") '
    NR==1 {
      if ($2 !~ /^[0-9]+$/) next
    }
    {
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

  local n_bins n_snps mean_per_bin
  n_bins=$(( $(wc -l < "${out}") - 1 ))
  n_snps=$(awk 'NR>1 {s+=$5} END {print s+0}' "${out}")
  mean_per_bin=$(awk -v n="${n_snps}" -v b="${n_bins}" 'BEGIN{ if (b>0) printf "%d", n/b; else print 0 }')
  qc_log "📊 summary  ${n_bins} bins · ${n_snps} SNPs · mean ${mean_per_bin}/bin"

  qc_preview_file "output" "${out}" 3

  t1=$(date +%s)
  elapsed=$(( t1 - t0 ))
  qc_banner_close "Q01" "${chr}" "${elapsed}"
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

qc_log "🏁 Q01 DONE. Tracks in ${QC_TRACKS}/"
qc_log_grep_tips
