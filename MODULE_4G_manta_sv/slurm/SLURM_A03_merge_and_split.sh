#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 127
#SBATCH --mem=237G
#SBATCH -t 2-00:00:00
#SBATCH -J manta_merge
#SBATCH -o logs/03_merge_split.%j.out
#SBATCH -e logs/03_merge_split.%j.err
# =============================================================================
# 03_merge_and_split.sh — Convert INV per-sample → Merge → Subset → Split
# =============================================================================
#
# CRITICAL: convertInversion.py must run PER-SAMPLE, NOT on the merged VCF.
#
# After bcftools merge, the VCF ID column contains compound IDs (semicolon-
# separated from multiple samples), e.g.:
#   MantaBND:148628:0:1:0:0:0:0;MantaBND:59701:0:1:0:0:0:1;...
#
# convertInversion.py's scanVcf() uses vcfRec.vid to build the invMateDict,
# and later does invMateDict[mateId]["CIPOS"]. With compound IDs, the mate
# lookup returns an empty string "" instead of a dict → TypeError crash:
#   TypeError: string indices must be integers, not str
#
# Solution: convert inversions per-sample BEFORE merging, then merge the
# already-converted VCFs. This is the standard approach (nf-core, CGAP, etc).
#
# Pipeline order:
#   1. Per-sample: convertInversion.py on each diploidSV.vcf.gz
#   2. bcftools merge all 226 converted VCFs
#   3. Subset to 81 unrelated
#   4. Split into per-type catalogs
# =============================================================================
set -euo pipefail
source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${SCRIPT_DIR}/../00_module4g_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
set -a
source "${CONFIG}"
set +a

CONVERT_INV_PY3="${SCRIPT_DIR}/convertInversion_py3.py"
mv_init_dirs
mv_log "=== STEP 3: Convert INV (per-sample) → Merge → Split by SV type ==="

mv_check_file "$SAMPLES_ALL" "Sample list (all)"
mv_check_file "$SAMPLES_UNRELATED" "Sample list (unrelated)"
mv_check_file "${CONVERT_INV_PY3}" "Python3-compatible convertInversion.py"

# Locate tools
SAMTOOLS_BIN=$(command -v samtools)
PYTHON2_BIN=$(command -v python2 || command -v python)
mv_log "samtools: ${SAMTOOLS_BIN}"
mv_log "python2: ${PYTHON2_BIN}"
mv_log "convertInversion.py: ${CONVERT_INV_PY3}"

N_ALL=$(wc -l < "$SAMPLES_ALL")
N_UNREL=$(wc -l < "$SAMPLES_UNRELATED")
mv_log "Full: ${N_ALL} | Unrelated: ${N_UNREL}"

# ── STEP 3A: Per-sample inversion conversion ───────────────────────────────
# convertInversion.py works correctly on single-sample VCFs where each BND
# record has a clean single MantaBND:... ID. It crashes on merged VCFs
# because bcftools merge concatenates IDs with semicolons.
# ─────────────────────────────────────────────────────────────────────────────
mv_log "--- 3A: Per-sample BND → INV conversion ---"

DIR_CONVERTED="${OUTDIR}/01b_per_sample_converted"
mkdir -p "${DIR_CONVERTED}"

convert_one_sample() {
  local sid="$1"
  local raw_vcf="${DIR_PERSAMPLE}/${sid}/results/variants/diploidSV.vcf.gz"
  local conv_vcf="${DIR_CONVERTED}/${sid}.diploidSV.inv_converted.vcf.gz"
  local logf="${DIR_LOGS}/convert_inv_${sid}.log"

  # Skip if already done and valid
  if [[ -f "${conv_vcf}" && -f "${conv_vcf}.tbi" ]]; then
    local ncheck
    ncheck=$(bcftools view -H "${conv_vcf}" 2>/dev/null | head -1 | wc -l)
    if [[ "${ncheck}" -gt 0 ]]; then
      echo "[$(date '+%T')] ${sid}: already converted, skipping"
      return 0
    else
      # Previous run produced corrupt output — redo
      rm -f "${conv_vcf}" "${conv_vcf}.tbi"
    fi
  fi

  if [[ ! -f "${raw_vcf}" ]]; then
    echo "[$(date '+%T')] ${sid}: ERROR — diploidSV.vcf.gz missing" >&2
    return 1
  fi

  echo "[$(date '+%T')] ${sid}: converting inversions..."

  # IMPORTANT: convertInversion.py was written for Python 2 where bytes == str.
  # Under Python 3, gzip.open('rb') returns bytes but the script does string
  # operations (.strip(), .split('\t'), .startswith('#')) → TypeError.
  # Even the plain-text path uses open(file, 'rb') which returns bytes in Py3.
  #
  # Workaround: decompress to plain .vcf, then use our local Py3-patched
  # copy (convertInversion_py3.sh) which opens files in text mode.
  local tmp_vcf="${DIR_CONVERTED}/${sid}.tmp.vcf"

  zcat "${raw_vcf}" > "${tmp_vcf}"

  python3 "${CONVERT_INV_PY3}" \
    "${SAMTOOLS_BIN}" \
    "${REF}" \
    "${tmp_vcf}" \
    2>>"${logf}" \
    | bgzip -c > "${conv_vcf}"

  # Clean up temp
  rm -f "${tmp_vcf}"

  # Verify output is valid
  if [[ ! -s "${conv_vcf}" ]]; then
    echo "[$(date '+%T')] ${sid}: ERROR — empty output" >&2
    return 1
  fi

  tabix -f -p vcf "${conv_vcf}"

  local n_total n_inv
  n_total=$(bcftools view -H "${conv_vcf}" 2>/dev/null | wc -l)
  n_inv=$(bcftools view -H "${conv_vcf}" 2>/dev/null \
    | awk -F'\t' '{if($8 ~ /SVTYPE=INV/) print}' | wc -l)
  echo "[$(date '+%T')] ${sid}: ${n_total} SVs total, ${n_inv} INV"
}

export -f convert_one_sample
export DIR_PERSAMPLE DIR_CONVERTED DIR_LOGS SAMTOOLS_BIN REF
# Use our local Py3-compatible convertInversion script
CONVERT_INV_PY3="${SCRIPT_DIR}/convertInversion_py3.py"
mv_check_file "${CONVERT_INV_PY3}" "Python3-compatible convertInversion.py"
export CONVERT_INV_PY3

parallel -j "${MANTA_PARALLEL}" --joblog "${DIR_LOGS}/convert_inv_parallel.log" \
  convert_one_sample {} \
  :::: "${SAMPLES_ALL}"

# Verify
FAIL=0
while IFS= read -r sid; do
  conv="${DIR_CONVERTED}/${sid}.diploidSV.inv_converted.vcf.gz"
  [[ -f "${conv}" && -f "${conv}.tbi" ]] || { mv_err "Missing converted VCF: ${sid}"; ((FAIL++)) || true; }
done < "${SAMPLES_ALL}"
[[ ${FAIL} -eq 0 ]] || mv_die "${FAIL} samples failed inversion conversion"

# Count total inversions across all samples
TOTAL_INV=$(for f in "${DIR_CONVERTED}"/*.inv_converted.vcf.gz; do
  bcftools view -H "$f" 2>/dev/null | awk -F'\t' '$8 ~ /SVTYPE=INV/'
done | wc -l)
mv_log "  Total INV records across all samples: ${TOTAL_INV}"

# ── STEP 3B: Build VCF list for merge (using converted VCFs) ───────────────
mv_log "--- 3B: Collecting converted per-sample VCFs ---"
VCF_LIST="${DIR_MERGED}/vcf_list.txt"
> "${VCF_LIST}"
while IFS= read -r sid; do
  conv="${DIR_CONVERTED}/${sid}.diploidSV.inv_converted.vcf.gz"
  echo "${conv}" >> "${VCF_LIST}"
done < "${SAMPLES_ALL}"

N_VCFS=$(wc -l < "${VCF_LIST}")
mv_log "  ${N_VCFS} converted VCFs to merge"

# ── STEP 3C: bcftools merge 226 ────────────────────────────────────────────
mv_log "--- 3C: bcftools merge (226 cohort) ---"
COHORT_226="${DIR_MERGED}/cohort_226.ALL.inv_converted.vcf.gz"

bcftools merge \
  -m none \
  --threads "${THREADS}" \
  -O z -o "${COHORT_226}" \
  --file-list "${VCF_LIST}"

tabix -f -p vcf "${COHORT_226}"
N_MERGED=$(bcftools view -H "${COHORT_226}" | wc -l)
mv_log "  226 merged cohort: ${N_MERGED} SV records"

# Quick type breakdown
mv_log "  Type breakdown:"
bcftools view -H "${COHORT_226}" | awk -F'\t' '{
  n=split($8,a,";")
  for(i=1;i<=n;i++) if(a[i] ~ /^SVTYPE=/) {sub(/^SVTYPE=/,"",a[i]); print a[i]; break}
}' | sort | uniq -c | sort -rn | while read cnt svt; do
  mv_log "    ${svt}: ${cnt}"
done

# ── STEP 3D: Subset to 81 unrelated ────────────────────────────────────────
mv_log "--- 3D: Subset to 81 unrelated ---"
COHORT_81="${DIR_SUBSET81}/cohort_81.ALL.vcf.gz"
bcftools view -S "${SAMPLES_UNRELATED}" --force-samples \
  -O z -o "${COHORT_81}" "${COHORT_226}"
tabix -f -p vcf "${COHORT_81}"

# Trim monomorphic sites in the 81 subset
COHORT_81_TRIM="${DIR_SUBSET81}/cohort_81.ALL.trimmed.vcf.gz"
bcftools view -i 'COUNT(GT="alt")>0' -O z -o "${COHORT_81_TRIM}" "${COHORT_81}"
tabix -f -p vcf "${COHORT_81_TRIM}"
N_81=$(bcftools view -H "${COHORT_81_TRIM}" | wc -l)
mv_log "  81 polymorphic: ${N_81} SVs"

# ── STEP 3E: Split by SV type ──────────────────────────────────────────────
mv_log "--- 3E: Splitting by SV type ---"

split_by_type() {
  local input_vcf="$1"
  local label="$2"
  local outdir="$3"

  mkdir -p "${outdir}"

  # DEL
  bcftools view -i 'INFO/SVTYPE="DEL"' -O z \
    -o "${outdir}/${label}.DEL.vcf.gz" "${input_vcf}"
  tabix -f -p vcf "${outdir}/${label}.DEL.vcf.gz"

  # DUP (Manta reports as DUP:TANDEM — match both)
  bcftools view -i 'INFO/SVTYPE="DUP" || INFO/SVTYPE="DUP:TANDEM"' -O z \
    -o "${outdir}/${label}.DUP.vcf.gz" "${input_vcf}"
  tabix -f -p vcf "${outdir}/${label}.DUP.vcf.gz"

  # INV (converted from BND by convertInversion.py — now has SVTYPE=INV)
  bcftools view -i 'INFO/SVTYPE="INV"' -O z \
    -o "${outdir}/${label}.INV.vcf.gz" "${input_vcf}"
  tabix -f -p vcf "${outdir}/${label}.INV.vcf.gz"

  # BND — true translocations (remaining after INV extraction)
  bcftools view -i 'INFO/SVTYPE="BND"' -O z \
    -o "${outdir}/${label}.BND.vcf.gz" "${input_vcf}"
  tabix -f -p vcf "${outdir}/${label}.BND.vcf.gz"

  # INS — all insertions first, then split into small vs large
  bcftools view -i 'INFO/SVTYPE="INS"' -O z \
    -o "${outdir}/${label}.INS_all.vcf.gz" "${input_vcf}"
  tabix -f -p vcf "${outdir}/${label}.INS_all.vcf.gz"

  # INS_small: fully assembled insertions.
  # These have SVLEN set (>0) and the full inserted sequence in ALT.
  bcftools view -i 'INFO/SVTYPE="INS" && INFO/SVLEN>0' -O z \
    -o "${outdir}/${label}.INS_small.vcf.gz" "${input_vcf}"
  tabix -f -p vcf "${outdir}/${label}.INS_small.vcf.gz"

  # INS_large: incompletely assembled insertions.
  # These have LEFT_SVINSSEQ and/or RIGHT_SVINSSEQ in INFO.
  # ALT is <INS>, SVLEN typically absent or 0.
  bcftools view -i 'INFO/SVTYPE="INS" && (INFO/LEFT_SVINSSEQ!="." || INFO/RIGHT_SVINSSEQ!=".")' -O z \
    -o "${outdir}/${label}.INS_large.vcf.gz" "${input_vcf}"
  tabix -f -p vcf "${outdir}/${label}.INS_large.vcf.gz"

  # Count everything
  for svtype in DEL DUP INV BND INS_all INS_small INS_large; do
    local f="${outdir}/${label}.${svtype}.vcf.gz"
    if [[ -f "$f" ]]; then
      local n
      n=$(bcftools view -H "$f" 2>/dev/null | wc -l)
      mv_log "  ${label} ${svtype}: ${n}"
    fi
  done
}

mv_log "  Splitting 226 cohort:"
split_by_type "${COHORT_226}" "catalog_226" "${DIR_SPLIT}"

mv_log "  Splitting 81 trimmed:"
split_by_type "${COHORT_81_TRIM}" "catalog_81" "${DIR_SPLIT}"

# ── STEP 3F: PASS-filtered final catalogs ──────────────────────────────────
mv_log "--- 3F: Building PASS-filtered final catalogs ---"

for cohort_label in "catalog_226" "catalog_81"; do
  for svtype in DEL DUP INV BND INS_small INS_large; do
    src="${DIR_SPLIT}/${cohort_label}.${svtype}.vcf.gz"
    dst="${DIR_FINAL}/${cohort_label}.${svtype}.PASS.vcf.gz"
    if [[ -f "${src}" ]]; then
      if [[ "${STRICT_REQUIRE_PASS}" -eq 1 ]]; then
        bcftools view -f PASS -i "QUAL>=${STRICT_MIN_QUAL}" \
          -O z -o "${dst}" "${src}"
      else
        bcftools view -i "QUAL>=${STRICT_MIN_QUAL}" \
          -O z -o "${dst}" "${src}"
      fi
      tabix -f -p vcf "${dst}"
      n=$(bcftools view -H "${dst}" 2>/dev/null | wc -l)
      mv_log "  ${cohort_label} ${svtype} PASS: ${n}"
    fi
  done
done

# ── STEP 3G: GT matrices + BEDs for each type ─────────────────────────────
mv_log "--- 3G: Building GT matrices and BED files ---"

for cohort_label in "catalog_226" "catalog_81"; do
  for svtype in DEL DUP INV INS_small INS_large; do
    vcf="${DIR_FINAL}/${cohort_label}.${svtype}.PASS.vcf.gz"
    [[ -f "${vcf}" ]] || continue

    n=$(bcftools view -H "${vcf}" 2>/dev/null | wc -l)
    [[ ${n} -gt 0 ]] || continue

    # GT matrix
    matrix="${DIR_FINAL}/${cohort_label}.${svtype}.PASS.GT_matrix.tsv"
    {
      samples_line=$(bcftools query -l "${vcf}" | tr '\n' '\t' | sed 's/\t$//')
      echo -e "CHROM\tPOS\tEND\tID\tSVLEN\tSVTYPE\t${samples_line}"
      bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVLEN\t%INFO/SVTYPE[\t%GT]\n' "${vcf}"
    } > "${matrix}"

    # BED
    bed="${DIR_FINAL}/${cohort_label}.${svtype}.PASS.bed"
    bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVLEN\t%INFO/SVTYPE\n' "${vcf}" \
      | awk 'BEGIN{OFS="\t"}{s=$2;e=$3;if(s>e){t=s;s=e;e=t}; if(e==s)e=s+1; print $1,s,e,$4,$5,$6}' \
      > "${bed}"
    mv_log "  ${cohort_label} ${svtype}: ${n} → matrix + BED"
  done

  # BND: output breakpoint table (no BED since no meaningful END)
  bnd_vcf="${DIR_FINAL}/${cohort_label}.BND.PASS.vcf.gz"
  if [[ -f "${bnd_vcf}" ]]; then
    n=$(bcftools view -H "${bnd_vcf}" 2>/dev/null | wc -l)
    if [[ ${n} -gt 0 ]]; then
      bnd_tsv="${DIR_FINAL}/${cohort_label}.BND.PASS.breakpoints.tsv"
      {
        samples_line=$(bcftools query -l "${bnd_vcf}" | tr '\n' '\t' | sed 's/\t$//')
        echo -e "CHROM\tPOS\tID\tALT\tMATEID\tEVENT\t${samples_line}"
        bcftools query -f '%CHROM\t%POS\t%ID\t%ALT\t%INFO/MATEID\t%INFO/EVENT[\t%GT]\n' "${bnd_vcf}"
      } > "${bnd_tsv}"
      mv_log "  ${cohort_label} BND: ${n} breakpoints"
    fi
  fi
done

mv_log "=== STEP 3 COMPLETE ==="
