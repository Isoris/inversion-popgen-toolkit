#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 127
#SBATCH --mem=237G
#SBATCH -t 2-00:00:00
#SBATCH -J tra_geno
#SBATCH -o logs/03_merge_geno.%j.out
#SBATCH -e logs/03_merge_geno.%j.err
set -euo pipefail
source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${SCRIPT_DIR}/../00_module4f_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
set -a
source "${CONFIG}"
set +a

dv_init_dirs
dv_log "=== STEPS 3–8: Merge → Regenotype → Cohort → Germline filter (TRA) ==="

dv_check_file "$REF" "Reference"
dv_check_file "$EXCL_BED" "Exclusion BED"
dv_check_file "$SAMPLES_ALL" "Sample list (all)"
dv_check_file "$SAMPLES_UNRELATED" "Sample list (unrelated)"
[[ -x "${DELLY_BIN}" ]] || dv_die "DELLY not executable"

N_ALL=$(wc -l < "$SAMPLES_ALL")
N_UNREL=$(wc -l < "$SAMPLES_UNRELATED")
dv_log "Full: ${N_ALL} | Unrelated: ${N_UNREL}"

# ── STEP 3: delly merge ─────────────────────────────────────────────────────
dv_log "--- STEP 3: delly merge ---"

DISC_LIST="${DIR_SITES}/disc_bcf_list.txt"
find "${DIR_DISC}" -maxdepth 1 -name '*.disc.TRA.bcf' | sort > "${DISC_LIST}"

while IFS= read -r bcf; do
  [[ -f "${bcf}.csi" ]] || bcftools index -f "$bcf"
done < "${DISC_LIST}"

SITES_BCF="${DIR_SITES}/sites.TRA.bcf"
"${DELLY_BIN}" merge -o "${SITES_BCF}" $(cat "${DISC_LIST}")
bcftools index -f "${SITES_BCF}"

N_SITES=$(bcftools view -H "${SITES_BCF}" | wc -l)
dv_log "  Merged TRA sites: ${N_SITES}"

# ── STEP 4: Regenotype ──────────────────────────────────────────────────────
dv_log "--- STEP 4: Regenotype ${N_ALL} samples ---"

run_genotype() {
  local sid="$1"
  local bam="${DIR_MARKDUP}/${sid}.markdup.bam"
  local out="${DIR_GENO}/${sid}.geno.TRA.bcf"
  local logf="${DIR_LOGS}/geno_${sid}.log"

  if [[ -f "${out}" && -f "${out}.csi" ]]; then
    echo "[$(date '+%T')] ${sid}: skip"; return 0
  fi

  echo "[$(date '+%T')] ${sid}: regenotyping..."
  "${DELLY_BIN}" call \
    -t BND \
    -g "${REF}" \
    -x "${EXCL_BED}" \
    -v "${SITES_BCF}" \
    -h "${DELLY_THREADS_PER_CALL}" \
    -o "${out}" \
    "${bam}" 2>>"${logf}"

  bcftools index -f "${out}" 2>>"${logf}"
}

export -f run_genotype
export DIR_MARKDUP REF EXCL_BED DELLY_BIN DELLY_THREADS_PER_CALL SITES_BCF DIR_GENO DIR_LOGS

parallel -j "${DELLY_PARALLEL}" --joblog "${DIR_LOGS}/genotype_parallel.log" \
  run_genotype {} :::: "${SAMPLES_ALL}"

FAIL=0
while IFS= read -r sid; do
  [[ -f "${DIR_GENO}/${sid}.geno.TRA.bcf" ]] || { dv_err "Missing: ${sid}"; ((FAIL++)) || true; }
done < "$SAMPLES_ALL"
[[ ${FAIL} -eq 0 ]] || dv_die "${FAIL} failed regenotyping."

# ── STEP 5: bcftools merge 226 ──────────────────────────────────────────────
dv_log "--- STEP 5: bcftools merge (226 cohort) ---"
ls "${DIR_GENO}"/*.geno.TRA.bcf > "${DIR_MERGED}/geno_list.txt"

COHORT_226="${DIR_MERGED}/cohort_226.TRA.merged.bcf"
bcftools merge -m id -O b -o "${COHORT_226}" --threads "${THREADS}" \
  $(cat "${DIR_MERGED}/geno_list.txt")
bcftools index -f "${COHORT_226}"
dv_log "  226 cohort: $(bcftools view -H "${COHORT_226}" | wc -l) TRAs"

# ── STEP 6: Subset to 81 ────────────────────────────────────────────────────
dv_log "--- STEP 6: Subset 81 unrelated ---"
COHORT_81="${DIR_SUBSET81}/cohort_81.TRA.merged.bcf"
bcftools view -S "${SAMPLES_UNRELATED}" --force-samples -O b -o "${COHORT_81}" "${COHORT_226}"
bcftools index -f "${COHORT_81}"

COHORT_81_TRIM="${DIR_SUBSET81}/cohort_81.TRA.trimmed.bcf"
bcftools view -i 'COUNT(GT="alt")>0' -O b -o "${COHORT_81_TRIM}" "${COHORT_81}"
bcftools index -f "${COHORT_81_TRIM}"
dv_log "  81 subset: $(bcftools view -H "${COHORT_81_TRIM}" | wc -l) polymorphic TRAs"

# ── STEP 7: delly filter -f germline ────────────────────────────────────────
dv_log "--- STEP 7: Germline filter (81 unrelated) ---"
GERMLINE_81="${DIR_FILTERED}/cohort_81.TRA.germline.bcf"

if "${DELLY_BIN}" filter -f germline -o "${GERMLINE_81}" "${COHORT_81_TRIM}" \
  2>"${DIR_LOGS}/delly_filter_germline.log"; then
  bcftools index -f "${GERMLINE_81}"
  dv_log "  Germline PASS: $(bcftools view -H "${GERMLINE_81}" | wc -l) TRAs"
else
  dv_log "  WARNING: germline filter returned non-zero (common with family structure)"
fi

# ── Strict filter expression ────────────────────────────────────────────────
STRICT_EXPR='INFO/PRECISE=1 && QUAL>=300 && INFO/PE>=3'

# ── STEP 8: Final catalogs + GT matrices ────────────────────────────────────
dv_log "--- STEP 8: Final catalogs ---"
dv_log "  Strict filter: ${STRICT_EXPR}"

# 226 raw catalog
FINAL_226_RAW_VCF="${DIR_FINAL}/catalog_226.TRA.raw.vcf.gz"
bcftools view -O z -o "${FINAL_226_RAW_VCF}" "${COHORT_226}"
tabix -f -p vcf "${FINAL_226_RAW_VCF}"

# 226 strict catalog
FINAL_226_VCF="${DIR_FINAL}/catalog_226.TRA.vcf.gz"
if [[ "${STRICT_REQUIRE_PASS}" -eq 1 ]]; then
  bcftools view -f PASS -i "${STRICT_EXPR}" -O z -o "${FINAL_226_VCF}" "${COHORT_226}"
else
  bcftools view -i "${STRICT_EXPR}" -O z -o "${FINAL_226_VCF}" "${COHORT_226}"
fi
tabix -f -p vcf "${FINAL_226_VCF}"

MATRIX_226="${DIR_FINAL}/catalog_226.TRA.GT_matrix.tsv"
{
  echo -e "CHROM\tPOS\tCHR2\tID\tPOS2\t$(bcftools query -l "${FINAL_226_VCF}" | tr '\n' '\t' | sed 's/\t$//')"
  bcftools query -f '%CHROM\t%POS\t%INFO/CHR2\t%ID\t%INFO/POS2[\t%GT]\n' "${FINAL_226_VCF}"
} > "${MATRIX_226}"
dv_log "  226 strict catalog: $(tail -n +2 "${MATRIX_226}" | wc -l) TRAs"

# 81 source
if [[ -f "${GERMLINE_81}" ]]; then
  INPUT_81="${GERMLINE_81}"
  LABEL_81="germline"
else
  INPUT_81="${COHORT_81_TRIM}"
  LABEL_81="trimmed"
fi

# 81 raw catalog
FINAL_81_RAW_VCF="${DIR_FINAL}/catalog_81.TRA.${LABEL_81}.raw.vcf.gz"
bcftools view -O z -o "${FINAL_81_RAW_VCF}" "${INPUT_81}"
tabix -f -p vcf "${FINAL_81_RAW_VCF}"

# 81 strict catalog
FINAL_81_PASS="${DIR_FINAL}/catalog_81.TRA.${LABEL_81}.PASS.vcf.gz"
if [[ "${STRICT_REQUIRE_PASS}" -eq 1 ]]; then
  bcftools view -f PASS -i "${STRICT_EXPR}" -O z -o "${FINAL_81_PASS}" "${INPUT_81}"
else
  bcftools view -i "${STRICT_EXPR}" -O z -o "${FINAL_81_PASS}" "${INPUT_81}"
fi
tabix -f -p vcf "${FINAL_81_PASS}"

MATRIX_81="${DIR_FINAL}/catalog_81.TRA.${LABEL_81}.PASS.GT_matrix.tsv"
{
  echo -e "CHROM\tPOS\tCHR2\tID\tPOS2\t$(bcftools query -l "${FINAL_81_PASS}" | tr '\n' '\t' | sed 's/\t$//')"
  bcftools query -f '%CHROM\t%POS\t%INFO/CHR2\t%ID\t%INFO/POS2[\t%GT]\n' "${FINAL_81_PASS}"
} > "${MATRIX_81}"

# BEDs for annotation — BND breakpoints as ±500bp windows
for vcf in "${FINAL_226_VCF}" "${FINAL_81_PASS}"; do
  bed="${vcf%.vcf.gz}.bed"
  bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/CHR2\t%INFO/POS2\n' "${vcf}" \
    | awk 'BEGIN{OFS="\t"}{
        bp=500
        s1=$2-bp; if(s1<0) s1=0; e1=$2+bp
        print $1,s1,e1,$3,"breakpoint1"
        s2=$5-bp; if(s2<0) s2=0; e2=$5+bp
        print $4,s2,e2,$3,"breakpoint2"
      }' > "${bed}"
done

# ── TRA-specific: filter to inter-chromosomal events only ───────────────────
dv_log "--- Filtering to inter-chromosomal (TRA) events ---"
for vcf_file in "${FINAL_226_VCF}" "${FINAL_81_PASS}"; do
  if [[ -f "${vcf_file}" ]]; then
    tmp="${vcf_file}.tmp.vcf.gz"
    bcftools view -i 'INFO/CHR2!=CHROM' -O z -o "${tmp}" "${vcf_file}"
    mv "${tmp}" "${vcf_file}"
    tabix -f -p vcf "${vcf_file}"
    dv_log "  $(basename "${vcf_file}"): $(bcftools view -H "${vcf_file}" | wc -l) inter-chromosomal events"
  fi
done

dv_log "=== STEPS 3–8 COMPLETE ==="
