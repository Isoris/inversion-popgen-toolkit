#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 127
#SBATCH --mem=237G
#SBATCH -t 2-00:00:00
#SBATCH -J delly_dup_geno
#SBATCH -o logs/03_merge_geno.%j.out
#SBATCH -e logs/03_merge_geno.%j.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${SCRIPT_DIR}/../00_module4c_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }

set -a
source "${CONFIG}"
set +a

dv_init_dirs
dv_log "=== STEPS 3–8: Merge → Regenotype → Cohort → Germline filter (${SVTYPE}) ==="

dv_check_file "$REF" "Reference"
dv_check_file "$EXCL_BED" "Exclusion BED"
dv_check_file "$SAMPLES_ALL" "Sample list (all)"
dv_check_file "$SAMPLES_UNRELATED" "Sample list (unrelated)"
[[ -x "${DELLY_BIN}" ]] || dv_die "DELLY not executable"

# STEP 3: merge discovered sites
dv_log "--- STEP 3: delly merge ---"
DISC_LIST="${DIR_SITES}/disc_bcf_list.txt"
find "${DIR_DISC}" -maxdepth 1 -name "*.disc.${SVTYPE}.bcf" | sort > "${DISC_LIST}"

while IFS= read -r bcf; do
  [[ -f "${bcf}.csi" ]] || bcftools index -f "$bcf"
done < "${DISC_LIST}"

SITES_BCF="${DIR_SITES}/sites.${SVTYPE}.bcf"
"${DELLY_BIN}" merge -o "${SITES_BCF}" $(cat "${DISC_LIST}")
bcftools index -f "${SITES_BCF}"

N_SITES=$(bcftools view -H "${SITES_BCF}" | wc -l)
dv_log "Merged ${SVTYPE} sites: ${N_SITES}"

# STEP 4: regenotype all
dv_log "--- STEP 4: Regenotype all samples ---"

run_genotype() {
  local sid="$1"
  local bam="${DIR_MARKDUP}/${sid}.markdup.bam"
  local out="${DIR_GENO}/${sid}.geno.${SVTYPE}.bcf"
  local logf="${DIR_LOGS}/geno_${sid}.log"

  if [[ -f "${out}" && -f "${out}.csi" ]]; then
    echo "[$(date '+%T')] ${sid}: skip"
    return 0
  fi

  "${DELLY_BIN}" call \
    -t "${SVTYPE}" \
    -g "${REF}" \
    -x "${EXCL_BED}" \
    -v "${SITES_BCF}" \
    -h "${DELLY_THREADS_PER_CALL}" \
    -o "${out}" \
    "${bam}" \
    2>>"${logf}"

  bcftools index -f "${out}" 2>>"${logf}"
}

export -f run_genotype
export DIR_MARKDUP REF EXCL_BED DELLY_BIN DELLY_THREADS_PER_CALL SITES_BCF DIR_GENO DIR_LOGS SVTYPE

parallel -j "${DELLY_PARALLEL}" --joblog "${DIR_LOGS}/genotype_parallel.log" \
  run_genotype {} :::: "${SAMPLES_ALL}"

FAIL=0
while IFS= read -r sid; do
  [[ -f "${DIR_GENO}/${sid}.geno.${SVTYPE}.bcf" ]] || { dv_err "Missing regenotype output for ${sid}"; ((FAIL++)) || true; }
done < "$SAMPLES_ALL"
[[ ${FAIL} -eq 0 ]] || dv_die "${FAIL} failed regenotyping."

# STEP 5: merge cohort 226
dv_log "--- STEP 5: bcftools merge 226 cohort ---"
ls "${DIR_GENO}"/*.geno."${SVTYPE}".bcf > "${DIR_MERGED}/geno_list.txt"

COHORT_226="${DIR_MERGED}/cohort_226.${SVTYPE}.merged.bcf"
bcftools merge -m id -O b -o "${COHORT_226}" --threads "${THREADS}" \
  $(cat "${DIR_MERGED}/geno_list.txt")
bcftools index -f "${COHORT_226}"

dv_log "226 cohort: $(bcftools view -H "${COHORT_226}" | wc -l) ${SVTYPE}"

# STEP 6: subset to 81
dv_log "--- STEP 6: Subset 81 unrelated ---"
COHORT_81="${DIR_SUBSET81}/cohort_81.${SVTYPE}.merged.bcf"
bcftools view -S "${SAMPLES_UNRELATED}" --force-samples -O b -o "${COHORT_81}" "${COHORT_226}"
bcftools index -f "${COHORT_81}"

COHORT_81_TRIM="${DIR_SUBSET81}/cohort_81.${SVTYPE}.trimmed.bcf"
bcftools view -i 'COUNT(GT="alt")>0' -O b -o "${COHORT_81_TRIM}" "${COHORT_81}"
bcftools index -f "${COHORT_81_TRIM}"

dv_log "81 subset: $(bcftools view -H "${COHORT_81_TRIM}" | wc -l) polymorphic ${SVTYPE}"

# STEP 7: germline filter
dv_log "--- STEP 7: Germline filter ---"
GERMLINE_81="${DIR_FILTERED}/cohort_81.${SVTYPE}.germline.bcf"

if "${DELLY_BIN}" filter -f germline -o "${GERMLINE_81}" "${COHORT_81_TRIM}" \
  2>"${DIR_LOGS}/delly_filter_germline.log"; then
  bcftools index -f "${GERMLINE_81}"
  dv_log "Germline output: $(bcftools view -H "${GERMLINE_81}" | wc -l) ${SVTYPE}"
else
  dv_log "WARNING: germline filter returned non-zero"
fi

# STEP 8: final catalogs
dv_log "--- STEP 8: Final catalogs ---"

STRICT_EXPR='INFO/SVTYPE="'${SVTYPE}'" && QUAL>='${STRICT_MIN_QUAL}' && INFO/PE>='${STRICT_MIN_PE}
if [[ "${STRICT_REQUIRE_PRECISE}" -eq 1 ]]; then
  STRICT_EXPR="${STRICT_EXPR} && INFO/PRECISE=1"
fi

FINAL_226_RAW_VCF="${DIR_FINAL}/catalog_226.${SVTYPE}.raw.vcf.gz"
bcftools view -i 'INFO/SVTYPE="'${SVTYPE}'"' -O z -o "${FINAL_226_RAW_VCF}" "${COHORT_226}"
tabix -f -p vcf "${FINAL_226_RAW_VCF}"

FINAL_226_VCF="${DIR_FINAL}/catalog_226.${SVTYPE}.vcf.gz"
if [[ "${STRICT_REQUIRE_PASS}" -eq 1 ]]; then
  bcftools view -f PASS -i "${STRICT_EXPR}" -O z -o "${FINAL_226_VCF}" "${COHORT_226}"
else
  bcftools view -i "${STRICT_EXPR}" -O z -o "${FINAL_226_VCF}" "${COHORT_226}"
fi
tabix -f -p vcf "${FINAL_226_VCF}"

MATRIX_226="${DIR_FINAL}/catalog_226.${SVTYPE}.GT_matrix.tsv"
{
  echo -e "CHROM\tPOS\tEND\tID\tSVLEN\t$(bcftools query -l "${FINAL_226_VCF}" | tr '\n' '\t' | sed 's/\t$//')"
  bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVLEN[\t%GT]\n' "${FINAL_226_VCF}"
} > "${MATRIX_226}"

if [[ -f "${GERMLINE_81}" ]]; then
  INPUT_81="${GERMLINE_81}"
  LABEL_81="germline"
else
  INPUT_81="${COHORT_81_TRIM}"
  LABEL_81="trimmed"
fi

FINAL_81_RAW_VCF="${DIR_FINAL}/catalog_81.${SVTYPE}.${LABEL_81}.raw.vcf.gz"
bcftools view -i 'INFO/SVTYPE="'${SVTYPE}'"' -O z -o "${FINAL_81_RAW_VCF}" "${INPUT_81}"
tabix -f -p vcf "${FINAL_81_RAW_VCF}"

FINAL_81_PASS="${DIR_FINAL}/catalog_81.${SVTYPE}.${LABEL_81}.PASS.vcf.gz"
if [[ "${STRICT_REQUIRE_PASS}" -eq 1 ]]; then
  bcftools view -f PASS -i "${STRICT_EXPR}" -O z -o "${FINAL_81_PASS}" "${INPUT_81}"
else
  bcftools view -i "${STRICT_EXPR}" -O z -o "${FINAL_81_PASS}" "${INPUT_81}"
fi
tabix -f -p vcf "${FINAL_81_PASS}"

MATRIX_81="${DIR_FINAL}/catalog_81.${SVTYPE}.${LABEL_81}.PASS.GT_matrix.tsv"
{
  echo -e "CHROM\tPOS\tEND\tID\tSVLEN\t$(bcftools query -l "${FINAL_81_PASS}" | tr '\n' '\t' | sed 's/\t$//')"
  bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVLEN[\t%GT]\n' "${FINAL_81_PASS}"
} > "${MATRIX_81}"

for vcf in "${FINAL_226_VCF}" "${FINAL_81_PASS}"; do
  bed="${vcf%.vcf.gz}.bed"
  bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVLEN\n' "${vcf}" \
    | awk 'BEGIN{OFS="\t"}{s=$2;e=$3;if(e=="."||e==""){e=s+1} if(s>e){t=s;s=e;e=t} print $1,s,e,$4,$5}' > "${bed}"
done

dv_log "=== STEPS 3–8 COMPLETE ==="
