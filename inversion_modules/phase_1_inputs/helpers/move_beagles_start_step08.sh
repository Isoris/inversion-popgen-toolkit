BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
BEAGLEDIR="${BASE}/het_roh/01_inputs_check/beagle_by_chr"
INVDIR="${BASE}/inversion_localpca_v7"
SCRIPTDIR="${BASE}/inversion_codebase_v8.4/MODULE_5A2_Discovery_Core"

mkdir -p "${INVDIR}/04_dosage_by_chr" "${INVDIR}/05_local_pca" "${INVDIR}/06_mds_candidates"

for beagle in "${BEAGLEDIR}"/main_qcpass.C_gar_LG*.beagle.gz; do
  python3 "${SCRIPTDIR}/steps/STEP_A01_beagle_to_dosage_by_chr.py" \
    --beagle "${beagle}" \
    --outdir "${INVDIR}/04_dosage_by_chr"
done
