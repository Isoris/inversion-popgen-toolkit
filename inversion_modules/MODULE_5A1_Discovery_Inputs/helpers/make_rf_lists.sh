BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
INVDIR="${BASE}/inversion_localpca_v7"

mkdir -p "${INVDIR}/chunks/rf_files"

while read -r chr; do
  printf "%s\n" "$chr" > "${INVDIR}/chunks/rf_files/${chr}.rf.txt"
done < "${BASE}/het_roh_workflow/chr.list"

find "${INVDIR}/chunks/rf_files" -name "*.rf.txt" | sort > "${INVDIR}/chunk_rf.list"

wc -l "${INVDIR}/chunk_rf.list"
head "${INVDIR}/chunk_rf.list"
