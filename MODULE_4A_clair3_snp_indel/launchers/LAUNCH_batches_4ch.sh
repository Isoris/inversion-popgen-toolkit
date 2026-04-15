cd /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/MODULE_4A_SNP_INDEL50_Clair3

mkdir -p meta/chr_batches_4
rm -f meta/chr_batches_4/*.txt

awk '
{
    file=sprintf("meta/chr_batches_4/chr_batch_%02d.txt", int((NR-1)/4)+1)
    print $0 >> file
}
' meta/chromosome_list.noLG01.txt

for f in meta/chr_batches_4/*.txt; do
    echo "=== $f ==="
    cat "$f"
done
