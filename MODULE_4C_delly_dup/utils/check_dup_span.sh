bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\n' 07_final_catalogs/catalog_226.DUP.vcf.gz \
| awk 'BEGIN{OFS="\t"}{
    len=$3-$2;
    if(len<0) len=-len;
    print $1,$2,$3,$4,len
}' > dup_span_from_end.tsv

awk '
{a[NR]=$5}
END{
  if(NR==0){print "no variants"; exit}
  asort(a)
  print "n=" NR
  print "min=" a[1]
  print "median=" a[int((NR+1)/2)]
  print "max=" a[NR]
}' dup_span_from_end.tsv
