# create_tabix_gtex_sigpairs.sh
# by Brian Chen
# This script takes GTEx v8 eQTL significant pair files and outputs index file for tabix
# $1 is path to tissue file

module load vcftools
module load htslib

tissueName=$(echo "${1##*/}" | cut -f 1,2 -d '.')
echo $tissueName

zcat $1 | sed '1D' | tr '_' '\t' | awk 'BEGIN {OFS = "\t"} {$2 = ($2 - 1)"\t"$2; print;}' | vcf-sort -p 8 | bgzip -c > $tissueName".v8.signif_variant_gene_pairs.tab.gz"

echo $tissueName".v8.signif_variant_gene_pairs.tab.gz created"

#create index file
tabix -p bed $tissueName".v8.signif_variant_gene_pairs.tab.gz" 

