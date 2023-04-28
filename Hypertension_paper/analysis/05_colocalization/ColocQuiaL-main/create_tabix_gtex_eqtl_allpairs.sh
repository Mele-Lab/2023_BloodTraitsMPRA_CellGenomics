# create_tabix_gtex_allpairs.sh
# by Brian Chen
# This script takes GTEx v8 eQTL all pair files and outputs index file for tabix
# $1 is path to tissue file

module load htslib
module load vcftools

tissueName=$(echo "${1##*/}" | cut -f 1,2 -d '.')
echo $tissueName

#zcat $1 | sed '1D' | tr '_' '\t' | awk 'BEGIN {OFS = "\t"} {t = $1; $1 = $2; $2 = $3 - 1; $4 = t"\t"$4; print;}' | vcf-sort -p 8 | bgzip -c > hg19/$tissueName".tab.gz"

#echo $tissueName".tab.gz created"

##gzip -v -d $tissueName".tab.gz" 

##rm $tissueName".tab.gz"

cat $tissueName".tab" | vcf-sort -p 8  | bgzip -c > $tissueName".tab.gz"

#bgzip -c $tissueName"sorted.tab" > $tissueName".tab.gz"

#create index file
tabix -p bed $tissueName".tab.gz" 

