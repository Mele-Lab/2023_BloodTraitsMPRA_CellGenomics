 #!/bin/sh

#BSUB -cwd /gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/03_EpiMapOverlap
#BSUB -J bedtools[1]%100
#BSUB -e bedtools.%J_%I.err
#BSUB -o bedtools.%J_%I.out
#BSUB -q bsc_ls
#BSUB -W 24:00
#BSUB -n 2
#BSUB -M 7000
#tot 378


module purge && module load java/1.8.0u66 intel/2017.4 impi/2017.4 MKL/2017.4 gcc/7.2.0 OPENSSL/1.1.1c PYTHON/3.7.4_pip SAMTOOLS/1.9 subread GATK/4.1.2.0 BEDTOOLS/2.25.0 R/3.6.3

#CPU = 10 
CPU=4
#set -e

file='TSS_map_hg10_gencode_release34.bed'
sample='closest_TSS'
data='/gpfs/projects/bsc83/Data/Hypertension/EpiMap/'
SNPs='/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/data/design/snp_coord_master.hg19.simple.chr.sorted.bed'

#zcat ${data}${file} | tail -n +2 | sort -k1,1V -k2,2n -k3,3n  > ${data}${file}.sort.bed
#sort -k1,1V -k2,2n -k3,3n $SNPs > ${SNPs}_EpiMap.bed
sort -k1,1 -k2,2n -k3,3n ${data}${file} > ${data}${file}.sort.bed

bedtools closest -b ${data}${file}.sort.bed -a ${SNPs} -d > ${sample}_SNPS_overlap.bed
