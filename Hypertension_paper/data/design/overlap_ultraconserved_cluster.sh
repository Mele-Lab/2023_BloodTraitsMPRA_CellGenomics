#!/bin/bash

#SBATCH -D /gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/data/design/
#SBATCH --job-name="phastcons intersect"
#SBATCH --error=bedtools.err
#SBATCH --output=bedtools.out
#SBATCH --qos=bsc_ls
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH -n 2
#SBATCH --mail-type=end
#SBATCH --mail-user=winona.oliveros@bsc.es

#tot 378


module purge && module load impi/2017.4 intel/2017.4 bedtools/2.25.0

fantom='/gpfs/projects/bsc83/Data/PhyloP/'

sort -k1,1 -k2,2n snp_coord_master.hg19.simple.chr.sorted_1based.bed > snp_coord_master.hg19.simple.chr_sorted_1based.bed

bedtools intersect -wa -wb \
    -a snp_coord_master.hg19.simple.chr_sorted_1based.bed \
    -b ${fantom}ultraconserved_element_phastcons_filt_2.sorted.bed \
    -names ultraconserved \
    -sorted > overlap_ultraconserved.bed