#!/bin/bash

#SBATCH -D /gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/data/design/
#SBATCH --job-name="chip-atlas intersect"
#SBATCH --error=bedtools.err
#SBATCH --output=bedtools.out
#SBATCH --qos=debug
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH -n 2
#SBATCH --mail-type=end
#SBATCH --mail-user=winona.oliveros@bsc.es

#tot 378


module purge && module load impi/2017.4 intel/2017.4 bedtools/2.25.0

Hypertension='/gpfs/projects/bsc83/Data/Hypertension/'

#sort -k1,1 -k2,2n snp_coord_master.hg19.simple.chr.sorted_1based.bed > snp_coord_master.hg19.simple.chr_sorted_1based.bed
sort -k1,1 -k2,2n ${Hypertension}Oth.Neu.05.AllAg.AllCell.bed > ${Hypertension}Oth.Neu.05.AllAg.AllCell.sorted.bed


bedtools intersect -wa -wb \
    -a snp_coord_master.hg19.simple.chr_sorted_1based.bed \
    -b ${Hypertension}Oth.Neu.05.AllAg.AllCell.sorted.bed \
    -names Chip_Atlas \
    -sorted > overlap_chip_atlas.new.Neuronal.bed