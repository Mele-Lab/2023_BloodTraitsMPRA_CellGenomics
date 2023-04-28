#!/bin/sh

#SBATCH --job-name="bedtools intersect"
#SBATCH --chdir=/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/03_EpiMapOverlap
#SBATCH --output=bedtools.%A_%a.out
#SBATCH --error=bedtools.%A_%a.err
#SBATCH --time=03:00:00
#SBATCH --array=1
#SBATCH --qos=bsc_ls
#SBATCH --constraint='medmem'
#SBATCH -N 1


module load bedtools
#CPU = 10 
CPU=4
#set -e

#SLURM_ARRAY_TASK_ID=1
echo $LSB_JOBINDEX
########## ARRAY STUFF ##########################################
#SLURM_ARRAY_TASK_ID=1

SEEDFILE='files_dhs_all.txt'

export file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SEEDFILE | cut -f1)
export sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SEEDFILE | cut -f2)

data='/gpfs/projects/bsc83/Data/Hypertension/EpiMap/'
#SNPs='/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/data/design/snp_coord_master.hg19.simple.chr.sorted.bed'
SNPs='/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/data/design/sentinel_coordinates_block.bed'

#zcat ${data}${file} | tail -n +2 | awk '{print $1"\t"$2"\t"$3"\t"$4}' | sort -k1,1 -k2,2n -k3,3n  > ${data}${file}.sort.bed
file='epimap_dhs_summit_200bp_hg19.sort.bed'
sort -k1,1V -k2,2n -k3,3n $SNPs > ${SNPs}_EpiMap.bed
sort -k1,1 -k2,2n -k3,3n ${data}${file} > ${data}${file}.sortF.bed

bedtools intersect -wa -wb \
    -a ${data}${file}.sortF.bed \
    -b ${SNPs}_EpiMap.bed \
    -sorted > ${sample}_loci_overlap.bed

