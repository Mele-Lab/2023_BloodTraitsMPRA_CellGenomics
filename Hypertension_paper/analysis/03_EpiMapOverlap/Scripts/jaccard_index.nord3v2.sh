#!/bin/sh

#SBATCH --job-name="bedtools jaccard"
#SBATCH --chdir=/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/03_EpiMapOverlap
#SBATCH --output=jaccard.%A_%a.out
#SBATCH --error=jaccard.%A_%a.err
#SBATCH --time=03:00:00
#SBATCH --array=1-6
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

SEEDFILE='jaccard_comparisons.txt'

export file1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SEEDFILE | cut -f1)
export file2=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SEEDFILE | cut -f2)

data='/gpfs/projects/bsc83/Data/Hypertension/'

zcat ${data}${file1} | sort -k 7n | head -n 50000 | sort -k1,1 -k2,2n > ${data}${file1}.sortF.top.bed
zcat ${data}${file2} | sort -k 7n | head -n 50000 | sort -k1,1 -k2,2n > ${data}${file2}.sortF.top.bed

#zcat ${data}${file1} | sort -k1,1 -k2,2n > ${data}${file1}.sortF.bed
#zcat ${data}${file2} | sort -k1,1 -k2,2n > ${data}${file2}.sortF.bed

bedtools jaccard -a ${data}${file1}.sortF.top.bed -b ${data}${file2}.sortF.top.bed > ${file2}_${file1}_jaccard.top.txt

