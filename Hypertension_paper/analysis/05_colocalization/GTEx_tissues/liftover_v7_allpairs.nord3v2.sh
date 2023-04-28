#!/bin/bash
#SBATCH --job-name="liftover all_pairs"
#SBATCH --chdir=/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/06_colocalization/GTEx_tissues
#SBATCH --output=logs/liftover.%A_%a.out
#SBATCH --error=logs/liftover.%A_%a.err
#SBATCH --time=03:00:00
#SBATCH --array=2-46
#SBATCH --qos=bsc_ls
#SBATCH -N 2
#SBATCH -n 2

export tissue_file_all=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../ColocQuiaL-main/GTEx_v8_Tissue_Summary_with_filenames_all.csv | cut -d ',' -f5)
export tissue_file_sig=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../ColocQuiaL-main/GTEx_v8_Tissue_Summary_with_filenames_all.csv | cut -d ',' -f6)

export tissue=$(echo ${tissue_file_sig} | cut -d '.' -f1)

export PATH=/apps/LIFTOVER/linux.x86_64.v369:$PATH   

zcat hg19/${tissue_file_all} | awk '{print "chr"$1"\t"$2-1"\t"$3-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > hg19/${tissue}_temp_hg19.bed

liftOver hg19/${tissue}_temp_hg19.bed /gpfs/projects/bsc83/utils/hg19ToHg38.over.chain hg19/${tissue}_temp_hg38.bed hg19/${tissue}_temp_hg19.unmapped -bedPlus=3 -tab