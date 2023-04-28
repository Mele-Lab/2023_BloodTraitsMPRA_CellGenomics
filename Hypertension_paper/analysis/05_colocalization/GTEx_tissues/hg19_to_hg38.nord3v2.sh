#!/bin/bash
#SBATCH --job-name="hg19 to hg38 all_pairs"
#SBATCH --chdir=/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/06_colocalization/GTEx_tissues
#SBATCH --output=logs/hg19_to38.%A_%a.out
#SBATCH --error=logs/hg19_to38.%A_%a.err
#SBATCH --time=03:00:00
#SBATCH --array=2-46
#SBATCH --qos=bsc_ls
#SBATCH --constraint='highmem'
#SBATCH --cpus-per-task=10

export tissue_file_all=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../ColocQuiaL-main/GTEx_v8_Tissue_Summary_with_filenames_all.csv | cut -d ',' -f5)
export tissue_file_sig=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../ColocQuiaL-main/GTEx_v8_Tissue_Summary_with_filenames_all.csv | cut -d ',' -f6)

export tissue=$(echo ${tissue_file_sig} | cut -d '.' -f1)

module load R 

Rscript hg19_to_hg38_all_pairs.R ${tissue}

cat ${tissue}.allpairs.tab | vcf-sort -p 8  | bgzip -c > ${tissue}.allpairs.tab.gz

#gzip ${tissue}.allpairs.tab