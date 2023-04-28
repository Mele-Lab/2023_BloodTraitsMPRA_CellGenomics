#!/bin/bash
#SBATCH --job-name="summarize_qtl_results"
#SBATCH --output=PP_COLOC_Results_Summary.%A_%a.out
#SBATCH --error=PP_COLOC_Results_Summary.%A_%a.err
#SBATCH --time=10:00:00
#SBATCH --qos=bsc_ls
#SBATCH -N 1
#SBATCH --cpus-per-task=10

module load R

/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/06_colocalization//ColocQuiaL-main//summarize_qtl_results.sh

#/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/06_colocalization//ColocQuiaL-main//summarize_qtl_results.susie.sh
