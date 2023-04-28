#!/bin/bash
#SBATCH --job-name="summarize_qtl_results"
#SBATCH --chdir=/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/06_colocalization/
#SBATCH --output=logs/TRAITNAME_COLOC_Results_Summary.%A_%a.out
#SBATCH --error=logs/TRAITNAME_COLOC_Results_Summary.%A_%a.err
#SBATCH --time=10:00:00
#SBATCH --qos=bsc_ls
#SBATCH -N 2
#SBATCH -n 2

#Run cod to collect the QTL colocalization results after they have all been generated.


module load R

COLOCQUIAL_DIR/summarize_qtl_results.sh
