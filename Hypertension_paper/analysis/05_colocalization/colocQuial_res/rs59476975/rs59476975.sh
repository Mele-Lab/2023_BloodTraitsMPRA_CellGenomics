#!/bin/bash
#SBATCH --job-name="ColocQuiaL"
#SBATCH --output=logs/rs59476975_DBP_COLOC.%A_%a.out
#SBATCH --error=logs/rs59476975_DBP_COLOC.%A_%a.err
#SBATCH --time=10:00:00
#SBATCH --qos=bsc_ls
#SBATCH -N 1
#SBATCH --cpus-per-task=10

export PATH=/apps/LIFTOVER/linux.x86_64.v369:$PATH   

module load R
module load htslib
module load vcftools

export LD_LIBRARY_PATH=/gpfs/apps/NORD3/R/4.1.0/INTEL/lib64/R/lib:$LD_LIBRARY_PATH

Rscript coloc_analysis.R