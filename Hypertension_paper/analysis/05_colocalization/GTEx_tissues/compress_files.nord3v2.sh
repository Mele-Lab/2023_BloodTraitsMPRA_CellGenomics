#!/bin/bash
#SBATCH --job-name="compress"
#SBATCH --chdir=/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/06_colocalization/GTEx_tissues
#SBATCH --output=logs/compress.out
#SBATCH --error=logs/compress.err
#SBATCH --time=03:00:00
#SBATCH --qos=bsc_ls
#SBATCH --constraint='medmem'
#SBATCH --cpus-per-task=10


#cp Artery_Coronary.allpairs.tab.gz Artery_Coronary.allpairs.tab
#cp Artery_Aorta.allpairs.tab.gz Artery_Aorta.allpairs.tab

gzip Artery_Aorta.allpairs.tab
gzip Artery_Coronary.allpairs.tab 
