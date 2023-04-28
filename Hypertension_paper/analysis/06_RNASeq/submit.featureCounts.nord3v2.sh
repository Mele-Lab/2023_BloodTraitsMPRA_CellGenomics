#!/bin/bash
#SBATCH --job-name="hisat2"
#SBATCH --chdir=/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/RNASeq/
#SBATCH --output=logs/feature_counts.PE.%A_%a.err
#SBATCH --error=logs/feature_counts.PE.%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --array=21-22
#SBATCH --qos=bsc_ls
#SBATCH -N 2
#SBATCH -n 2

/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/RNASeq/feature_counts.sh PE.samples.txt feature_counts hisat2/
