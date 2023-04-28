#!/bin/bash

#SBATCH --job-name=hisat2
#SBATCH --workdir=/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/RNASeq/
#SBATCH --array=21-22
#SBATCH --error=/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/RNASeq/logs/hisat2.PE.%A_%a.err
#SBATCH --output=/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/RNASeq/logs/hisat2.PE.%A_%a.out
#SBATCH --qos=bsc_ls
#SBATCH --time=03:00:00
#SBATCH --nodes=2

/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/RNASeq/hisat2.PE.sh PE.samples.txt hisat2 \
trimmed_reads/less/