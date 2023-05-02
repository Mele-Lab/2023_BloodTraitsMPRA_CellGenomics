#!/bin/bash

#SBATCH --job-name=fimo
#SBATCH --workdir=/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/data/03_fimo/
#SBATCH --error=fimo.PE.err
#SBATCH --output=fimo.PE.out
#SBATCH --qos=bsc_ls
#SBATCH --time=4:00:00
#SBATCH --nodes=1

#CPU = 10 
CPU=4

module load meme

#fimo --o FIMO_results/hum_tf Kaia_FIMO/human_curated_tfs.txt /gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/data/Hypertension__pooled.index.fa ### index for MPRA seqs

##### Script to run FIMO on MPRA sequences
fimo --o FIMO_results/hum_tf_PE Kaia_FIMO/human_curated_tfs.txt /gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/luciferase/luciferase_haplotypes_forFIMO.txt ### luciferase seqs