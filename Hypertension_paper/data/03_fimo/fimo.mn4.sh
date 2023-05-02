#!/bin/bash

#SBATCH --job-name=fimo
#SBATCH --workdir=/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/data/
#SBATCH --error=fimo.err
#SBATCH --output=fimo.out
#SBATCH --qos=bsc_ls
#SBATCH --time=48:00:00
#SBATCH --nodes=12

#CPU = 10 
CPU=4

module load meme

fimo --o FIMO_results Kaia_FIMO/human_cisbp_pfms.txt Hypertension__pooled.index.fa