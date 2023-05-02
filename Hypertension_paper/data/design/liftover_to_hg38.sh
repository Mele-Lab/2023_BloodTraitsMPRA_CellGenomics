#!/bin/bash

#SBATCH --job-name=liftover
#SBATCH --workdir=/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/data/design/
#SBATCH --error=liftover.err
#SBATCH --output=liftover.out
#SBATCH --qos=bsc_ls
#SBATCH --time=48:00:00
#SBATCH --nodes=12

/gpfs/projects/bsc83/utils/liftOver snp_coord_master.hg19.simple.chr.sortedplus1.bed /gpfs/projects/bsc83/utils/hg19ToHg38.over.chain snp_coord_master.hg38.simple.chr.sorted.bed Unmapped.liftover.txt