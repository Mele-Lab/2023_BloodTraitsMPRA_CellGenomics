#!/bin/bash

#SBATCH --job-name=bedops
#SBATCH --workdir=/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/data/design/
#SBATCH --error=bedops.err
#SBATCH --output=bedops.out
#SBATCH --qos=bsc_ls
#SBATCH --time=48:00:00
#SBATCH --nodes=12

module load gcc/7.1.0 bedops

data='/gpfs/projects/bsc83/Data/gene_annotation/gencode/release_34/'

closest-features --closest --dist snp_coord_master.hg19.simple.chr.sorted.bed ${data}genes_sorted.hg19.bed > clostest_gene_snps.txt 