#!/bin/bash

#SBATCH --job-name=fastqc
#SBATCH --workdir=/gpfs/projects/bsc83/Projects/Breast/Fastq/primeEditing/XT6B10I/DEL20463.20220531/220510_A00481_0412_BHM25YDSX3/
#SBATCH --error=/gpfs/projects/bsc83/Projects/Breast/Fastq/primeEditing/XT6B10I/DEL20463.20220531/fastqc.err
#SBATCH --output=/gpfs/projects/bsc83/Projects/Breast/Fastq/primeEditing/XT6B10I/DEL20463.20220531/fastqc.out
#SBATCH --qos=bsc_ls
#SBATCH --time=10:00:00
#SBATCH --nodes=2


module load java fastqc/0.11.5

fastqc 220510_A00481_0412_BHM25YDSX3/7*