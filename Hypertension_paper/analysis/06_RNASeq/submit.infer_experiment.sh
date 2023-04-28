#!/bin/bash

#SBATCH --job-name=infer_experiment
#SBATCH --workdir=/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/RNASeq/
#SBATCH --array=21-22
#SBATCH --error=logs/infer_experiment.%A_%a.err
#SBATCH --output=logs/infer_experiment.%A_%a.out
#SBATCH --qos=bsc_ls
#SBATCH --time=02:00:00
#SBATCH --nodes=1


/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/RNASeq/infer_experiment.sh PE.samples.txt hisat2/
