#!/bin/bash

#SBATCH --job-name=fimo
#SBATCH --workdir=/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/data/03_fimo/
#SBATCH --error=logos/fimo.%J_%I.err
#SBATCH --output=logos/fimo.%J_%I.out
#SBATCH --qos=bsc_ls
#SBATCH --array=23-49
#SBATCH --time=48:00:00
#SBATCH --nodes=12


module load meme

SEEDFILE='/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/data/03_fimo/Tf_interested_disrupted.txt'

export a=$(awk "NR==$SLURM_ARRAY_TASK_ID" $SEEDFILE)
IFS="	"
set $a
sample=$1

ceqlogo -i Kaia_FIMO/human_curated_tfs.txt -m $sample -o logos/new/$sample.eps -f EPS