#!/bin/bash

# User input files ####
infile=$1 # File with fastq filenames ($1) and inpath ($2)  
outdir=$2 # Directory to save fastqc results

# Variables ####
a=($(sed -n "${SLURM_ARRAY_TASK_ID}p" ${infile}))
sample=${a[1]}
#inpath=${a[2]}


gene_annotation=/gpfs/projects/bsc83/Data/gene_annotation/RSeQC/hg38_Gencode_V28.bed


# Output directory ####
outpath=${outdir}/${sample}
#mkdir -p ${outpath}

# Infer strandness
module load python/3.6.1
infer_experiment.py -r ${gene_annotation}  -s 1000000 -i ${outpath}/${sample}.bam > ${outpath}/${sample}.inferred_strandness.txt
