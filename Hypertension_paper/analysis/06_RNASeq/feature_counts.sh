#!/bin/bash

# User input files ####
infile=$1 # File with fastq filenames ($1) and inpath ($2)  
outdir=$2 # Directory to save fastqc results
indir=$3

# Variables ####
a=($(sed -n "${SLURM_ARRAY_TASK_ID}p" ${infile}))
sample=${a[1]}
#inpath=${a[2]}

inpath=${indir}/${sample}
# Output directory ####
outpath=${outdir}/${sample}
mkdir -p ${outpath}

annotation=/gpfs/projects/bsc83/Data/gene_annotation/gencode/release_26/gencode.v26.basic.annotation.gtf

module load subread

featureCounts -p -a ${annotation} -t exon -g gene_id -o ${outpath}/${sample}.2.counts.txt ${inpath}/${sample}.2.bam
