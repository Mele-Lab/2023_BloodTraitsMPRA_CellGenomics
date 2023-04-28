#!/bin/bash

# User input files ####
infile=$1 # File with fastq filenames ($1) and inpath ($2)  
outdir=$2 # Directory to save mapping results
indir=$3 # Directory with fastq files

# Variables ####
a=($(sed -n "${SLURM_ARRAY_TASK_ID}p" ${infile}))
sample=${a[0]}
#inpath=${a[2]}


inpath=${indir}
# Output directory ####
outpath=${outdir}/${sample}/
mkdir -p ${outpath}

sample_name_L1=${sample}_R1.fastq.gz
sample_name_L2=${sample}_R2.fastq.gz

assembly=/gpfs/projects/bsc83/Data/hisat2_indeces/GRCh38/GRCh38

# Output directory ####
#novel_ss=${outpath}/${sample_name}.novel_ss.txt
alin_summary=${outpath}/${sample}.hisat2_summary.txt
mkdir -p ${outpath}

# Map PE RNAseq reads
module load gcc/7.2.0 hisat2 samtools 

/gpfs/projects/bsc83/utils/hisat2-2.1.0/hisat2 -x ${assembly} -1 ${inpath}/${sample_name_L1} -2 ${inpath}/${sample_name_L2} --phred33  --time --summary-file ${alin_summary}  | samtools sort -T ${TMPDIR}/${sample}.tmp  -O bam -o ${outpath}/${sample}.bam -


samtools index ${outpath}/${sample}.bam


