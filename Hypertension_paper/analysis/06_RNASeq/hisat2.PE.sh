#!/bin/bash

# User input files ####
infile=$1 # File with fastq filenames ($1) and inpath ($2)  
outdir=$2 # Directory to save mapping results
indir=$3 # Directory with fastq files

# Variables ####
a=($(sed -n "${SLURM_ARRAY_TASK_ID}p" ${infile}))
sample=${a[0]}
name=${a[1]}

echo ${sample}
echo ${name}


inpath=${indir}
# Output directory ####
outpath=${outdir}/${name}/
mkdir -p ${outpath}

sample_name_L1=${sample}_L002_R1_001_val_1.fq.gz
sample_name_L2=${sample}_L002_R2_001_val_2.fq.gz

#sample_name_L1_2=${sample}_L001_R1_001_val_1.fq.gz
#sample_name_L2_2=${sample}_L001_R2_001_val_2.fq.gz

assembly=/gpfs/projects/bsc83/Data/hisat2_indeces/GRCh38/GRCh38

# Output directory ####
alin_summary=${outpath}/${name}.2.hisat2_summary.txt
mkdir -p ${outpath}

# Map PE RNAseq reads
#module load gcc/7.2.0 hisat2 samtools 
module load htslib samtools

#if [ ${stranded} == "Yes" ]
#then
/gpfs/projects/bsc83/utils/hisat2-2.1.0/hisat2 -x ${assembly} -1 ${inpath}/${sample_name_L1} -2 ${inpath}/${sample_name_L2} --phred33  --time --summary-file ${alin_summary}  | samtools sort -T ${TMPDIR}/${name}.tmp  -O bam -o ${outpath}/${name}.2.bam -
#else
#	hisat2-align-s -x ${assembly} -1  ${inpath}/${sample}.1.fastq -2 ${inpath}/${sample}.2.fastq --phred33  --time --summary-file ${alin_summary} | samtools sort -T ${TMPDIR}/${sample_name}.tmp -O bam -o ${outpath}/${sample_name}.bam -
#fi
#module unload gcc/7.2.0
#module load gcc/7.1.0 samtools/1.5

samtools index ${outpath}/${name}.2.bam


