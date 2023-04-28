#!/bin/sh
#SBATCH --job-name=cutadapt
#SBATCH --workdir=/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/RNASeq/
#SBATCH --cpus-per-task=10   
#SBATCH --ntasks=1
#SBATCH --qos=bsc_ls
#SBATCH --time=10:00:00
#SBATCH --output=/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/RNASeq/logs/cutadap.PE.%A_%a.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/RNASeq/logs/cutadapt.PE.%A_%a.err
#SBATCH --array=21-22
# normal queue bsc_ls

#tot 378

#CPU = 10 
CPU=4
#set -e

#SLURM_ARRAY_TASK_ID=1
echo $SLURM_ARRAY_TASK_ID
########## ARRAY STUFF ##########################################
#SLURM_ARRAY_TASK_ID=1

SEEDFILE='PE.samples.txt'

a=($(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SEEDFILE}))
sample=${a[0]}
name=${a[1]}

sample_name_L1=${sample}_L001_R1_001.fastq.gz
sample_name_L2=${sample}_L001_R2_001.fastq.gz
#sample_name_L1_2=${sample}_L001_R1_001.fastq.gz
#sample_name_L2_2=${sample}_L001_R2_001.fastq.gz

#################################################################3
inpath='/gpfs/projects/bsc83/Projects/Breast/Fastq/8VQ67JA/DEL19783.20220303/220228_A00481_0367_BH7MY5DSX3/'
inpath='/gpfs/projects/bsc83/Projects/Breast/Fastq/primeEditing/XT6B10I/DEL20463.20220531/220510_A00481_0412_BHM25YDSX3/'

OUTDIR='trimmed_reads/'less

mkdir -p $OUTDIR

module load python/3.6.1

/gpfs/projects/bsc83/Projects/Breast/barcode_counting_materials/programs/trim_galore_zip/trim_galore -q 20 -o ${OUTDIR} --paired ${inpath}${sample_name_L1} ${inpath}${sample_name_L2}

