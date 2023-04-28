#!/bin/bash
# Runs colocquial.R 
#This version is designed to use an LSF job submission to parallelize the coloc jobs
# Input files needed in directory: qtl_config.sh (modified) 

module load R
module load bedtools

echo "modules loaded"

#source the config file for this run (parameters such as GWAS file location and column names) and the setup config (paths to dependency files for the pipeline)
source /gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/06_colocalization/colocQuial_res/qtl_config.sh
source $setup_config_sh

#replace "TRAITNAME" with the $trait in out and err file names for summary file bsub
sed "s/TRAITNAME/$trait/" $colocquial_dir/ColocQuiaL-main/summarize_results.bsub | sed "s|COLOCQUIAL_DIR|$colocquial_dir/ColocQuiaL-main/|" > ./summarize_results.sh

#run bsub to collect all of the COLOC results into 1 file
sbatch summarize_results.sh 
