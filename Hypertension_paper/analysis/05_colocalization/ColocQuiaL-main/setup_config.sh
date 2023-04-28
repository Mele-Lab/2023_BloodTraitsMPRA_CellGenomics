#this bash config file allows you to set the paths to all the dependency file you need to run the qtl colocalizer pipeline
#this config file is required for running the colocalizer pipeline on multiple loci at a time

#path to the directory where the ColocQuiaL code is saved locally
colocquial_dir="/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/06_colocalization/"

#provide the path to plink refernce files to be used for plink commands
plink_bfile=""

#provide the path to the plink ped files for the list of individuals you wish to use in your LD reference panel
plink_keep=""

#provide the ID of the bsub queue you wish to submit your ColocQuiaL jobs to
bsub_queue="bsc_ls"
