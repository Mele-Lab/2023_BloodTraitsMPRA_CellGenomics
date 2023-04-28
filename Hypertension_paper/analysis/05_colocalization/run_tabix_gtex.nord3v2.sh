#!/bin/bash
#SBATCH --job-name="create_tabix"
#SBATCH --chdir=/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/06_colocalization/GTEx_tissues
#SBATCH --output=logs/create_tabix.sig.%A_%a.out
#SBATCH --error=logs/create_tabix.sig.%A_%a.err
#SBATCH --time=03:00:00
#SBATCH --array=2-46
#SBATCH --qos=bsc_ls
#SBATCH -N 2
#SBATCH -n 2


#export tissue=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../GTEx_all_pairs_file.txt| cut -f1)
export tissue_file_all=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../ColocQuiaL-main/GTEx_v8_Tissue_Summary_with_filenames_all.csv | cut -d ',' -f5)
export tissue_file_sig=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../ColocQuiaL-main/GTEx_v8_Tissue_Summary_with_filenames_all.csv | cut -d ',' -f6)

export tissue=$(echo ${tissue_file_sig} | cut -d '.' -f1)

#ln -s /gpfs/projects/bsc83/Data/GTEx/v8/cis_QTLs/cis_eQTLs/GTEx_Analysis_v8_eQTL/${tissue}.signif_variant_gene_pairs.txt.gz GTEx_tissues/

echo ${tissue}

bash ../ColocQuiaL-main/create_tabix_gtex_eqtl_allpairs.sh /gpfs/projects/bsc83/Data/GTEx/v7/cis_QTLs/cis_eQTLs/GTEx_Analysis_v7_eQTL/${tissue}.allpairs.txt.gz

