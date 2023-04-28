#this R config file allows you to set the paths to all the dependency file you need to run the ColocQuiaL pipeline
#this config file is required for running the colocalizer pipeline on one locus

#provide the path to plink reference files to be used for plink commands
plink_bfile = ""
#provide the path to the plink ped files for the list of individuals you wish to use in your LD reference panel
plink_keep = ""
#provide path to GRCh37 to GRCh38 variant hash table directory
hash_table_dir = "/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/data/Coloq/"
#LD matrix path
path_ld = '/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/06_colocalization/sentinel_LD/'

#provide path to eqtl tissue table
eQTL_tissue_table = "/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/06_colocalization/ColocQuiaL-main/GTEx_v8_Tissue_Summary_with_filenames_all.csv"

#provide path to significant eqtl data tabix directory NOTE: don't include last slash
eQTL_sig_qtl_tabix_dir = "/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/06_colocalization/GTEx_tissues/"
#column number of column in significant eqtl files that contains geneID
eQTL_sig_geneID_col = 7 

#provide path to the all eqtl data tabix directory
eQTL_all_qtl_tabix_dir ="/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/06_colocalization/GTEx_tissues/"

#provide the header to the eqtl data
eQTL_all_header = c("chrom_b38", "chromStart_b38", "chromEnd_b38", "eGeneID", "A1_eqtl", "A2_eqtl", "build", "tss_distance", "ma_samples", "ma_count", "maf", "pvalue_eQTL", "slope", "slope_se")
#name of column in header representing geneID
eQTL_all_geneID = "eGeneID" 
#name of column in header representing chromosome
eQTL_all_chrom = "chrom_b38"
#name of column in header representing end coordinate 
eQTL_all_chromEnd = "chromEnd_b38"
#name of column in header representing p-value
eQTL_all_pvalue = "pvalue_eQTL"

#intron information
sQTL_all_intron_chr = "intron_chr" 
sQTL_all_intron_bp_first = "intron_bp_first"
sQTL_all_intron_bp_end = "intron_bp_end"
sQTL_all_intron_clu = "intron_clu"


#path to GRCh37 to GRCh38 liftOver chain file
liftOver_chain = "/gpfs/projects/bsc83/utils/hg19ToHg38.over.chain"
#path to recombination rate directory (expected to be 1K genome chromosome recombination rate files or in the same format as these files: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130507_omni_recombination_rates/
recomb_rate_data="/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/06_colocalization/CEU/CEU"
