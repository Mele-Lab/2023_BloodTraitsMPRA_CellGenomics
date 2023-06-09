#ID for user's trait of interest. (Can be any string)
trait = "PP" 
#path to the input files
traitFilePath = "/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/06_colocalization/PP_Evangelou_parsed.txt" 
#column IDs from trait file
trait_A1col = "Tested_Allele" 
trait_A2col = "Other_Allele" 
trait_SNPcol = "SNP" 
trait_CHRcol = "CHR" 
trait_BPcol = "POS" 
trait_Pcol = "P" 
trait_Ncol = "N" 
trait_MAFcol = "Freq_Tested_Allele" 
#trait info not in the input file
#traitType is set either to "cc" or "quant"
traitType = "quant" 
#This is the proportion of samples that are cases in a case control GWAS, if you are using a quantitative trait this should be set to "". traitProp = cases / case + controls
traitProp =  #look this up
#locus information for running coloc. 
chrom = 2
colocStart = 216287093
colocStop = 216306635
#reference genome build: "hg19" or "hg38"
build = "hg19"
lead_SNP = "rs1250259" 
#"eqtl" or "sqtl"
qtlType = "eqtl"
#plink parameters
clump_P1 = 0.0000001
clump_KB = 1000
clump_R2 = 0.2
#config file with paths to the qtl data and plink files
#setup_config_R = "/project/voight_GWAS/bychen9/eQTL_colocalizer/setup_config.R"
setup_config_R = "/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/06_colocalization/ColocQuiaL-main/setup_config.R"
