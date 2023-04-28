#get LD matrix ######
library(LDlinkR)
library(ieugwasr)

##### read SNPs for colocalization ######
DBP <- read.table('../../analysis/06_colocalization/DBP_Evangelou_parsed.txt', sep = '\t',
                  header = T)
SBP <- read.table('../../analysis/06_colocalization/SBP_Evangelou_parsed.txt', sep = '\t',
                  header = T)
PP <- read.table('../../analysis/06_colocalization/PP_Evangelou_parsed.txt', sep = '\t',
                  header = T)

#### get LD matrix

sentinel_snps = read.table('../../data/design/sentinel_trait_ld_snps.txt', sep = '\t', header = T)
head(sentinel_snps)
sentinel_snps = sentinel_snps[sentinel_snps$Trait %in% c('DBP','PP','SBP'),]

#### iterate through sentinel snps ##### 
sentinel_analysis <- unique(sentinel_snps$sentinel[sentinel_snps$sentinel %in% c(DBP$SNP, SBP$SNP, PP$SNP)])

for (sentinel in sentinel_analysis) {
  snps_ld <- sentinel_snps$snp[sentinel_snps$sentinel == sentinel]
  snps_ld <- snps_ld[snps_ld %in% c(DBP$SNP, SBP$SNP, PP$SNP)]
  
  if (length(snps_ld)<2) {
    snps_ld <- c(snps_ld,snps_ld)
  }
  
  print(sentinel)
  
  # LD_2 <- LDmatrix(snps = snps_ld, 
  #                    pop = c("CEU","TSI","FIN","GBR","IBS"), 
  #                    token = "c2ed83367cc2",
  #                    genome_build = "grch38")
  
  LD_plink <- ld_matrix(
    snps_ld,
    plink_bin = genetics.binaRies::get_plink_binary(),
    bfile = "Downloads/ieugwasr/EUR"
  )
  
  
  print("done LD")
  
  write.table(LD_plink, paste0('../../analysis/06_colocalization/sentinel_LD/',
                         sentinel,'_LD_matrix.txt'), quote = F, col.names = T, row.names = T, sep = '\t')
  
}


