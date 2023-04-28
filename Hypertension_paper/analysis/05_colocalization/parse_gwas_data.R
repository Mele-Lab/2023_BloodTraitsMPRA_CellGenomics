###### modify summary statistics --------

### add rsID -----
### dbp
DBP <- read.table('../../data/Nature_Gen_SNPs/evangelou/Evangelou_30224653_DBP.txt.gz',
                  sep = ' ', header = T)
head(DBP)

sentinel_snps = read.table('../../data/design/sentinel_trait_ld_snps.txt', sep = '\t', header = T)
head(sentinel_snps)
sentinel_snps = sentinel_snps[sentinel_snps$Trait %in% c('DBP','PP','SBP'),]

DBP$change <- paste0(gsub(':SNP','',DBP$MarkerName))
head(DBP)

DBP_filt <- DBP[DBP$change %in% gsub(':[A-Z]+:[A-Z]+$','',sentinel_snps$change[sentinel_snps$Trait == 'DBP']),]
table(sentinel_snps$Trait)

rm(DBP)

### SBP 
SBP <- read.table('../../data/Nature_Gen_SNPs/evangelou/Evangelou_30224653_SBP.txt.gz',
                  sep = ' ', header = T)
head(SBP)

SBP$change <- paste0(gsub(':SNP','',SBP$MarkerName))
head(SBP)

SBP_filt <- SBP[SBP$change %in% gsub(':[A-Z]+:[A-Z]+$','',sentinel_snps$change[sentinel_snps$Trait == 'SBP']),]
table(sentinel_snps$Trait)

rm(SBP)

#### PP 
PP <- read.delim('../../data/Nature_Gen_SNPs/evangelou/Evangelou_30224653_PP.txt.gz',
                  sep = ' ', header = F, skip = 1)
head(PP)

PP$change <- paste0(gsub(':SNP','',PP$V1))
colnames(PP) <- colnames(SBP_filt)
head(PP)

PP_filt <- PP[PP$change %in% gsub(':[A-Z]+:[A-Z]+$','',sentinel_snps$change[sentinel_snps$Trait == 'PP']),]
table(sentinel_snps$Trait)

rm(PP)

#### ad rsid ----
library(dplyr)
library(tidyr)

sentinel_snps$info <- gsub(':[A-Z]+:[A-Z]+$','',sentinel_snps$change)

DBP_filt <-  DBP_filt %>%
  separate(change, c("CHR", "POS"), ":")
head(DBP_filt)

DBP_filt$info <- gsub(':SNP','',DBP_filt$MarkerName)

DBP_filt_m <- merge(DBP_filt, sentinel_snps[,c('info','snp','E')])

SBP_filt <-  SBP_filt %>%
  separate(change, c("CHR", "POS"), ":")
head(SBP_filt)

SBP_filt$info <- gsub(':SNP','',SBP_filt$MarkerName)

SBP_filt_m <- merge(SBP_filt, sentinel_snps[,c('info','snp','E')])

PP_filt <-  PP_filt %>%
  separate(change, c("CHR", "POS"), ":")
head(PP_filt)

PP_filt$info <- gsub(':SNP','',PP_filt$MarkerName)

PP_filt_m <- merge(PP_filt, sentinel_snps[,c('info','snp','E')])

#### reorganize and rename columns -----
DBP_filt_m <-  DBP_filt_m %>%
  separate(E, c("Ref_Allele", "Alt_Allele"), ":")
head(DBP_filt_m)
DBP_filt_m <- DBP_filt_m[,c('CHR','POS','snp','Allele1','Allele2','Freq1','Effect','StdErr','P','N_effective',
                            "Ref_Allele", "Alt_Allele")]
colnames(DBP_filt_m) <- c('CHR','POS','SNP','Tested_Allele','Other_Allele','Freq_Tested_Allele','BETA','SE','P',
                          'N','Ref_Allele','Alt_Allele')
DBP_filt_m$Tested_Allele <- toupper(DBP_filt_m$Tested_Allele)
DBP_filt_m$Other_Allele <- toupper(DBP_filt_m$Other_Allele)

###SBP
SBP_filt_m <-  SBP_filt_m %>%
  separate(E, c("Ref_Allele", "Alt_Allele"), ":")
head(SBP_filt_m)
SBP_filt_m <- SBP_filt_m[,c('CHR','POS','snp','Allele1','Allele2','Freq1','Effect','StdErr','P','N_effective',
                            "Ref_Allele", "Alt_Allele")]
colnames(SBP_filt_m) <- c('CHR','POS','SNP','Tested_Allele','Other_Allele','Freq_Tested_Allele','BETA','SE','P',
                          'N','Ref_Allele','Alt_Allele')
SBP_filt_m$Tested_Allele <- toupper(SBP_filt_m$Tested_Allele)
SBP_filt_m$Other_Allele <- toupper(SBP_filt_m$Other_Allele)

##PP 
PP_filt_m <-  PP_filt_m %>%
  separate(E, c("Ref_Allele", "Alt_Allele"), ":")
head(PP_filt_m)
PP_filt_m <- PP_filt_m[,c('CHR','POS','snp','Allele1','Allele2','Freq1','Effect','StdErr','P','N_effective',
                            "Ref_Allele", "Alt_Allele")]
colnames(PP_filt_m) <- c('CHR','POS','SNP','Tested_Allele','Other_Allele','Freq_Tested_Allele','BETA','SE','P',
                          'N','Ref_Allele','Alt_Allele')
PP_filt_m$Tested_Allele <- toupper(PP_filt_m$Tested_Allele)
PP_filt_m$Other_Allele <- toupper(PP_filt_m$Other_Allele)

write.table(DBP_filt_m, '../../analysis/06_colocalization/DBP_Evangelou_parsed.txt',
            sep = '\t', quote = F, col.names = T, row.names = F)
write.table(SBP_filt_m, '../../analysis/06_colocalization/SBP_Evangelou_parsed.txt',
            sep = '\t', quote = F, col.names = T, row.names = F)
write.table(PP_filt_m, '../../analysis/06_colocalization/PP_Evangelou_parsed.txt',
            sep = '\t', quote = F, col.names = T, row.names = F)


#### change beta values according to ref and alternative allele -----
DBP <- read.delim('../../analysis/06_colocalization/DBP_Evangelou_parsed.txt')
SBP <- read.delim('../../analysis/06_colocalization/SBP_Evangelou_parsed.txt')
PP <- read.delim('../../analysis/06_colocalization/PP_Evangelou_parsed.txt')

DBP$BETA[DBP$Tested_Allele == DBP$Ref_Allele] <- -(DBP$BETA[DBP$Tested_Allele == DBP$Ref_Allele])
SBP$BETA[SBP$Tested_Allele == SBP$Ref_Allele] <- -(SBP$BETA[SBP$Tested_Allele == SBP$Ref_Allele])
PP$BETA[PP$Tested_Allele == PP$Ref_Allele] <- -(PP$BETA[PP$Tested_Allele == PP$Ref_Allele])

write.table(DBP, '../../analysis/06_colocalization/DBP_Evangelou_parsed.txt',
            sep = '\t', quote = F, col.names = T, row.names = F)
write.table(SBP, '../../analysis/06_colocalization/SBP_Evangelou_parsed.txt',
            sep = '\t', quote = F, col.names = T, row.names = F)
write.table(PP, '../../analysis/06_colocalization/PP_Evangelou_parsed.txt',
            sep = '\t', quote = F, col.names = T, row.names = F)


#### correlated GWAS with MPRA data #####
### read differential results ####
CM_results <- read.table('../../analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_with_significance.CM.new_back.005.txt', sep='\t', header = T)
VSMC_results <- read.table('../../analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_with_significance.VSMC.new_back.005.txt', sep='\t', header = T)
### index 
index <- read.table('../../data/design/Hypertension__pooled.index.txt', sep = '\t', header = T)
head(index)

#### merge results with index ####
index['dupe_info'] <- gsub('\\..*','',index$tile_id)

CM_results = merge(CM_results,index[index$tile_type != 'RANDOM',c("name","dupe_info")], by="dupe_info")
VSMC_results = merge(VSMC_results,index[index$tile_type != 'RANDOM',c("name","dupe_info")], by="dupe_info")

CM_results <- CM_results %>% distinct()
VSMC_results <- VSMC_results %>% distinct()

CM_results
VSMC_results

CM_results <- CM_results %>%
  separate(name, c("pos", "snp_info","info"), "__")
head(CM_results)
VSMC_results <- VSMC_results %>%
  separate(name, c("pos", "snp_info","info"), "__")
head(VSMC_results)

DBP <- read.table('../../analysis/06_colocalization/DBP_Evangelou_parsed.txt',
                  sep = '\t', header = T)
SBP <- read.table('../../analysis/06_colocalization/SBP_Evangelou_parsed.txt',
                  sep = '\t', header = T)
PP <- read.table('../../analysis/06_colocalization/PP_Evangelou_parsed.txt',
                  sep = '\t', header = T)
head(DBP)

all_gwas <- rbind(DBP,SBP,PP)

snps_gwas <- merge(VSMC_results, all_gwas, by.x='snp_info',by.y='SNP')
snps_gwas$P <- -log10(as.numeric(snps_gwas$P))
snps_gwas$fdr_comp <- -log10(as.numeric(snps_gwas$fdr_comp))
snps_gwas$logFC_comp <- abs(as.numeric(snps_gwas$logFC_comp))
snps_gwas$BETA <- abs(as.numeric(snps_gwas$BETA))

library("ggpubr")
ggscatter(snps_gwas, x = "fdr_comp", y = "P", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "MPRA abs(logFC)", ylab = "GWAS abs(BETA)")


## how many sentinel are reg? ####

CM_results <- read.table('../../analysis/CM_eqtls_rep_cons_dhs_tbs_scores.txt',
                         sep='\t', header = T)
VSMC_results <- read.table('../../analysis/VSMC_eqtls_rep_cons_dhs_tbs_scores.txt', 
                           sep='\t', header = T)

sentinel_reg <- nrow(VSMC_results[VSMC_results$snp_info == VSMC_results$sentinel,])
non_sent_reg <- nrow(VSMC_results[VSMC_results$snp_info != VSMC_results$sentinel,])
dat <- data.frame(
  "Reg" = c(sentinel_reg, non_sent_reg),
  "non Reg" = c(135-sentinel_reg, 4475-non_sent_reg),
  row.names = c("sentinel", "LD"),
  stringsAsFactors = FALSE
)

dat

fisher.test(dat)  
