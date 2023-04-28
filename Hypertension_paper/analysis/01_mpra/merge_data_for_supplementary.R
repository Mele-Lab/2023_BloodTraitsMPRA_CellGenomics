library(dplyr)
library(tidyr)

table_phil_CM <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/2021_04_30_Results_with_all_info_FINAL.cm.txt', header = T, sep = '\t')
table_phil_VSMC <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/2021_04_30_Results_with_all_info_FINAL.vsmc.txt', header = T, sep = '\t')

head(table_phil_CM)

### read LD-link results ####
ldlink <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/PDLink.results.134Sentinel.txt', 
            sep = '\t', header = T)

### read repeat element ####
repeat_masker <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/snps_repeat_masker.info.txt', sep=' ', header=T)
repeat_masker$query <- gsub('_.*','',repeat_masker$query)
repeat_masker <- repeat_masker %>% distinct()

### read fimo results ####
tile_motifs = read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/tiles_tfbs.best_motifs.csv', sep='\t', header = T)

### merge with LDLink ####
head(ldlink)
table_phil_CM_ld <- merge(table_phil_CM, ldlink[,c('sentinel','RS_Number','Correlated_Alleles','RegulomeDB')], by.x = 'snp_info', by.y = 'RS_Number', all.x = T)
table_phil_VSMC_ld <- merge(table_phil_VSMC, ldlink[,c('sentinel','RS_Number','Correlated_Alleles','RegulomeDB')], by.x = 'snp_info', by.y = 'RS_Number', all.x = T)

### merge with repeat element ####
head(repeat_masker)
vals_significance <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/02_activs/CM_vals.significance.5.txt', sep = '\t', header = T)
head(vals_significance)
vals_significance <- vals_significance %>%
  separate(name, c("pos", "snp_info","info"), "__")
head(vals_significance)
oligo_info <- vals_significance[,c('snp_info','pos')]
colnames(oligo_info) <- c('snp_info','query')

repeat_masker <- merge(repeat_masker, oligo_info, all.x = T)
repeat_masker <- repeat_masker %>% distinct()

table_phil_CM_ld_repeat <- merge(table_phil_CM_ld, repeat_masker[,c('family','beginR','snp_info')], by.x = 'snp_info', by.y = 'snp_info', all.x = T)
table_phil_VSMC_ld_repeat <- merge(table_phil_VSMC_ld, repeat_masker[,c('family','beginR','snp_info')], by.x = 'snp_info', by.y = 'snp_info', all.x = T)

### get TF info --> same TF and TFBS between pairs ####
head(tile_motifs)
tile_motifs <- tile_motifs %>%
  separate(sequence_name, c("pos", "snp_info","info"), "__")
head(tile_motifs)
motifs_wildtype <- tile_motifs[tile_motifs$tile_type == 'WILDTYPE_BUT',]
motifs_snps <- tile_motifs[tile_motifs$tile_type == 'WILDTYPE_SNP',]

snps <- list()
sameTF <- list()
sameTFBS <- list()
#rbind.data.frame(same_tf)
for (snp in unique(tile_motifs$snp_info)) {
  print(snp)
  wildtype_tfbs <- motifs_wildtype[motifs_wildtype$snp_info == snp,]
  snp_tfbs <- motifs_snps[motifs_snps$snp_info == snp,]
  int_tf <- intersect(wildtype_tfbs$motif_alt_id, snp_tfbs$motif_alt_id)
  n_same_tf <- length(int_tf)
  if (n_same_tf >= 1) {
    tfbs_merged <- merge(as.data.frame(table(wildtype_tfbs$motif_alt_id)),as.data.frame(table(snp_tfbs$motif_alt_id)), by='Var1')
    tfbs_merged$diff <-  tfbs_merged$Freq.x - tfbs_merged$Freq.y
    n_same_tfb <- sum(tfbs_merged$Freq.y[tfbs_merged$diff == 0])
  }else{
    n_same_tfb <- 0
  }
  snps <- append(snps, snp)
  sameTF <- append(sameTF, n_same_tf)
  sameTFBS <- append(sameTFBS, list(int_tf))
}
snps_list <- unlist(snps)
names(sameTF) <- snps_list
names(sameTFBS) <- snps_list
results <- list(as.data.frame(unlist(snps)), as.data.frame(unlist(sameTF)))
results_df <- do.call('cbind.data.frame', results)
colnames(results_df) <- c('snp_info','n_same_tfbs')

#### merge with TFBS info ####
table_phil_CM_ld_repeat_tf <- merge(table_phil_CM_ld_repeat, results_df, by.x = 'snp_info', by.y = 'snp_info', all.x = T)
table_phil_VSMC_ld_repeat_tf <- merge(table_phil_VSMC_ld_repeat, results_df, by.x = 'snp_info', by.y = 'snp_info', all.x = T)

### save results #####
write.table(table_phil_CM_ld_repeat_tf, '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/2021_04_30_Results_with_all_info_FINAL.cm.info.txt',
            quote = F, row.names = F, col.names = T, sep = '\t')
write.table(table_phil_VSMC_ld_repeat_tf, '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/2021_04_30_Results_with_all_info_FINAL.vsmc.info.txt',
            quote = F, row.names = F, col.names = T, sep = '\t')

