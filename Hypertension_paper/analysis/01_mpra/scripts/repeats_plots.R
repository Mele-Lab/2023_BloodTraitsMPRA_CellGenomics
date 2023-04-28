##### repetitive activity and logFC ------
library(dplyr)
library(tidyr)
##### data ------
### enrichment of repeat families ####
### RepeatMasker ####
sentinel_snps = read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/sentinel_trait_ld_snps.txt', sep = '\t', header = T)
head(sentinel_snps)
sentinel_snps = sentinel_snps[sentinel_snps$Trait %in% c('DBP','PP','SBP'),]

repeat_masker <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/snps_repeat_masker.info.txt', sep=' ', header=T)
repeat_masker$query <- gsub('_.*','',repeat_masker$query)
repeat_masker <- repeat_masker %>% distinct()

sentinel = merge(sentinel_snps,repeat_masker, by.x= 'coord', by.y= 'query')

vals_significance_cm <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/02_activs/CM_vals.significance.5.txt', sep = '\t', header = T)
head(vals_significance_cm)
vals_significance_cm <- vals_significance_cm %>%
  separate(name, c("pos", "snp_info","info"), "__")
head(vals_significance_cm)

vals_significance_vsmc <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/02_activs/VSMC_vals.significance.4.txt', sep = '\t', header = T)
head(vals_significance_vsmc)
vals_significance_vsmc <- vals_significance_vsmc %>%
  separate(name, c("pos", "snp_info","info"), "__")
head(vals_significance_vsmc)

CM_results <- CM_results %>%
  separate(name, c("pos", "snp_info","info"), "__")
head(CM_results)
VSMC_results <- VSMC_results %>%
  separate(name, c("pos", "snp_info","info"), "__")
head(VSMC_results)

###### merge activity with repetitive -----
vals_significance_cm$Repetitive <- 'No Repeat'
vals_significance_cm$Repetitive[vals_significance_cm$snp_info %in% sentinel$snp] <- 'Repetitive'

vals_significance_vsmc$Repetitive <- 'No Repeat'
vals_significance_vsmc$Repetitive[vals_significance_vsmc$snp_info %in% sentinel$snp] <- 'Repetitive'

p <- ggboxplot(vals_significance_vsmc, x = "sig", y = "VSMC_log",
               color = "Repetitive", palette = "jco")
p + stat_compare_means(aes(group = Repetitive) , label = "p.format") + xlab('') + ylab('MPRA activity')
  
##### merge differential activity vs repetitive ------
CM_results$Repetitive <- 'No Repeat'
CM_results$Repetitive[CM_results$snp %in% sentinel$snp] <- 'Repetitive'

VSMC_results$Repetitive <- 'No Repeat'
VSMC_results$Repetitive[VSMC_results$snp %in% sentinel$snp] <- 'Repetitive'

p <- ggboxplot(VSMC_results, x = "Repetitive", y = "abs_ef",
               color = "Repetitive", palette = "jco")
p + stat_compare_means(aes(group = Repetitive) , label = "p.format") + xlab('') + ylab('abs(logFC)')

mu <- ddply(VSMC_results, "Repetitive", summarise, grp.mean=mean(abs_ef, na.rm = TRUE))
head(mu)
