## Compare luciferase results with MPRA ----
library(dplyr)
library(plyr)
library(tidyr)
library(ggplot2)

#### read data -----
### read differential results ####
CM_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_with_significance.CM.new_back.005.txt', sep='\t', header = T)
VSMC_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_with_significance.VSMC.new_back.005.txt', sep='\t', header = T)

### index 
index <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/Hypertension__pooled.index.txt', sep = '\t', header = T)
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

### luciferase snps ----
luciferase <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/luciferase_snps.txt',
                         sep = '\t', header = F)
colnames(luciferase) <- c('snp_info','sentinel','gene')

##### merge info ----
luci_stats <- merge(luciferase, CM_results[,c('fdr_comp','logFC_comp','snp_info')])
luci_stats <- merge(luci_stats, VSMC_results[,c('fdr_comp','logFC_comp','snp_info')], by = 'snp_info')

### swipe logfc of reversed snps -----
luci_stats[luci_stats$snp_info == 'rs319686',c('logFC_comp.x','logFC_comp.y')] <- -luci_stats[luci_stats$snp_info == 'rs319686',c('logFC_comp.x','logFC_comp.y')]
luci_stats[luci_stats$snp_info == 'rs6773208',c('logFC_comp.x','logFC_comp.y')] <- -luci_stats[luci_stats$snp_info == 'rs319686',c('logFC_comp.x','logFC_comp.y')]
luci_stats[luci_stats$snp_info == 'rs11712464',c('logFC_comp.x','logFC_comp.y')] <- -luci_stats[luci_stats$snp_info == 'rs319686',c('logFC_comp.x','logFC_comp.y')]
luci_stats[luci_stats$snp_info == 'rs12330363',c('logFC_comp.x','logFC_comp.y')] <- -luci_stats[luci_stats$snp_info == 'rs319686',c('logFC_comp.x','logFC_comp.y')]

luci_stats <- luci_stats[order(luci_stats$gene),]
colnames(luci_stats) <- c('snp_info','sentinel','gene', 'CM_fdr_comp','CM_logFC_comp', 'VSMC_fdr_comp','VSMC_logFC_comp')
library(reshape2)
luci_stats_fdr <- luci_stats[,c('snp_info','sentinel','gene', 'CM_fdr_comp','VSMC_fdr_comp')]
luci_stats_logfc <- luci_stats[,c('snp_info','sentinel','gene', 'CM_logFC_comp','VSMC_logFC_comp')]

colnames(luci_stats_fdr) <- c('snp_info','sentinel','gene','CM','VSMC')
colnames(luci_stats_logfc) <- c('snp_info','sentinel','gene','CM','VSMC')
head(luci_stats_logfc)

#### read luci logfc ----
luci_fc_results <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/luciferase/FC_RefSeq-vs_AltSeq_Wini.txt',
                              header = T, sep = '\t')
luci_fc_results$id <- rownames(luci_fc_results)
luci_fc_results_m <- melt(luci_fc_results,id.vars=c("id"))

luci_fc_results_m <- luci_fc_results_m[grep('5\\.3',luci_fc_results_m$variable),]

### merge luci with mpra -----
sumed <- luci_stats_logfc[,c('gene','CM','VSMC')] %>% 
  dplyr::group_by(gene) %>% 
  dplyr::summarise(CM = sum(CM), VSMC = sum(VSMC))
sumed

mean_val <- luci_stats_logfc %>% 
  dplyr::group_by(gene) %>% 
  dplyr::summarise(CM = median(CM), VSMC = median(VSMC))
mean_val

luci_fc_results_m$gene <- gsub('\\..*','',luci_fc_results_m$variable)
luci_fc_results_m$value <- gsub(',','\\.',luci_fc_results_m$value)
luci_fc_results_m$value <- as.numeric(luci_fc_results_m$value)

luci_fc_results_mean <- luci_fc_results_m %>% 
  dplyr::group_by(gene) %>% 
  dplyr::summarise(luci = mean(value))
luci_fc_results_mean
luci_fc_results_mean$logFC <- log2(luci_fc_results_mean$luci)

all <- merge(mean_val, luci_fc_results_mean)

#### plot correlation ----
library("ggpubr")
ggscatter(all, x = "CM", y = "logFC", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "MPRA CMs", ylab = "Luciferase logFC")

rownames(all) <- all$gene
range <- max(abs(as.matrix(all[,c(2,3,5)])))
pheatmap::pheatmap(as.matrix(all[,c(2,3,5)]), scale = 'column')

#### merge with significance MPRA -----
colnames(luci_stats_fdr) <- c('snp_info','sentinel','gene','CM','VSMC')
colnames(luci_stats_logfc) <- c('snp_info','sentinel','gene','CM','VSMC')

### plot boxplot of mpra and median -----
head(all)
head(luci_fc_results_m)
head(luci_stats_logfc)

ggplot(luci_stats_logfc, aes(x=gene, y=CM))+
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) +
  theme_minimal() + geom_point(data=all,aes(x=gene, y=logFC), color='red')

