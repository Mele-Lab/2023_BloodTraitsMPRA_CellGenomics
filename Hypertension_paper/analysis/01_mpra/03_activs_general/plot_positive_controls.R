library(dplyr)
library(plyr)
library(tidyr)
library(ggplot2)

##### Plot positive controls #######

activities <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/02_activs/alpha_per_elem.quantification.VSMC.4.txt', sep = '\t')
head(activities)

vals_significance <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/02_activs/CM_vals.significance.5.txt', sep = '\t', header = T)
head(vals_significance)
vals_significance <- vals_significance %>%
  separate(name, c("pos", "snp_info","info"), "__")

controls <- c('rs6445040','rs9270898','rs60002611','rs117104239','rs4360494','rs62012628')
controls_data <- vals_significance[vals_significance$snp_info %in% controls,]
head(controls_data)

counts_dir = "~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/01_counts/"
mpranalyze_dir = paste0(counts_dir,"/mpranalyze_files")

### index ####
index = read.table("~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/Hypertension__pooled.index.txt", sep="\t", header = T)
head(index)

### VSMC data ####
VSMC_data = read.table(paste0(counts_dir,"/VSMC__all_counts_final.4.txt"), sep="\t", header = T)
head(VSMC_data)

### CM data ####
CM_data = read.table(paste0(counts_dir,"/CM__all_counts_final.5.txt"), sep="\t", header = T)
head(CM_data)

colnames(VSMC_data) = c("barcode", "dna_1", "VSMC_rep1", "VSMC_rep2", "VSMC_rep3","VSMC_rep6")
colnames(CM_data) = c("barcode", "dna_1", "CM_rep1", "CM_rep2", "CM_rep3", "CM_rep4","CM_rep5")

df_VSMC = merge(VSMC_data, index[,c("barcode", "element", "tile_type", "tile_id","snp","name")], by="barcode")
df_CM = merge(CM_data, index[,c("barcode", "element", "tile_type", "tile_id","snp","name")], by="barcode")
head(df_VSMC)


df_CM$median_rep = rowMeans(df_CM[,c("CM_rep1", "CM_rep2", "CM_rep3", "CM_rep4","CM_rep5")], na.rm = T)
df_VSMC$median_rep = rowMeans(df_VSMC[,c("VSMC_rep1", "VSMC_rep2", "VSMC_rep3","VSMC_rep6")], na.rm = T)

df_CM$dupe_info = gsub('\\..*','',df_CM$tile_id)
df_VSMC$dupe_info = gsub('\\..*','',df_VSMC$tile_id)

df_CM$tile_id <- df_CM$dupe_info
positive_controls_CM <- merge(df_CM, controls_data, by = c("snp","tile_id","tile_type","dupe_info"), all.x = FALSE)

dupe_snp <- positive_controls_VSMC[,c('dupe_info','snp')]
dupe_snp <- dupe_snp[dupe_snp$snp!='none',]
dupe_snp <- dupe_snp %>% distinct()

df_VSMC$tile_id <- df_VSMC$dupe_info
positive_controls_VSMC <- merge(df_VSMC, controls_data, by = c("snp","tile_id","tile_type","dupe_info"))

positive_controls_CM$activity = positive_controls_CM$median_rep/positive_controls_CM$dna_1
positive_controls_VSMC$activity = positive_controls_VSMC$median_rep/positive_controls_VSMC$dna_1

#log10(median_rep)
ggplot(positive_controls_VSMC, aes(x = dupe_info, y = activity, fill=as.factor(is_snp))) +
  #geom_jitter(aes(color=tile_type),alpha=0.4)+
  geom_point(aes(color=as.factor(is_snp)),pch = 21, alpha=0.3, position = position_jitterdodge())+
  geom_boxplot(aes(alpha=factor(sig)),outlier.shape = NA)+
  ylab("Ratio RNA/DNA") + 
  xlab("") + 
  theme(legend.position="top", legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white", color = NA)) +
  labs(color="", fill = "")+
  scale_color_discrete(name = "", labels = c("REF", "ALT"))+
  scale_fill_discrete(name = "", labels = c("REF", "ALT"))+
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(axis.text.x = element_text(angle=45, hjust = 1))+
  #scale_color_brewer(palette="Set1") + 
  scale_x_discrete(breaks = dupe_snp$dupe_info, labels = dupe_snp$snp)+
  #stat_compare_means(aes(group = tile_type), label = "p.signif", hide.ns = T)+
  ylim(c(0,4))

ggplot(positive_controls_CM, aes(x = tile_type, y = activity)) +
  #geom_jitter(aes(color=tile_type),alpha=0.4)+
  geom_point(aes(color=tile_type),pch = 21, alpha=0.3, position = position_jitterdodge())+
  geom_violin(outlier.shape = NA, alpha=0.5)+
  ylab("Ratio RNA/DNA") + 
  xlab("") + 
  scale_color_discrete(name = "", labels = c("REF", "ALT"))+
  scale_fill_discrete(name = "", labels = c("REF", "ALT"))+
  theme(legend.position="top", legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white", color = NA)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(axis.text.x = element_text(angle=45, hjust = 1))+
  ylim(0,10)
#scale_color_brewer(palette="Set1") + 

#### BARPLOT CONTROLS VSMC AND CM ####
all_controls <- merge(controls_data, controls_data_cm)

ggplot(controls_data_vsmc, aes(x = snp_info, y = VSMC, fill=as.factor(is_snp))) +
  geom_bar(aes(alpha=factor(sig), fill=as.factor(is_snp), color= as.factor(is_snp)),stat = 'identity',
           position="dodge")+
  ylab("MPRA Activity") + 
  xlab("") + 
  theme(legend.position="top", legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white", color = NA)) +
  labs(color="", fill = "")+
  scale_fill_discrete(name = "", labels = c("REF", "ALT"))+
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(axis.text.x = element_text(angle=45, hjust = 1))

ggplot(controls_data, aes(x = snp_info, y = CM, fill=as.factor(is_snp))) +
  geom_bar(aes(alpha=factor(sig), fill=as.factor(is_snp), color= as.factor(is_snp)),stat = 'identity',
           position="dodge")+
  ylab("MPRA Activity") + 
  xlab("") + 
  theme(legend.position="top", legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white", color = NA)) +
  labs(color="", fill = "")+
  scale_fill_discrete(name = "", labels = c("REF", "ALT"))+
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(axis.text.x = element_text(angle=45, hjust = 1))
  #ylim(c(0,0.8))

random_vsmc <- vals_significance[vals_significance$tile_type == 'RANDOM',]
ggplot(random_vsmc, aes(x = tile_type, y = VSMC)) +
  geom_violin(fill='dark grey')+
  ylab("MPRA Activity") + 
  xlab("") + 
  theme(legend.position="top", legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white", color = NA)) +
  labs(color="", fill = "")+
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(axis.text.x = element_text(angle=45, hjust = 1))+
  ylim(c(0,0.6))

