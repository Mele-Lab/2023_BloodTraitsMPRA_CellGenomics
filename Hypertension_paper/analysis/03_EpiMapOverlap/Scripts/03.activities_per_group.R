### Script to plot activities of varints comparing genomic locations
library(tidyr)
library(dplyr)
library(reshape2)

tissues <- c('CARDIAC.MYOCYT_','CARD.MUSCL_','BRN.VASC.SMTH.MUSC_','CORONARY.ATY2_','SM.MUSCLE_')

enh_prom_elements <- readRDS('../../../analysis/03_EpiMapOverlap/enh_prom_elements_overlap_per_tissue.rds')
sentinel_snps_merged <- read.table('../../../analysis/03_EpiMapOverlap/sentinel_snps_closest_element.txt',
            sep = '\t', header = T)
enh_prom_repr <- readRDS('~../../../analysis/03_EpiMapOverlap/proportions_enh_prom_repr.rds')

#### annotate elements #####
for (tissue in tissues) {
  sentinel_snps_merged[[tissue]] <- sentinel_snps_merged$DHS
  sentinel_snps_merged[[tissue]][sentinel_snps_merged$change %in% enh_prom_repr[[tissue]]$Enh] <- 'Enhancer'
  sentinel_snps_merged[[tissue]][sentinel_snps_merged$change %in% enh_prom_repr[[tissue]]$Prom] <- 'Promoter'
  sentinel_snps_merged[[tissue]][sentinel_snps_merged$change %in% enh_prom_repr[[tissue]]$Repressed] <- 'Repressed'
  }

annotation <- melt(sentinel_snps_merged[,c('change',tissues)], id.vars = 'change')

#### Do boxplots with abs(logFC) and MPRAact #####
#### MPRA activity ####

vals_significance_cm <- read.table('../../../data/02_activs/CM_vals.significance.5.txt', sep = '\t', header = T)
head(vals_significance_cm)
vals_significance_cm <- vals_significance_cm %>%
  separate(name, c("pos", "snp_info","info"), "__")
head(vals_significance_cm)

vals_significance_vsmc <- read.table('../../../data/02_activs/VSMC_vals.significance.4.txt', sep = '\t', header = T)
head(vals_significance_vsmc)
vals_significance_vsmc <- vals_significance_vsmc %>%
  separate(name, c("pos", "snp_info","info"), "__")
head(vals_significance_vsmc)

vals_significance_cm_dhs <- merge(vals_significance_cm, annotation, by.x = 'info', by.y = 'change')
vals_significance_vsmc_dhs <- merge(vals_significance_vsmc, annotation, by.x = 'info', by.y = 'change')

vals_significance_cm_dhs$value <- factor(vals_significance_cm_dhs$value, levels = c('No','DHS_EpiMap','Enhancer','Promoter','Repressed'))
vals_significance_vsmc_dhs$value <- factor(vals_significance_vsmc_dhs$value, levels = c('No','DHS_EpiMap','Enhancer','Promoter','Repressed'))


vals_significance_cm$group <- 'Not open'
vals_significance_cm$group[vals_significance_cm$info %in% sentinel_snps_merged$change[sentinel_snps_merged$DHS == 'DHS_EpiMap']] <- 'DHS_EpiMap' 
vals_significance_cm$group[vals_significance_cm$info %in% sentinel_snps_merged$change[sentinel_snps_merged$Enhancer != 'No']] <- 'Enhancer' 
vals_significance_cm$group[vals_significance_cm$info %in% sentinel_snps_merged$change[sentinel_snps_merged$Promoter != 'No']] <- 'Promoter' 
#vals_significance_cm$group[vals_significance_cm$info %in% DHS_merged$V8[DHS_merged$dyadic != 'No']] <- 'Dyadic' 

table(vals_significance_cm$group)
vals_significance_cm <- vals_significance_cm[vals_significance_cm$info %in% sentinel_snps_merged$change,]

vals_significance_vsmc$group <- 'Not open'
vals_significance_vsmc$group[vals_significance_vsmc$info %in% sentinel_snps_merged$change[sentinel_snps_merged$DHS == 'DHS_EpiMap']] <- 'DHS_EpiMap' 
vals_significance_vsmc$group[vals_significance_vsmc$info %in% sentinel_snps_merged$change[sentinel_snps_merged$Enhancer != 'No']] <- 'Enhancer' 
vals_significance_vsmc$group[vals_significance_vsmc$info %in% sentinel_snps_merged$change[sentinel_snps_merged$Promoter != 'No']] <- 'Promoter' 
#vals_significance_cm$group[vals_significance_cm$info %in% DHS_merged$V8[DHS_merged$dyadic != 'No']] <- 'Dyadic' 

table(vals_significance_vsmc$group)
vals_significance_vsmc <- vals_significance_vsmc[vals_significance_vsmc$info %in% sentinel_snps_merged$change,]

vals_significance_cm$group <- factor(vals_significance_cm$group, levels = c('Not open','DHS_EpiMap','Enhancer','Promoter'))
vals_significance_vsmc$group <- factor(vals_significance_vsmc$group, levels = c('Not open','DHS_EpiMap','Enhancer','Promoter'))


foo <- pairwise.wilcox.test(vals_significance_cm$CM, vals_significance_cm$group, p.adjust.method="BH")
foo

foo <- pairwise.wilcox.test(vals_significance_vsmc$VSMC, vals_significance_vsmc$group, p.adjust.method="BH")
foo

# 
library(ggpubr)
 #c("Other", "Other DHS"), c("Other", "Enhancer Core"), c("Other", "Dyadic"),
my_comparisons <- list( 
   c('Not open','Promoter'), c('Promoter','DHS_EpiMap'),c('Promoter','Enhancer'))
# 
ggplot(vals_significance_vsmc, aes(x = group, y = VSMC_log)) +
  geom_violin(aes(fill=group), notch=F, outlier.size=0.5, outlier.alpha = 0.5) +
  ggtitle('') +
  #stat_summary(fun.y=mean, geom="point", shape=23, size=4)+
  xlab('MPRA Activity')+
  ylab('')+
  #ylim(c(-10,10))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'))+
  scale_fill_brewer(palette="Reds") + theme_classic() +
  stat_compare_means(comparisons = my_comparisons, label.y = c(10, 12, 14,16))
# 

ggplot(vals_significance_vsmc, aes(x = VSMC_log, y = group, fill=group)) +
  geom_density_ridges(rel_min_height = 0.005, quantile_lines = TRUE, alpha = 0.75,
                      quantiles = 2, scale = 4) + 
  xlim(c(-3,2.5)) + theme_classic() + 
  theme(legend.position = 'None')+
  scale_fill_brewer(palette="Reds") +
  ylab('')

ggplot(vals_significance_cm, aes(CM_log, color=group, 
                 group=group)) + 
  stat_ecdf(geom = "step")+ xlim(c(-3,3))+
  labs(title="CM MPRA activity",
       y = "Empirical Cumulative density", x="MPRA Activity")+
  theme_classic()


#install.packages("ggridges")
library(ggridges)
for (tissue in tissues) {
  foo <- pairwise.wilcox.test(vals_significance_vsmc_dhs$VSMC[vals_significance_vsmc_dhs$variable == tissue], 
                              vals_significance_vsmc_dhs$value[vals_significance_vsmc_dhs$variable == tissue], p.adjust.method="none")
  print(tissue)
  print(foo)
}


ggplot(vals_significance_cm_dhs, aes(x = CM_log, y = value, fill=value, height = stat(density))) +
  geom_density_ridges(rel_min_height = 0.005, quantile_lines = TRUE, alpha = 0.5,
                      quantiles = 2, stat = "density", scale = 4) + 
  facet_grid(variable ~ .) +
  xlim(c(-3,2.5)) + theme_classic() + 
  theme(legend.position = 'right')+
  scale_fill_brewer(palette="Set1") +
  ylab('')

ggplot(vals_significance_vsmc_dhs, aes(x = VSMC_log, y = value, fill=value, height = stat(density))) +
  geom_density_ridges(rel_min_height = 0.005, quantile_lines = TRUE, alpha = 0.5,
                      quantiles = 2, stat = "density", scale = 4) + 
  facet_grid(variable ~ .) +
  xlim(c(-3,3)) + theme_classic() + 
  theme(legend.position = 'right')+
  scale_fill_brewer(palette="Set1") +
  ylab('')#650x570

ggplot(vals_significance_cm_dhs, aes(CM_log, color=value, 
                                 group=value)) + 
  stat_ecdf(geom = "step")+ xlim(c(-2.5,2))+
  facet_grid(variable ~ .) + 
  labs(title="CM MPRA activity",
       y = "Empirical Cumulative density", x="MPRA Activity")+
  theme_classic()

ggplot(vals_significance_cm_dhs, aes(x = log2(CM), y = value, fill=value)) +
  geom_violin(alpha=0.8,
    # Notch?
    notch=T,
    notchwidth = 0.8,
    # custom outliers
    outlier.shape =NA)+
  facet_grid(variable ~ .)+
  xlim(c(-8,5)) + 
  theme_classic() + 
  theme(legend.position = 'right')+
  scale_fill_brewer(palette="Set1") +
  ylab('')

### index 
index <- read.table('../../../data/design/Hypertension__pooled.index.txt', sep = '\t', header = T)
head(index)

### read differential results ####
CM_results <- read.table('../../../analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_with_significance.CM.new_back.005.txt', sep='\t', header = T)
VSMC_results <- read.table('../../../analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_with_significance.VSMC.new_back.005.txt', sep='\t', header = T)

#### merge results with index ####
index['dupe_info'] <- gsub('\\..*','',index$tile_id)

CM_results = merge(CM_results,index[index$tile_type != 'RANDOM',c("name","dupe_info")], by="dupe_info")
VSMC_results = merge(VSMC_results,index[index$tile_type != 'RANDOM',c("name","dupe_info")], by="dupe_info")

CM_results <- CM_results %>% distinct()
VSMC_results <- VSMC_results %>% distinct()

CM_results$info <- gsub('.*__','',CM_results$name)
VSMC_results$info <- gsub('.*__','',VSMC_results$name)


CM_results$group <- 'Not open'
CM_results$group[CM_results$coord %in% sentinel_snps_merged$change[sentinel_snps_merged$DHS =='DHS_EpiMap']] <- 'DHS_EpiMap' 
CM_results$group[CM_results$coord %in% sentinel_snps_merged$change[sentinel_snps_merged$Enhancer != 'No']] <- 'Enhancer' 
CM_results$group[CM_results$coord %in% sentinel_snps_merged$change[sentinel_snps_merged$Promoter != 'No']] <- 'Promoter' 
#CM_results$group[CM_results$info %in% DHS_merged$V8[DHS_merged$dyadic != 'No']] <- 'Dyadic' 

table(CM_results$group)
CM_results <- CM_results[CM_results$coord %in% sentinel_snps_merged$change,]

VSMC_results$group <- 'Not open'
VSMC_results$group[VSMC_results$coord %in% sentinel_snps_merged$change[sentinel_snps_merged$DHS =='DHS_EpiMap']] <- 'DHS_EpiMap' 
VSMC_results$group[VSMC_results$coord %in% sentinel_snps_merged$change[sentinel_snps_merged$Enhancer != 'No']] <- 'Enhancer' 
VSMC_results$group[VSMC_results$coord %in% sentinel_snps_merged$change[sentinel_snps_merged$Promoter != 'No']] <- 'Promoter' 
#VSMC_results$group[VSMC_results$info %in% DHS_merged$V8[DHS_merged$dyadic != 'No']] <- 'Dyadic' 
VSMC_results <- VSMC_results[VSMC_results$coord %in% sentinel_snps_merged$change,]
table(VSMC_results$group)

CM_results$group <- factor(CM_results$group, levels = c('Not open','DHS_EpiMap','Enhancer','Promoter'))
VSMC_results$group <- factor(VSMC_results$group, levels = c('Not open','DHS_EpiMap','Enhancer','Promoter'))


foo <- pairwise.wilcox.test(abs(CM_results$logFC_comp), CM_results$group, p.adjust.method="none")
foo

library(ggpubr)
#c("Other", "Other DHS"), c("Other", "Enhancer Core"), c("Other", "Dyadic"),
my_comparisons <- list( 
  c('Other','Promoter'), c('Promoter','DHS_EpiMap'),c('Promoter','Enhancer'))

ggplot(VSMC_results, aes(x = group, y = abs(logFC_comp))) +
  geom_violin(aes(fill=group), notch=F, outlier.size=0.5, outlier.alpha = 0.5) +
  ggtitle('') +
  #stat_summary(fun.y=mean, geom="point", shape=23, size=4)+
  ylab('MPRA Activity')+
  xlab('')+
  #ylim(c(0,9))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'))+
  scale_fill_brewer(palette="Reds") + theme_classic() +
  stat_compare_means(comparisons = my_comparisons, label.y = c(4, 5, 6,7))


install.packages("ggridges")
library(ggridges)

ggplot(CM_results, aes(x = abs(logFC_comp), y = group, fill=group)) +
  geom_density_ridges(rel_min_height = 0.005, quantile_lines = TRUE, alpha = 0.75,
                      quantiles = 2, scale = 4) + 
  xlim(c(-0.5,3)) + theme_classic() + 
  theme(legend.position = 'None')+
  scale_fill_brewer(palette="Reds") +
  ylab('')

ggplot(CM_results, aes(x = abs(logFC_comp), y= group, fill = group)) +
  geom_density() +
  scale_fill_brewer(palette = "Reds")

ggplot(CM_results, aes(abs(logFC_comp), color=group, 
                                     group=group)) + 
  stat_ecdf(geom = "step")+ xlim(c(-2.5,2))+
  labs(title="CM MPRA effect size",
       y = "Empirical Cumulative density", x="MPRA Activity")+
  theme_classic()

foo <- pairwise.wilcox.test(CM_results$logFC_comp, CM_results$group, p.adjust.method="none")
foo


CM_results_dhs <- merge(CM_results, annotation, by.x = 'coord', by.y = 'change')
VSMC_results_dhs <- merge(VSMC_results, annotation, by.x = 'coord', by.y = 'change')

CM_results_dhs$value <- factor(CM_results_dhs$value, levels = c('Not open','DHS_EpiMap','Enhancer','Promoter','Repressed'))
VSMC_results_dhs$value <- factor(VSMC_results_dhs$value, levels = c('Not open','DHS_EpiMap','Enhancer','Promoter','Repressed'))

#install.packages("ggridges")
library(ggridges)

CM_results_dhs$new <- CM_results_dhs$value
CM_results_dhs$new[CM_results_dhs$value %in% c('Enhancer','Promoter')] <- 'Regulatory'

for (tissue in tissues) {
  foo <- pairwise.wilcox.test(abs(VSMC_results_dhs$logFC_comp[VSMC_results_dhs$variable == tissue]), 
                              VSMC_results_dhs$value[VSMC_results_dhs$variable == tissue], p.adjust.method="none")
  print(tissue)
  print(foo)
}

ggplot(CM_results_dhs, aes(x = abs(logFC_comp), y = value, fill=value)) +
  geom_density_ridges(rel_min_height = 0.005, quantile_lines = TRUE, alpha = 0.5,
                      quantiles = 2, scale = 4) + 
  facet_grid(variable ~ .) +
  xlim(c(0,2.5)) + theme_classic() + 
  theme(legend.position = 'right')+
  scale_fill_brewer(palette="Set1") +
  ylab('')

ggplot(VSMC_results_dhs[VSMC_results_dhs$fdr_comp < 0.05,], aes(x = abs(logFC_comp), y = value, fill=value)) +
  geom_density_ridges(rel_min_height = 0.005, quantile_lines = TRUE, alpha = 0.5,
                      quantiles = 2, scale = 4) + 
  facet_grid(variable ~ .) +
  xlim(c(0,2.5)) + theme_classic() + 
  theme(legend.position = 'right')+
  scale_fill_brewer(palette="Set1") +
  ylab('')

ggplot(CM_results_dhs, aes(x = abs(logFC_comp), y = variable, fill=value)) +
  geom_boxplot(alpha=0.8,
               # Notch?
               notch=FALSE,
               notchwidth = 0.8,
               # custom outliers
               outlier.shape =NA)+
  xlim(c(0,2.5)) + theme_classic() + 
  theme(legend.position = 'right')+
  scale_fill_brewer(palette="Set1") +
  ylab('')

ggplot(CM_results_dhs[!is.na(CM_results_dhs$value),], aes(abs(logFC_comp), color=value, 
                                     group=value)) + 
  stat_ecdf(geom = "step")+ xlim(c(-0.5,3))+
  facet_grid(variable ~ .) + 
  labs(title="CM MPRA effect size",
       y = "Empirical Cumulative density", x="MPRA Activity")+
  theme_classic()

