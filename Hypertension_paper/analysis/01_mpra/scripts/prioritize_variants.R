##### create score for putative causal variants -------
library(dplyr)
library(tidyr)

#### read data ------

##### info ----
CM_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/CM_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt', sep='\t', header=T)
VSMC_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/VSMC_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt', sep='\t', header=T)

head(CM_results)

##### DHS -----
DHS <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/03_EpiMapOverlap/ALL_DHS_Epimap_core_SNPS_overlap.bed', header = F, sep = '\t')

### how many are enhancers? ####
#### these data come from the EpiMap repository #####
enhancers <- read.table('~/marenostrum/Data/Hypertension/EpiMap/ENH_masterlist_indices_0indexed.tsv', header = F)
map_row_id <- read.table('~/marenostrum/Data/Hypertension/EpiMap/masterlist_DHSs_733samples_WM20180608_all_coords_hg19_r25_e100_names.core.srt.tsv', header = T)
old_new_id <- read.table('~/marenostrum/Data/Hypertension/EpiMap/masterlist_DHSs_733samples_WM20180608_all_chunkIDs2indexIDs.txt', header = F)
head(enhancers)
head(map_row_id)
head(DHS)
head(old_new_id)

enhancers_old_id <- map_row_id$name[map_row_id$cls %in% enhancers$V1]
new_id_enhancers <- old_new_id$V1[old_new_id$V2 %in% enhancers_old_id]

DHS$enhancers <- 'No'
DHS$enhancers[DHS$V4 %in% new_id_enhancers] <- 'Enhancer'

### how many are promoters? ####
#### these data come from the EpiMap repository #####
promoters <- read.table('~/marenostrum/Data/Hypertension/EpiMap/PROM_masterlist_indices_0indexed.tsv', header = F)

promoters_old_id <- map_row_id$name[map_row_id$cls %in% promoters$V1]
new_id_promoters <- old_new_id$V1[old_new_id$V2 %in% promoters_old_id]

DHS$promoters <- 'No'
DHS$promoters[DHS$V4 %in% new_id_promoters] <- 'Promoter'

length(unique(DHS$V8[DHS$promoters != 'No']))
#[1] 225
length(unique(DHS$V4[DHS$promoters != 'No']))
#[1] 509

### how many are dyadic? ####
#### these data come from the EpiMap repository #####
dyadic <- read.table('~/marenostrum/Data/Hypertension/EpiMap/DYADIC_masterlist_indices_0indexed.tsv', header = F)

dyadic_old_id <- map_row_id$name[map_row_id$cls %in% dyadic$V1]
new_id_dyadic <- old_new_id$V1[old_new_id$V2 %in% dyadic_old_id]

DHS$dyadic <- 'No'
DHS$dyadic[DHS$V4 %in% new_id_dyadic] <- 'Dyadic'

#### merge with sentinel info ####
CM_results$change <- paste0(CM_results$chr, ':',CM_results$start, ':',CM_results$ref,':',CM_results$alt)

CM_results$DHS <- 'No'
CM_results$DHS[CM_results$change %in% DHS$V8] <- 'DHS_EpiMap'

CM_results$Enhancer <- 'No'
CM_results$Enhancer[CM_results$change %in% DHS$V8[DHS$enhancers!='No']] <- 'Enhancer'

CM_results$Promoter <- 'No'
CM_results$Promoter[CM_results$change %in% DHS$V8[DHS$promoters!='No']] <- 'Promoter'

CM_results$Dyadic <- 'No'
CM_results$Dyadic[CM_results$change %in% DHS$V8[DHS$dyadic!='No']] <- 'Dyadic'

#VSMC
VSMC_results$change <- paste0(VSMC_results$chr, ':',VSMC_results$start, ':',VSMC_results$ref,':',VSMC_results$alt)

VSMC_results$DHS <- 'No'
VSMC_results$DHS[VSMC_results$change %in% DHS$V8] <- 'DHS_EpiMap'

VSMC_results$Enhancer <- 'No'
VSMC_results$Enhancer[VSMC_results$change %in% DHS$V8[DHS$enhancers!='No']] <- 'Enhancer'

VSMC_results$Promoter <- 'No'
VSMC_results$Promoter[VSMC_results$change %in% DHS$V8[DHS$promoters!='No']] <- 'Promoter'

VSMC_results$Dyadic <- 'No'
VSMC_results$Dyadic[VSMC_results$change %in% DHS$V8[DHS$dyadic!='No']] <- 'Dyadic'

### ultraconserved element ----
ultra_elements <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/overlap_ultraconserved.bed',
                             header=F, sep = '\t')
head(ultra_elements)

CM_results$ultraconserved <- 0
CM_results$ultraconserved[CM_results$change %in% ultra_elements$V4] <- 1
table(CM_results$ultraconserved)

VSMC_results$ultraconserved <- 0
VSMC_results$ultraconserved[VSMC_results$change %in% ultra_elements$V4] <- 1
table(VSMC_results$ultraconserved)

#### write preliminar table ----
write.table(CM_results,'~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/CM_eqtls_rep_cons_dhs_tbs.txt', sep = '\t',
            quote = F,col.names = T, row.names = F)

write.table(VSMC_results,'~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/VSMC_eqtls_rep_cons_dhs_tbs.txt', sep = '\t',
            quote = F,col.names = T, row.names = F)

#### read results -----
CM_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/CM_eqtls_rep_cons_dhs_tbs.txt', sep = '\t',
            header = T)

VSMC_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/VSMC_eqtls_rep_cons_dhs_tbs.txt', sep = '\t',
            header = T)
head(CM_results)

VSMC_results[,c('X','index','stat_comp','df.test_comp','df.dna_comp','df.rna.full_comp','df.rna.red_comp',
                'is_ctrl','tile_type','native_status','V7','V8','end')] <- NULL
CM_results[,c('X','index','stat_comp','df.test_comp','df.dna_comp','df.rna.full_comp','df.rna.red_comp',
              'is_ctrl','tile_type','native_status','V7','V8','end')] <- NULL

VSMC_results <- VSMC_results[,colnames(CM_results)]
eqtls_cols <- grep('eQTL',colnames(VSMC_results))

##### score number 2 ------

#### merge with sentinel info---
VSMC_results_s <- merge(VSMC_results, sentinel_snps[,c('sentinel','snp')], by.x = 'snp_info', by.y = 'snp')
CM_results_s <- merge(CM_results, sentinel_snps[,c('sentinel','snp')], by.x = 'snp_info', by.y = 'snp')

top_VSMC <- VSMC_results_s %>% group_by(sentinel) %>% top_n(-1, rankingFDR)
top_CM <- CM_results_s %>% group_by(sentinel) %>% top_n(-1, rankingFDR)

get_groups_2 <- function(x, output) {
  score <- 1
  if (x['ultraconserved'] == 1) {
    score <- score + 1
  }
  if ('YES' %in% (x[eqtls_cols])) {
    score <- score + 1
  }
  if (x['sentinel'] == x['snp_info']) {
    score <- score + 1
  }
  if (x['snp'] %in% top_CM$snp_info) {
    score <- score + 1
  }
  if (x['DHS'] == 'DHS_EpiMap') {
    score <- score + 1
  }
  if (!'' %in% (x[c('TF_WT','TF_SNP')])) {
    score <- score + 1
  }
  if (x['active'] == 'Active') {
    score <- score + 1
  }
  if (x['overlap'] == 'Differential Overlap') {
    score <- score + 1
  }
  if (x['overlap'] == 'Differential + Active Overlap') {
    score <- score + 1
  }
  return(score)

}

VSMC_results_s$score_groups_2 <-apply(VSMC_results_s, 1, get_groups_2)
CM_results_s$score_groups_2 <-apply(CM_results_s, 1, get_groups_2)

#### write scores table ----
write.table(CM_results_s,'~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/CM_eqtls_rep_cons_dhs_tbs_scores.txt', sep = '\t',
            quote = F,col.names = T, row.names = F)

write.table(VSMC_results_s,'~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/VSMC_eqtls_rep_cons_dhs_tbs_scores.txt', sep = '\t',
            quote = F,col.names = T, row.names = F)

## some plots ---
ggplot(VSMC_results_s, aes(x = score_groups_2)) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-1)+
  ylim(c(0,200))+ xlab('score 1-10') + ylab('Number of variants')+
  theme_minimal()


#### number of eQTLs per sentinel ----

get_n_eqtls <- function(x, output) {
  score <- 0
  if ('YES' %in% (x[eqtls_cols])) {
    score <- score + 1
  }
   return(score)
}

CM_results_s$isqtl <-apply(CM_results_s, 1, get_n_eqtls)
VSMC_results_s$isqtl <-apply(VSMC_results_s, 1, get_n_eqtls)

library(dplyr)
eqtls_sentinel_CM <- CM_results_s %>%
  group_by(sentinel) %>%
  dplyr::summarise(eqtls = sum(isqtl), number = n())

eqtls_sentinel_VSMC <- VSMC_results_s %>%
  group_by(sentinel) %>%
  dplyr::summarise(eqtls = sum(isqtl), number = n())

eqtls_sentinel_CM$freq <- eqtls_sentinel_CM$eqtls/eqtls_sentinel_CM$number
eqtls_sentinel_VSMC$freq <- eqtls_sentinel_VSMC$eqtls/eqtls_sentinel_VSMC$number

ggplot(eqtls_sentinel_CM, aes(x=number, y=eqtls)) +
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_minimal() +
  geom_text_repel(
    data=eqtls_sentinel_CM %>% filter(freq>=0.8), # Filter data first
    aes(label=sentinel)
  )


#### plot log2FC of prioritized variants ----
ggplot(CM_results_s, aes(x = abs(logFC_comp), y = as.factor(score_groups_2), fill=as.factor(score_groups_2))) +
  geom_density_ridges(rel_min_height = 0.005, quantile_lines = TRUE, alpha = 0.75,
                      quantiles = 2, scale = 4) +
  #xlim(c(-3,2.5)) +
  theme_classic() +
  theme(legend.position = 'None')+
  scale_fill_brewer(palette="Reds") +
  ylab('')

my_comparisons <- list(
  c('1','5'), c('1','4'),c('1','3'),c('1','2'))
my_comparisons <- list(
   c('1','4'),c('1','3'),c('1','2'))

##foldchange
ggplot(VSMC_results_s[VSMC_results_s$score_groups_2>0,], aes(x = as.factor(score_groups_2), y = abs(logFC_comp))) +
  geom_violin(aes(fill=as.factor(score_groups_2)), width=1) +
  geom_boxplot(width=0.1, color="black", alpha=0.2)+
  ggtitle('') +
  xlab('Group 0-5')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'))+
  scale_fill_brewer(palette="Reds") + theme_classic() +
  stat_compare_means(comparisons = my_comparisons, label.y = c(4, 4.3, 4.5,4.7))

CM_results_s$group <- 'Non prioritized'
CM_results_s$group[CM_results_s$score_groups_2 > 3] <- 'Prioritized'

ggplot(CM_results_s, aes(x = as.factor(group), y = abs(logFC_comp))) +
  geom_violin(aes(fill=as.factor(group)), width=1) +
  geom_boxplot(width=0.1, color="black", alpha=0.2)+
  ggtitle('') +
  xlab('Group 0-5')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'))+
  scale_fill_brewer(palette="Reds") + theme_classic() +
  stat_compare_means()


##pval
ggplot(CM_results_s[CM_results_s$score_groups_2>0,], aes(x = as.factor(score_groups_2), y = -log10(fdr_comp))) +
  geom_violin(aes(fill=as.factor(score_groups_2)), width=1) +
  geom_boxplot(width=0.1, color="black", alpha=0.2)+
  ggtitle('') +
  xlab('Group 0-5')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'))+
  scale_fill_brewer(palette="Reds") + theme_classic() +
  stat_compare_means(comparisons = my_comparisons, label.y = c(100, 110, 120,130))

ggplot(CM_results_s, aes(x = as.factor(group), y = -log10(fdr_comp))) +
  geom_violin(aes(fill=as.factor(group)), width=1) +
  geom_boxplot(width=0.1, color="black", alpha=0.2)+
  ggtitle('') +
  xlab('Group 0-5')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'))+
  scale_fill_brewer(palette="Reds") + theme_classic() +
  stat_compare_means()

#### number of differential TFBS -----
CM_results_s$group <- 'Non prioritized'
CM_results_s$group[CM_results_s$score_groups_2 > 3] <- 'Prioritized'

ggplot(CM_results_s, aes(x = as.factor(group), y = abs(dist))) +
  geom_violin(aes(fill=as.factor(group)), width=1) +
  geom_boxplot(width=0.1, color="black", alpha=0.2)+
  ggtitle('') +
  xlab('Group 0-5')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'))+
  scale_fill_brewer(palette="Reds") + theme_classic() +
  stat_compare_means()


### how many SNPs go in the same direction? per sentinel ------
head(CM_results_s)

fold_change <- CM_results_s %>%
  group_by(sentinel) %>%
  dplyr::summarise(positive = sum(logFC_comp > 0), negative = sum(logFC_comp < 0), total = n())

ggplot(fold_change[fold_change$total>20,], aes(sentinel), ylim(-180:180)) + 
  geom_bar(data = fold_change[fold_change$total>20,], 
           aes(y = positive), stat = "identity", position = "dodge", fill='dark blue') +
  geom_bar(data = fold_change[fold_change$total>20,], 
           aes(y = -negative), stat = "identity", position = "dodge", fill='dark red') + 
  geom_hline(yintercept = 0,colour = "grey90") + ylab('Nº reg variants') +
  xlab('')+
  coord_flip()+ theme_minimal()

last_plot() + 
  geom_text(data = fold_change[fold_change$total>20,], 
            aes(sentinel, positive, label=positive),
            position = position_dodge(width=0.9), vjust = -0.1, size=4) +
  geom_text(data = fold_change[fold_change$total>20,], 
            aes(sentinel, -negative, label=negative),
            position = position_dodge(width=0.9), vjust = -0, size=4) +
  coord_cartesian(ylim = c(-180, 180)) +
  coord_flip()


#### upset plot of number of functional annotations ######
CM_results_s <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/CM_eqtls_rep_cons_dhs_tbs_scores.txt', sep = '\t',
                           header = T)

VSMC_results_s <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/VSMC_eqtls_rep_cons_dhs_tbs_scores.txt', sep = '\t',
                             header = T)
ID <- 1:7

ggplot(data=CM_results_s, aes(x=score_groups_2)) +
  geom_bar( fill=wes_palette("GrandBudapest1")[4]) +
  geom_text(stat='count', aes(label=..count..), vjust=-0) + 
  theme_light() + 
  ylab('Nº of regulatory variants') +
  scale_x_continuous("Score", labels = as.character(ID), breaks = ID)

eqtls_cols <- grep('eQTL',colnames(CM_results_s))

library(ComplexHeatmap)
lt <- list(eQTL = CM_results_s$snp_info[rowSums(sapply(CM_results_s[,eqtls_cols],`==`,e2='YES'), na.rm = T)>=1],
           DHS = CM_results_s$snp_info[CM_results_s$DHS == 'DHS_EpiMap'],
           Enhancer= CM_results_s$snp_info[CM_results_s$Enhancer == 'Enhancer'],
           Promoter=CM_results_s$snp_info[CM_results_s$Promoter == 'Promoter'],
           Ultraconserved=CM_results_s$snp_info[CM_results_s$ultraconserved == 1],
           TFBS = CM_results_s$snp_info[CM_results_s$N_diff_motifs != 0])

m = make_comb_mat(lt)
UpSet(m, comb_order = order(comb_size(m), decreasing = T))

ht = draw(UpSet(m, 
                comb_order = order(comb_size(m), decreasing = T),
                right_annotation = upset_right_annotation(m,
                                                          add_numbers = TRUE),
                top_annotation = upset_top_annotation(m, add_numbers = TRUE)))
print(ht)
od = column_order(ht)
cs = comb_size(m)
decorate_annotation("intersection_size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = "bottom", gp = gpar(fontsize = 8))
})


