library(dplyr)
library(tidyr)
library(ggplot2)

### read fimo results ####
tile_motifs = read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/tiles_tfbs.best_motifs.csv', sep='\t', header = T)

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
  dif_tf <- setdiff(wildtype_tfbs$motif_alt_id, snp_tfbs$motif_alt_id)
  if (n_same_tf >= 1) {
    tfbs_merged <- merge(as.data.frame(table(wildtype_tfbs$motif_alt_id)),as.data.frame(table(snp_tfbs$motif_alt_id)), by='Var1')
    tfbs_merged$diff <-  tfbs_merged$Freq.x - tfbs_merged$Freq.y
    n_same_tfb <- sum(tfbs_merged$Freq.y[tfbs_merged$diff == 0])
  }else{
    n_same_tfb <- 0
  }
  snps <- append(snps, snp)
  sameTF <- append(sameTF, n_same_tf)
  sameTFBS <- append(sameTFBS, length(dif_tf))
}
snps_list <- unlist(snps)
names(sameTF) <- snps_list
names(sameTFBS) <- snps_list
results <- list(as.data.frame(unlist(snps)), as.data.frame(unlist(sameTF)), as.data.frame(unlist(sameTFBS)))
results_df <- do.call('cbind.data.frame', results)
colnames(results_df) <- c('snp_info','n_same_tfbs','diff_TF')

#### read diff TF ####
fimo_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/03_fimo/diff_TFBS.csv', sep='\t', header = T)

#### read sentinel SNPs ####
sentinel_snps = read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/sentinel_trait_ld_snps.txt', sep = '\t', header = T)
head(sentinel_snps)
sentinel_snps = sentinel_snps[sentinel_snps$Trait %in% c('DBP','PP','SBP'),]

### read differential results ####
CM_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/CM_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt', sep='\t', header=T)
VSMC_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/VSMC_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt', sep='\t', header=T)

### plot TF vs differential ####
fimo_results <- fimo_results[fimo_results$SNP %in% sentinel_snps$snp,]
results_df <- results_df[results_df$snp_info %in% sentinel_snps$snp,]

results_df$disrupted[results_df$diff_TF == 0] <- 'Same TF'
results_df$disrupted[results_df$diff_TF != 0] <- 'Disrupted TF'

diff_snps <- unique(CM_results$snp_info[CM_results$fdr_comp<=0.05])
diff_snps_vsmc <- unique(VSMC_results$snp_info[VSMC_results$fdr_comp<=0.05])

results_df$DiffCM <- 'Not Diff'
results_df$DiffCM[results_df$snp_info %in% diff_snps] <- 'Diff'  

results_df$DiffVSMC <- 'Not Diff'
results_df$DiffVSMC[results_df$snp_info %in% diff_snps_vsmc] <- 'Diff'  

results_df$ActiveCM <- 'Not Active'
results_df$ActiveVSMC <- 'Not Active'

active_snps <- unique(vals_significance_cm$snp_info[vals_significance_cm$CM_padj<0.05])
active_snps_vsmc <- unique(vals_significance_vsmc$snp_info[vals_significance_vsmc$VSMC_padj<0.05])

results_df$ActiveCM[results_df$snp_info %in% active_snps] <- 'Active'
results_df$ActiveVSMC[results_df$snp_info %in% active_snps_vsmc] <- 'Active'

results_df$DiffAct <- 'Not DiffAct'
results_df$DiffAct[results_df$DiffCM == 'Diff' & results_df$ActiveCM == 'Active'] <- 'DiffAct'  

results_df$DiffActVSMC <- 'Not DiffAct'
results_df$DiffActVSMC[results_df$DiffVSMC == 'Diff' & results_df$ActiveVSMC == 'Active'] <- 'DiffAct'  


cm_tf <- as.data.frame(table(results_df$disrupted, results_df$DiffActVSMC))
cm_tf$Perc <- ""
for (element in c('Disrupted TF','Same TF')){
  for (elem2 in c('DiffAct','Not DiffAct')) {
    cm_tf$Perc[cm_tf$Var1 == element & cm_tf$Var2 == elem2] <- cm_tf$Freq[cm_tf$Var1 == element & cm_tf$Var2 == elem2]/sum(cm_tf$Freq[cm_tf$Var2 == elem2])
    
  }
} 

ggplot(cm_tf[cm_tf$Var1 == 'Disrupted TF',], aes(x=Var2, y= as.numeric(Perc),fill=Var2)) + 
  geom_bar(stat="identity", position="dodge") + 
  scale_fill_brewer(palette = "Set2") + 
  xlab('')+
  ylim(c(0,0.8))+
  ylab('Prop of SNPs with disrupted TF')+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'), 
        axis.text = element_text(size = 10))

#### have variants in promoters higher numbers of TFs¿
### read differential results ####
CM_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/all_snps_info_risk_CM.txt', sep='\t', header=T)
VSMC_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/all_snps_info_risk_VSMC.txt', sep='\t', header=T)

tfbs <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/02_TF_analysis/turnover_results_tf_diff.txt')
head(tfbs)

index <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/Hypertension__pooled.index.txt', sep = '\t', header = T)
head(index)

tfbs <- merge(tfbs, index[,c("name", "element")] %>% distinct(), by.x = 'element_ref', by.y = 'element')
head(tfbs)

#### promoter annotations 
library(tidyr)
library(dplyr)
library(reshape2)
tissues <- c('CARDIAC.MYOCYT_','CARD.MUSCL_','BRN.VASC.SMTH.MUSC_','CORONARY.ATY2_','SM.MUSCLE_')

sentinel_snps_merged <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/03_EpiMapOverlap/sentinel_snps_closest_element.txt',
                                   sep = '\t', header = T)
enh_prom_repr <- readRDS('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/03_EpiMapOverlap/proportions_enh_prom_repr.rds')

#### annotate elements #####
for (tissue in tissues) {
  sentinel_snps_merged[[tissue]] <- sentinel_snps_merged$DHS
  sentinel_snps_merged[[tissue]][sentinel_snps_merged$change %in% enh_prom_repr[[tissue]]$Enh] <- 'Enhancer'
  sentinel_snps_merged[[tissue]][sentinel_snps_merged$change %in% enh_prom_repr[[tissue]]$Prom] <- 'Promoter'
  sentinel_snps_merged[[tissue]][sentinel_snps_merged$change %in% enh_prom_repr[[tissue]]$Repressed] <- 'Repressed'
}

CM_results$group <- 'Not open'
CM_results$group[CM_results$change %in% sentinel_snps_merged$change[sentinel_snps_merged$DHS == 'DHS_EpiMap']] <- 'DHS_EpiMap' 
CM_results$group[CM_results$change %in% sentinel_snps_merged$change[sentinel_snps_merged$Enhancer != 'No']] <- 'Enhancer' 
CM_results$group[CM_results$change %in% sentinel_snps_merged$change[sentinel_snps_merged$Promoter != 'No']] <- 'Promoter' 
#vals_significance_cm$group[vals_significance_cm$info %in% DHS_merged$V8[DHS_merged$dyadic != 'No']] <- 'Dyadic' 

table(CM_results$group)

VSMC_results$group <- 'Not open'
VSMC_results$group[CM_results$change %in% sentinel_snps_merged$change[sentinel_snps_merged$DHS == 'DHS_EpiMap']] <- 'DHS_EpiMap' 
VSMC_results$group[CM_results$change %in% sentinel_snps_merged$change[sentinel_snps_merged$Enhancer != 'No']] <- 'Enhancer' 
VSMC_results$group[CM_results$change %in% sentinel_snps_merged$change[sentinel_snps_merged$Promoter != 'No']] <- 'Promoter' 
#vals_significance_cm$group[vals_significance_cm$info %in% DHS_merged$V8[DHS_merged$dyadic != 'No']] <- 'Dyadic' 

table(VSMC_results$group)

tfbs_cm <- merge(tfbs, CM_results, by = 'name')
head(tfbs_cm)
tfbs_vsmc <- merge(tfbs, VSMC_results, by = 'name')
head(tfbs_vsmc)

ggplot(tfbs_cm, aes(x=group, y=total_motifs, fill=Activity)) + 
  geom_boxplot()+
  theme_bw()

p <- ggboxplot(tfbs_vsmc, x = "group", y = "total_motifs",
               color = "Activity", palette = "jco") + 
  stat_compare_means(aes(group = Activity), label = "p.format")
print(p)

my_comparisons <- list( c("Promoter", "Enhancer"), c("Promoter", "Not open"), c("DHS_EpiMap", "Promoter") )
p <- ggboxplot(tfbs_cm, x = "group", y = "total_motifs",
                palette = "jco") + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format")
print(p)

tfbs_cm$Active <- 'not active'
tfbs_cm$Active[tfbs_cm$Activity=='Active'] <- 'Active in CMs'
tfbs_cm$Active[tfbs_cm$name %in% tfbs_vsmc$name[tfbs_vsmc$Activity == 'Active']] <- 'Active in VSMCs'
tfbs_cm$Active[tfbs_cm$Active == 'Active in VSMCs' & tfbs_cm$Activity=='Active'] <- 'Active in both'
table(tfbs_cm$Active)

p <- ggboxplot(tfbs_cm, x = "group", y = "total_motifs",
               color = "Active", palette = "jco") + 
  stat_compare_means(aes(group = Active), label = "p.format")
print(p)


