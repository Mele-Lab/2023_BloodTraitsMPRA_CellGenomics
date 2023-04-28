#### Differential plots
library(dplyr)
library(plyr)
library(tidyr)
library(ggplot2)

### read counts ####
## CM 
counts_dir = "~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/01_counts"
mpranalyze_dir = paste0(counts_dir,"/mpranalyze_files")

CM_counts = read.table(paste0(counts_dir,"/CM__all_counts_final.5.txt"), sep="\t", header = T)
head(CM_counts)

### VSMC
VSMC_counts = read.table(paste0(counts_dir,"/VSMC__all_counts_final.4.txt"), sep="\t", header = T)
head(VSMC_counts)

### index 
index <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/Hypertension__pooled.index.txt', sep = '\t', header = T)
head(index)

colnames(VSMC_counts) = c("barcode", "dna_1", "VSMC_rep1", "VSMC_rep2", "VSMC_rep3","VSMC_rep6")
colnames(CM_counts)  = c("barcode", "dna_1", "CM_rep1", "CM_rep2", "CM_rep3", "CM_rep4","CM_rep5")


### read differential results ####
CM_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_with_significance.CM.new_back.005.txt', sep='\t', header = T)
VSMC_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_with_significance.VSMC.new_back.005.txt', sep='\t', header = T)

### merge counts with index ####
df_VSMC = merge(VSMC_counts,index[,c("barcode", "element", "tile_type", "tile_id","snp","name")], by="barcode")
df_CM = merge(CM_counts,index[,c("barcode", "element", "tile_type", "tile_id","snp","name")], by="barcode")

#df_VSMC['median_rep'] = rowMeans(df[,c("VSMC_rep1", "VSMC_rep2", "VSMC_rep3","VSMC_rep6")], na.rm = T)
#df_CM['median_rep'] = rowMeans(df[,c("CM_rep1", "CM_rep2", "CM_rep3", "CM_rep4")], na.rm = T)

### get barcodes per snp ####
barc_VSMC = df_VSMC %>% group_by(name) %>% dplyr::summarize(num=n())
barc_CM = df_CM %>% group_by(name) %>% dplyr::summarize(num=n())

write.table(barc_CM, '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/02_activs/barc_recovery_CM.txt',
            sep = '\t', row.names = F, col.names = T, quote = F)
write.table(barc_VSMC, '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/02_activs/barc_recovery_VSMC.txt',
            sep = '\t', row.names = F, col.names = T, quote = F)

#### merge results with index ####
index['dupe_info'] <- gsub('\\..*','',index$tile_id)

CM_results = merge(CM_results,index[index$tile_type != 'RANDOM',c("name","dupe_info")], by="dupe_info")
VSMC_results = merge(VSMC_results,index[index$tile_type != 'RANDOM',c("name","dupe_info")], by="dupe_info")

#### merge results with barcodes ####
CM_results = merge(CM_results,barc_CM)
VSMC_results = merge(VSMC_results,barc_VSMC)

CM_results <- CM_results %>% distinct()
VSMC_results <- VSMC_results %>% distinct()

CM_results
VSMC_results

#### plot volcano plot with barcode data #####

ggplot(VSMC_results[VSMC_results$tile_type == 'WILDTYPE_SNP_INDIV' & VSMC_results$fdr_comp < 0.05,], aes(x=logFC_comp,y=-log10(fdr_comp), color=num))+
  geom_point(data=VSMC_results[VSMC_results$tile_type == 'WILDTYPE_SNP_INDIV' & VSMC_results$fdr_comp > 0.05,], size = 0.5, color = 'grey')+
  geom_point(size = 1)+
  scale_color_gradient(low="blue", high="red")+
  xlim(c(-4,4))+
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  xlab('SNP effect size')+
  ylab('-log10(FDR)')+
  labs(color='N. of barcodes')+ 
  #geom_hline(yintercept = 0)+
  geom_hline(yintercept = 0.05)

### differential + active sentinel snps ####
#### plots of sentinel snps with 1 differential Active####
## Activity 
sentinel_snps = read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/sentinel_trait_ld_snps.txt', sep = '\t', header = T)
head(sentinel_snps)
sentinel_snps = sentinel_snps[sentinel_snps$Trait %in% c('DBP','PP','SBP'),]

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

mital_names <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/name_Mital_SNPS.txt', sep='\t', header=T)
#mital_details <- index$parse_details[index$name %in% mital_names$name]
mital_name <- unique(index$dupe_info[index$name %in% mital_names$name])

active_snps_cm <- unique(vals_significance_cm$snp_info[vals_significance_cm$CM_padj<0.05 & !gsub(';.*','',vals_significance_cm$parse_details) %in% mital_names$name])
active_snps_vsmc <- unique(vals_significance_vsmc$snp_info[vals_significance_vsmc$VSMC_padj<0.05 & !gsub(';.*','',vals_significance_vsmc$parse_details) %in% mital_names$name])

CM_results <- CM_results %>%
  separate(name, c("pos", "snp_info","info"), "__")
head(CM_results)
VSMC_results <- VSMC_results %>%
  separate(name, c("pos", "snp_info","info"), "__")
head(VSMC_results)

diff_snps <- unique(CM_results$snp_info[CM_results$fdr_comp<0.05 & CM_results$tile_type == 'WILDTYPE_SNP_INDIV' & !CM_results$dupe_info %in% mital_name])
diff_snps_vsmc <- unique(VSMC_results$snp_info[VSMC_results$fdr_comp<0.05 & VSMC_results$tile_type == 'WILDTYPE_SNP_INDIV' & !VSMC_results$dupe_info %in% mital_name])

sentinel_snps$activeDiff <- 'Not ActiveDiff'
sentinel_snps$activeDiff[sentinel_snps$snp %in% diff_snps_vsmc & sentinel_snps$snp %in% active_snps_vsmc] <- 'ActiveDiff'
sentinel_snps$activeDiffCM <- 'Not ActiveDiff'
sentinel_snps$activeDiffCM[sentinel_snps$snp %in% diff_snps & sentinel_snps$snp %in% active_snps_cm] <- 'ActiveDiff'
sentinel_snps$activeCM <- 'Not Active'
sentinel_snps$activeVSMC <- 'Not Active'
sentinel_snps$activeCM[sentinel_snps$snp %in% active_snps_cm] <- 'Active'
sentinel_snps$activeVSMC[sentinel_snps$snp %in% active_snps_vsmc] <- 'Active'
sentinel_snps$DiffCM <- 'Not Diff'
sentinel_snps$DiffVSMC <- 'Not Diff'
sentinel_snps$DiffCM[sentinel_snps$snp %in% diff_snps] <- 'Diff'
sentinel_snps$DiffVSMC[sentinel_snps$snp %in% diff_snps_vsmc] <- 'Diff'


sentinel <- unique(sentinel_snps$sentinel)

sentinel_snps_sent <- sentinel_snps[sentinel_snps$snp %in% sentinel,]

grouped <- sentinel_snps[,c('sentinel','Chr','Trait','activeDiff')] %>% distinct()

active_sentinel <- unique(grouped$sentinel[grouped$activeDiff == 'Active'])

sentinel_snps_sent$activeDiff[sentinel_snps_sent$sentinel %in% active_sentinel] <- 'ActiveDiff'
dist_snps <- as.data.frame(table(sentinel_snps_sent$activeDiff))

# Compute the position of labels
dist_snps <- dist_snps %>% 
  #arrange(desc(Var1)) %>%
  mutate(prop = Freq / sum(dist_snps$Freq) *100) %>%
  mutate(ypos = (prop)- 0.1*prop )

ggplot(dist_snps, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(width=1, color="white", stat = 'identity') +
  coord_polar("y")+
  geom_text(aes(y = ypos+1, label = paste0(round(Freq),'')), color = "white", size=5) +
  ylab("")+
  xlab("")+
  ggtitle("Sentinel SNPs with at least one LD differential SNP")+
  theme(legend.position="bottom") +
  theme_void() + 
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_fill_brewer(palette="Set2")

number_snps <- as.data.frame(table(sentinel_snps$sentinel[sentinel_snps$DiffVSMC == 'Diff' & sentinel_snps$activeVSMC == 'Active']))
number_snps <- as.data.frame(table(sentinel_snps$sentinel))
ggplot(number_snps, aes(x=Freq)) + 
  geom_density(adjust = 1/3, color='darkgrey')+
  xlab('N of LD SNPs')+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'))

hist(number_snps$Freq, col = "steelblue", frame = FALSE,
     breaks = 60) + abline(v=median(number_snps$Freq), col="red")

### enrichment of repeat families ####
### RepeatMasker ####
repeat_masker <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/snps_repeat_masker.info.txt', sep=' ', header=T)
repeat_masker$query <- gsub('_.*','',repeat_masker$query)
repeat_masker <- repeat_masker %>% distinct()

oligo_info <- vals_significance_cm[,c('dupe_info','pos')]
colnames(oligo_info) <- c('dupe_info','query')

sentinel = merge(sentinel_snps,repeat_masker, by.x= 'coord', by.y= 'query')

ggplot(sentinel[sentinel$DiffVSMC == 'Diff',], aes(x=beginR)) + 
  geom_bar()+
  coord_flip()+
  xlab('N of Diff SNPs')+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50')) +
  theme_minimal()

## enrichment of repeat family

get_col <- function(fdr){
  if(abs(fdr) < -log10(0.05)){
    col <- "light grey"
  }else{
    col <- col_numeric(rev(brewer.pal(11,"RdBu")),
                       domain = seq(-3,
                                    3,length.out = 11))(fdr)
  }
}

active_snps_cm <- unique(vals_significance_cm$pos[vals_significance_cm$CM_padj<0.05 & !vals_significance_cm$parse_details %in% mital_details])
active_snps_vsmc <- unique(vals_significance_vsmc$pos[vals_significance_vsmc$VSMC_padj<0.05 & !vals_significance_cm$parse_details %in% mital_details])

diff_snps <- unique(gsub('__.*','',CM_results$name[CM_results$fdr_comp<0.05 & CM_results$tile_type == 'WILDTYPE_SNP_INDIV' & !CM_results$info %in% mital_details$info]))
diff_snps_vsmc <- unique(gsub('__.*','',VSMC_results$name[VSMC_results$fdr_comp<0.05 & VSMC_results$tile_type == 'WILDTYPE_SNP_INDIV' & !VSMC_results$info %in% mital_details$info]))

library(UpSetR)
list_snps <- list(active_snps_cm[!is.na(active_snps_cm)],
                  active_snps_vsmc[!is.na(active_snps_vsmc)], 
                  diff_snps[!is.na(diff_snps)],
                  diff_snps_vsmc[!is.na(diff_snps_vsmc)],
                  repeat_masker$query)

names(list_snps) <- c('Active CMs', 'Active VSMCs', 'Diff CMs','Diff VSMCs','Repeat')

m2 = list_to_matrix(list_snps)
m2_mat <- make_comb_mat(m2)
UpSet(m2_mat)
UpSet(m2_mat, 
      comb_order = order(comb_size(m2_mat), decreasing = T),
      right_annotation = upset_right_annotation(m2_mat, gp = gpar(fill = c("#474448","#8D818C",'#A5A299','#B4B8C5','#E9EBF8')),
                                                add_numbers = TRUE),
      top_annotation = upset_top_annotation(m2_mat, add_numbers = TRUE))


ht = draw(UpSet(m2_mat[grep('....1',comb_name(m2_mat))], 
                comb_order = order(comb_size(m2_mat[grep('....1',comb_name(m2_mat))]), decreasing = T),
                right_annotation = upset_right_annotation(m2_mat[grep('....1',comb_name(m2_mat))], gp = gpar(fill = c("#474448","#8D818C",'#A5A299','#B4B8C5','#E9EBF8')),
                                                          add_numbers = TRUE),
                top_annotation = upset_top_annotation(m2_mat[grep('....1',comb_name(m2_mat))], add_numbers = TRUE)))
print(ht)
od = column_order(ht)
cs = comb_size(m2_mat[grep('....1',comb_name(m2_mat))])
#Intersection\nsize
decorate_annotation("intersection_size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = "bottom", gp = gpar(fontsize = 8))
})



# enrichment of DSE in particular types of AS events
my_fisher <- function(type){
  #             DS      Not DS
  # type
  # other_types
  print(type)
  repeat_fam1 <- repeat_masker$query[repeat_masker$beginR == type]
  repeat_other <- repeat_masker$query[repeat_masker$beginR != type]
  #type_counts_original <- length(repeat_fam1[repeat_fam1 %in% active_snps_cm & repeat_fam1 %in% diff_snps])
  #other_type_counts_original <- length(repeat_other[repeat_other %in% active_snps_cm & repeat_other %in% diff_snps])
  type_counts_original <- length(repeat_fam1[repeat_fam1 %in% diff_snps_vsmc])
  other_type_counts_original <- length(repeat_other[repeat_other %in% diff_snps_vsmc])
  type_counts_random <- length(repeat_fam1)-type_counts_original
  other_type_counts_random <- length(repeat_other)-other_type_counts_original
  
  m <- matrix(c(type_counts_original, type_counts_random, other_type_counts_original, other_type_counts_random), 2,2, byrow = T)
  print(m)
  m[is.na(m)] <- 0
  rownames(m) <- c("Intron", "Other")
  colnames(m) <- c("Original","Random")
  f <- fisher.test(m)
  print(f)
  return(list("f" = f,
              "overlap" = x11))
}
# Two-tailed Fisher test
families <- unique(repeat_masker$beginR)
fisher_results <- lapply(families, function(chr) my_fisher(chr))
names(fisher_results) <- families
# Plot
par(mfrow=c(2,2), mar=c(5,7,4,4), mgp=c(4,1,0))
odds_ratio <- list()
adj.P.Val <- list()
# odds ratio
odds_ratio <- sapply(families, function(chr) fisher_results[[chr]][['f']]$estimate )
# adj.P.Val
adj.P.Val <- p.adjust(sapply(families, function(chr) fisher_results[[chr]][['f']]$p.value), method = "BH")
print(paste0(sum(adj.P.Val<0.05)))
print(paste0("Depleted: ", sum(odds_ratio[which(adj.P.Val<0.05)]<1))) #  depleted
print(paste0("Enriched: ", sum(odds_ratio[which(adj.P.Val<0.05)]>1))) # enriched
# fdr
fdr <- -log10(adj.P.Val)
fdr[fdr>3] <- 3 # for nice color scale, min adj.P.Val 0.01
fdr[odds_ratio<1] <- -fdr[odds_ratio<1] # if depleted blue

names(odds_ratio) <- families
odds_ratio <- as.data.frame(odds_ratio)

odds_ratio$fam <- rownames(odds_ratio)
colnames(odds_ratio) <- c('odds_ratio','Family')
odds_ratio$pval <- adj.P.Val

threshold_OE <- odds_ratio$pval < 0.05
length(which(threshold_OE))
odds_ratio$thr <- threshold_OE 

ggplot(odds_ratio, aes(x=odds_ratio, y=Family)) + 
  geom_bar(aes(fill=thr),stat = 'identity') +
  ggtitle('location enrichment') +
  ylab('')+
  scale_fill_manual(values = c('grey','Dark red'))+
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'), 
        legend.position = 'none') + 
  geom_vline(aes(xintercept=1), linetype="dotted")

### genomic location of differential SNPs ####
head(sentinel_snps)

number_snps <- as.data.frame(table(sentinel_snps$sentinel[sentinel_snps$DiffCM== 'Diff']))
ggplot(number_snps, aes(x=Freq)) + 
  geom_density(adjust = 1/3, color='darkgrey')+
  xlab('N of LD SNPs')+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'))

hist(number_snps$Freq, col = "steelblue", frame = FALSE,
     breaks = 60)


##### upset plot differential and Active #####
sum(active_snps_vsmc %in% active_snps_cm)

library(ComplexHeatmap)
lt <- list(CM = unique(sentinel_snps$snp[sentinel_snps$activeCM == 'Active']),
           VSMC = unique(sentinel_snps$snp[sentinel_snps$activeVSMC == 'Active']))

m = make_comb_mat(lt)
UpSet(m, comb_order = order(comb_size(m), decreasing = T))

ht = draw(UpSet(m, 
                comb_order = order(comb_size(m), decreasing = T),
                right_annotation = upset_right_annotation(m, gp = gpar(fill = c("#662c91ff","#f15a29ff")),
                                                          add_numbers = TRUE),
                top_annotation = upset_top_annotation(m, add_numbers = TRUE)))
print(ht)
od = column_order(ht)
cs = comb_size(m)
decorate_annotation("intersection_size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = "bottom", gp = gpar(fontsize = 8))
})


###diff
library(ComplexHeatmap)
lt <- list(CM = unique(sentinel_snps$snp[sentinel_snps$DiffCM == 'Diff']),
           VSMC = unique(sentinel_snps$snp[sentinel_snps$DiffVSMC == 'Diff']))

m = make_comb_mat(lt)
UpSet(m, comb_order = order(comb_size(m), decreasing = T))

ht = draw(UpSet(m, 
                comb_order = order(comb_size(m), decreasing = T),
                right_annotation = upset_right_annotation(m, gp = gpar(fill = c("#662c91ff","#f15a29ff")),
                                                          add_numbers = TRUE),
                top_annotation = upset_top_annotation(m, add_numbers = TRUE)))
print(ht)
od = column_order(ht)
cs = comb_size(m)
decorate_annotation("intersection_size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = "bottom", gp = gpar(fontsize = 8))
})

###diffAct
library(ComplexHeatmap)
lt <- list(CM = unique(sentinel_snps$snp[sentinel_snps$DiffCM == 'Diff' & sentinel_snps$activeCM == 'Active']),
           VSMC = unique(sentinel_snps$snp[sentinel_snps$DiffVSMC == 'Diff' & sentinel_snps$activeVSMC == 'Active']))

m = make_comb_mat(lt)
UpSet(m, comb_order = order(comb_size(m), decreasing = T))

ht = draw(UpSet(m, 
                comb_order = order(comb_size(m), decreasing = T),
                right_annotation = upset_right_annotation(m, gp = gpar(fill = c("#662c91ff","#f15a29ff")),
                                                          add_numbers = TRUE),
                top_annotation = upset_top_annotation(m, add_numbers = TRUE)))
print(ht)
od = column_order(ht)
cs = comb_size(m)
decorate_annotation("intersection_size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = "bottom", gp = gpar(fontsize = 8))
})

library(karyoploteR)
sentinel_snps_pos <- sentinel_snps %>%
  separate(change, c("chr", "pos","ref","alt"), ":")
sentinel_snps_pos$chr <- paste0('chr',sentinel_snps_pos$chr)
head(sentinel_snps_pos)
sentine_diff <- sentinel_snps_pos[sentinel_snps_pos$DiffCM == 'Diff',]
sentinel_snps_pos$end <- sentinel_snps_pos$pos
vars <- toGRanges(sentinel_snps_pos[,c('chr','pos','end','snp')])

kp <- plotKaryotype()
kpAddBaseNumbers(kp)
kpPlotDensity(kp, vars)
kpPoints(kp, data=vars, y=0.5)
kpPlotMarkers(kp, chr=sentinel_snps_pos$chr, x=as.numeric(sentinel_snps_pos$pos), labels='')


library(ggbio)
autoplot(vars, layout = "karyogram")
autoplot(dn, layout = "karyogram", aes(color = exReg, fill = exReg), alpha = 0.5) +
  scale_color_discrete(na.value = "brown")

write.table(sentinel_snps_pos[,c('Chr','pos','Trait')], '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/list_snps_trait_simple.txt',
                         sep = '\t',quote = F, row.names = F, col.names = F)
write.table(sentinel_snps_pos[sentinel_snps_pos$DiffCM=='Diff',c('Chr','pos','Trait')], '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/list_snps_trait_simple_diffCM.txt',
            sep = '\t',quote = F, row.names = F, col.names = F)


###### info about indels -----

sentinel_snps <- sentinel_snps %>%
  separate(E, c("REF", "ALT"), ":")
head(sentinel_snps)

sentinel_snps$indel <- 'NO'
sentinel_snps$indel[nchar(sentinel_snps$REF) >1 | nchar(sentinel_snps$ALT) >1] <- 'YES'


##### TF indels differential ----
tf_disruption_analysis <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/02_TF_analysis/turnover_results_tf_diff.txt')
index <- index %>%
  separate(name, c("pos", "snp_info","info"), "__")
head(index)

### merge data
tf_disruption_analysis <- merge(tf_disruption_analysis, index[,c("snp_info",'element')], by.x = 'element_alt', by.y = 'element')
tf_disruption_analysis <- tf_disruption_analysis %>% distinct()

### merge with snp info 
tf_disruption_analysis_snp <- merge(tf_disruption_analysis, sentinel_snps, by.x = 'snp_info', by.y = 'snp')
tf_disruption_analysis_snp <- tf_disruption_analysis_snp %>% distinct()

#### plots
library(ggpubr)
# Change outline and fill colors by groups ("sex")
# Use a custom palette
ggdensity(tf_disruption_analysis_snp, x = "total_motifs",
          add = "mean", rug = TRUE,
          color = "indel", fill = "indel",
          palette = c("#0073C2FF", "#FC4E07"))

ggdensity(tf_disruption_analysis_snp, x = "perc_shared_motifs",
          add = "mean", rug = TRUE,
          color = "indel", fill = "indel",
          palette = c("#0073C2FF", "#FC4E07"))

wilcox.test(perc_shared_motifs ~ indel, data=tf_disruption_analysis_snp) 

ggdensity(tf_disruption_analysis_snp, x = "abs_delta_motifs",
          add = "mean", rug = TRUE,
          color = "indel", fill = "indel",
          palette = c("#0073C2FF", "#FC4E07"))

a <- ggplot(tf_disruption_analysis_snp, aes(x = perc_shared_motifs))
a + stat_ecdf(aes(color = indel,linetype = indel), 
              geom = "step", size = 1.5) +
  scale_color_manual(values = c("#00AFBB", "#E7B800"))+
  labs(y = "f(% shared motif)")

### differential 
ggdensity(tf_disruption_analysis_snp[tf_disruption_analysis_snp$DiffVSMC == 'Diff',], x = "perc_shared_motifs",
          add = "mean", rug = TRUE,
          color = "indel", fill = "indel",
          palette = c("#0073C2FF", "#FC4E07"))

ggdensity(tf_disruption_analysis_snp[tf_disruption_analysis_snp$DiffCM == 'Diff',], x = "abs_delta_motifs",
          add = "mean", rug = TRUE,
          color = "indel", fill = "indel",
          palette = c("#0073C2FF", "#FC4E07"))

wilcox.test(perc_shared_motifs ~ indel, data=tf_disruption_analysis_snp[tf_disruption_analysis_snp$DiffVSMC== 'Diff',]) 


