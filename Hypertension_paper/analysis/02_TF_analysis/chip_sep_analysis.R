##### Chip-seq analysis ------
#### read data ----
sentinel_snps = read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/sentinel_trait_ld_snps.txt', sep = '\t', header = T)
head(sentinel_snps)
sentinel_snps = sentinel_snps[sentinel_snps$Trait %in% c('DBP','PP','SBP'),]

###### chip-atlas ------
chip_atlas <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/overlap_chip_atlas.new.bed',
                        sep = '\t', header = F)
chip_atlas <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/overlap_chip_atlas.new.Kidney.bed',
                         sep = '\t', header = F)
chip_atlas <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/overlap_chip_atlas.new.Neuronal.bed',
                         sep = '\t', header = F)
head(chip_atlas)


#### chip-atlas 
sentinel_atlas <- merge(sentinel_snps, chip_atlas, by.x = 'change',by.y = 'V4', all.x = T)
head(sentinel_atlas)

head(sentinel_atlas)
sentinel_atlas$TF <- gsub('.*=','', gsub('%.*', '', sentinel_atlas$V8))

sentinel_atlas_table_A <- as.data.frame(table(sentinel_atlas$TF))
ggplot(sentinel_atlas_table_A[sentinel_atlas_table_A$Freq > 25,], aes(x=Var1, y=Freq)) +
  geom_bar(stat = 'identity') +
  theme_minimal() +
  coord_flip()


### which snps have more chip-seq? -----
CM_results_s <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/CM_eqtls_rep_cons_dhs_tbs_scores.txt', sep = '\t',
                           header = T)

VSMC_results_s <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/VSMC_eqtls_rep_cons_dhs_tbs_scores.txt', sep = '\t',
                             header = T)

sentinel_atlas_table_snp <- as.data.frame(table(sentinel_atlas$snp[!is.na(sentinel_atlas$TF)]))
sentinel_remap_n_chip <- merge(CM_results_s, sentinel_atlas_table_snp, by.x = 'snp_info',by.y = 'Var1', all.x = T)

hist(sentinel_remap_n_chip$Freq, col = "steelblue", frame = FALSE) +
  abline(v=median(sentinel_remap_n_chip$Freq, na.rm = T), col="red")

library(ggpubr)
ggscatter(sentinel_remap_n_chip, x = "Freq", y = "logFC_comp",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Nº Chip", ylab = "abs(logFC")

##### enrichment of mapping at reg variants #####
sentinel_atlas$reg_CM <- 'not reg'
sentinel_atlas$reg_CM[sentinel_atlas$snp %in% CM_results_s$snp_info] <- 'reg'
table(sentinel_atlas$reg_CM)

sentinel_atlas$reg_VSMC <- 'not reg'
sentinel_atlas$reg_VSMC[sentinel_atlas$snp %in% VSMC_results_s$snp_info] <- 'reg'
table(sentinel_atlas$reg_VSMC)

length(unique(sentinel_atlas$snp[!is.na(sentinel_atlas$TF) & sentinel_atlas$reg_CM == 'reg']))
length(unique(sentinel_atlas$snp[!is.na(sentinel_atlas$TF) & sentinel_atlas$reg_VSMC == 'reg']))

### fisher enrichment ####
bound_reg <- length(unique(sentinel_atlas$snp[!is.na(sentinel_atlas$TF) & sentinel_atlas$reg_CM == 'reg']))
bound_notreg <- length(unique(sentinel_atlas$snp[!is.na(sentinel_atlas$TF) & sentinel_atlas$reg_CM != 'reg']))
not_bound_reg <- length(unique(sentinel_atlas$snp[sentinel_atlas$reg_CM == 'reg'])) - bound_reg
not_bound_not_reg <- length(unique(sentinel_atlas$snp[sentinel_atlas$reg_CM != 'reg'])) - bound_notreg

bound_reg <- length(unique(sentinel_atlas$snp[!is.na(sentinel_atlas$TF) & sentinel_atlas$reg_VSMC == 'reg']))
bound_notreg <- length(unique(sentinel_atlas$snp[!is.na(sentinel_atlas$TF) & sentinel_atlas$reg_VSMC != 'reg']))
not_bound_reg <- length(unique(sentinel_atlas$snp[sentinel_atlas$reg_VSMC == 'reg'])) - bound_reg
not_bound_not_reg <- length(unique(sentinel_atlas$snp[sentinel_atlas$reg_VSMC != 'reg'])) - bound_notreg

m <- matrix(c(bound_reg, bound_notreg, not_bound_reg, not_bound_not_reg), 2,2, byrow = T)
print(m)
m[is.na(m)] <- 0
rownames(m) <- c("TF", "Other")
colnames(m) <- c("reg","not reg")
f <- fisher.test(m, alternative = 'greater')
print(f)

phyper(bound_reg-1, length(unique(sentinel_atlas$snp[sentinel_atlas$reg_CM == 'reg'])), length(unique(sentinel_atlas$snp[sentinel_atlas$reg_CM != 'reg'])), length(unique(sentinel_atlas$snp[!is.na(sentinel_atlas$TF)])), lower.tail = FALSE)

### are reg variants having more peaks?
frequencies <- as.data.frame(table(sentinel_atlas$snp))

#### activity and differential results -----
CM_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_ordered.CM.005.new_back.txt', sep='\t', header = T)
VSMC_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_ordered.VSMC.005.new_back.txt', sep='\t', header = T)

diff_snps <- unique(CM_results$snp[CM_results$fdr_comp<0.05])
diff_snps_vsmc <- unique(VSMC_results$snp[VSMC_results$fdr_comp<0.05])

frequencies$reg_CM <- 'not reg'
frequencies$reg_CM[frequencies$Var1 %in% diff_snps] <- 'reg'
table(frequencies$reg_CM)

frequencies$reg_VSMC <- 'not reg'
frequencies$reg_VSMC[frequencies$Var1 %in% diff_snps_vsmc] <- 'reg'
table(frequencies$reg_VSMC)

library(ggpubr)
p <- ggboxplot(frequencies, x = "reg_VSMC", y = "Freq",
               color = "reg_VSMC", palette = "jco",
               add = "jitter")
#  Add p-value
p + stat_compare_means()
# Change method
p + stat_compare_means(method = "t.test")


#### sig_TFBS_from_FIMO #####
sig_motifs_cm <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/TF_Analysis/sig_motifs.txt', sep = '\t', header = T)
tf <- unique(sig_motifs_cm$HGNC.symbol)
(tf %in% unique(sentinel_atlas$TF[!is.na(sentinel_atlas$TF)]))

head(sentinel_atlas)
sentinel_atlas <- sentinel_atlas[,c("change","sentinel","name","TF","reg_CM","reg_VSMC")] %>% distinct()

sentinel_atlas_table <- as.data.frame(table(sentinel_atlas$TF[!is.na(sentinel_atlas$TF) & sentinel_atlas$reg_CM == 'reg']))
sentinel_atlas_table_vsmc <- as.data.frame(table(sentinel_atlas$TF[!is.na(sentinel_atlas$TF) & sentinel_atlas$reg_VSMC == 'reg']))
sentinel_atlas_table_all <- as.data.frame(table(sentinel_atlas$TF[!is.na(sentinel_atlas$TF)]))

head(sentinel_atlas_table_vsmc)

write.table(sentinel_atlas_table$Var1, 
            '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/02_TF_analysis/Tf_bound_reg_CM.v2.txt',
            sep = '\n', quote = F, row.names = F, col.names = F)
write.table(sentinel_atlas_table_vsmc$Var1, 
            '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/02_TF_analysis/Tf_bound_reg_VSMC.v2.txt',
            sep = '\n', quote = F, row.names = F, col.names = F)
write.table(sentinel_atlas_table_A$Var1, 
            '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/02_TF_analysis/Tf_bound_reg_all.v2.txt',
            sep = '\n', quote = F, row.names = F, col.names = F)


library(VennDiagram)
venn.diagram(
  x = list(
    sentinel_atlas_table$Var1 , 
    sentinel_atlas_table_vsmc$Var1 
  ),
  category.names = c("TF in reg CMs" , "TF in reg VSMC"),
  filename = 'venn.Neurons.svg',
  output = TRUE ,
  imagetype="svg" ,
  height = 7 , 
  width = 7 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3))
)

phyper(91-1, 114, 66, 93, lower.tail = FALSE, log.p = FALSE)

colnames(sentinel_atlas_table) <- c('TF', 'CM')
colnames(sentinel_atlas_table_vsmc) <- c('TF','VSMC')

TF_chip <- merge(sentinel_atlas_table, sentinel_atlas_table_vsmc)
library(reshape2)
TF_chip <- TF_chip[TF_chip$CM > 10,]
mdata <- melt(TF_chip, id=c("TF"))
cvs_samples <- CVS_df[CVS_df$Var1 %in% mdata$TF,]
kid_samples <- Kidney_df[Kidney_df$Var1 %in% mdata$TF,]
neu_samples <- neural_df[neural_df$Var1 %in% mdata$TF,]
mdata <- mdata[mdata$TF %in% neu_samples$Var1,]

plot1 <- ggplot(mdata, aes(y=value, x=TF, fill=variable)) + 
  geom_bar(stat = 'identity',position="dodge") + 
  theme_classic() + 
  coord_flip() + xlab('') + ylab('Nº of regulatory variants w/ TF bound') +
  theme(axis.text.y=element_blank(), #remove x axis labels
        axis.ticks.y=element_blank()) #remove x axis ticks
plot2 <- ggplot(neu_samples, aes(y=Freq, x=Var1)) + 
  geom_bar(stat = 'identity',color='gray') + 
  theme_classic() + 
  coord_flip() + xlab('TF') + ylab('Nº of cell-types')
  

library(egg)
ggarrange(plot2, plot1, widths = c(1.5,2))
library(cowplot)
plot_grid(plot2, plot1, rel_widths = c(0.4,0.6))

my_fisher <- function(type){
  #             DS      Not DS
  # type
  # other_types
  print(type)
  chip_fam1 <- sentinel_atlas[sentinel_atlas$TF %in% c(type),]
  chip_other <- sentinel_atlas[!sentinel_atlas$TF %in% c(type) & !is.na(sentinel_atlas$TF),]
  #type_counts_original <- length(repeat_fam1[repeat_fam1 %in% active_snps_cm & repeat_fam1 %in% diff_snps])
  #other_type_counts_original <- length(repeat_other[repeat_other %in% active_snps_cm & repeat_other %in% diff_snps])
  type_counts_original <- nrow((chip_fam1[chip_fam1$reg_VSMC %in% c('reg'),]))
  other_type_counts_original <- nrow((chip_other[chip_other$reg_VSMC %in% c('reg'),]))
  type_counts_random <- nrow(chip_fam1)-type_counts_original
  other_type_counts_random <- nrow(chip_other)-other_type_counts_original
  
  m <- matrix(c(type_counts_original, type_counts_random, other_type_counts_original, other_type_counts_random), 2,2, byrow = T)
  print(m)
  m[is.na(m)] <- 0
  rownames(m) <- c("TF", "Other")
  colnames(m) <- c("reg","not reg")
  f <- fisher.test(m)
  print(f)
  return(list("f" = f,
              "overlap" = x11))
}
# Two-tailed Fisher test
families <- unique(sentinel_atlas$TF)[!is.na(unique(sentinel_atlas$TF))]
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
odds_ratio$pval <- sapply(families, function(chr) fisher_results[[chr]][['f']]$p.value)
odds_ratio$adjpval <- adj.P.Val

threshold_OE <- odds_ratio$pval < 0.05
length(which(threshold_OE))
odds_ratio$thr <- threshold_OE 

ggplot(odds_ratio[odds_ratio$pval<0.4,], aes(x=odds_ratio, y=Family)) + 
  geom_bar(aes(fill=thr),stat = 'identity') +
  ggtitle('location enrichment') +
  ylab('')+
  scale_fill_manual(values = c('grey','Dark red'))+
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'), 
        legend.position = 'none') + 
  geom_vline(aes(xintercept=1), linetype="dotted")


### more bound in reg? 

###### more differential in repeats? -------

#### read data -----
library(GenomicRanges)
sentinel_snps = read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/sentinel_trait_ld_snps.txt', sep = '\t', header = T)
head(sentinel_snps)
sentinel_snps = sentinel_snps[sentinel_snps$Trait %in% c('DBP','PP','SBP'),]
sentinel_snps <- sentinel_snps %>% distinct()

#### activity and differential results -----
CM_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_ordered.CM.005.new_back.txt', sep='\t', header = T)
VSMC_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_ordered.VSMC.005.new_back.txt', sep='\t', header = T)

diff_snps <- unique(CM_results$snp[CM_results$fdr_comp<0.05])
diff_snps_vsmc <- unique(VSMC_results$snp[VSMC_results$fdr_comp<0.05])

sentinel_snps$DiffCM <- 'Not Diff'
sentinel_snps$DiffVSMC <- 'Not Diff'
sentinel_snps$DiffCM[sentinel_snps$snp %in% diff_snps] <- 'Diff'
sentinel_snps$DiffVSMC[sentinel_snps$snp %in% diff_snps_vsmc] <- 'Diff'

chip_atlas <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/overlap_chip_atlas.new.bed',
                         sep = '\t', header = F)
head(chip_atlas)

### we need to get info in regions -----
head(sentinel_snps)
sentinel_snps_cord <- as.list(sentinel_snps$coord)
sentinel_snps_gr = lapply(sentinel_snps_cord, function (x) {res=strsplit(x, ':')}) %>%
  unlist %>%
  matrix(ncol=3, byrow=T) %>%
  as.data.frame %>%
  dplyr::select(chrom=V1, start=V2, end=V3) %>%
  makeGRangesFromDataFrame


chip_gr = chip_atlas %>%
  dplyr::select(chrom=V5, start=V6, end=V7) %>% distinct() %>% 
  makeGRangesFromDataFrame

diff_cord <- as.list(sentinel_snps$coord[sentinel_snps$DiffCM == 'Diff'])
diff_gr = lapply(diff_cord, function (x) {res=strsplit(x, ':')}) %>%
  unlist %>%
  matrix(ncol=3, byrow=T) %>%
  as.data.frame %>%
  dplyr::select(chrom=V1, start=V2, end=V3) %>%
  makeGRangesFromDataFrame

##### make test -------
library(regioneR)
pt <- permTest(A=diff_gr, B=chip_gr, randomize.function=resampleRegions, universe=sentinel_snps_gr,
               evaluate.function=numOverlaps, ntimes=1000,verbose=FALSE)
summary(pt)
pt
plot(pt)


#### compare kidney, neuronal, cvd
###### chip-atlas ------
chip_atlas <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/overlap_chip_atlas.new.bed',
                         sep = '\t', header = F)
chip_atlas_kidney <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/overlap_chip_atlas.new.Kidney.bed',
                         sep = '\t', header = F)
chip_atlas_neuronal <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/overlap_chip_atlas.new.Neuronal.bed',
                         sep = '\t', header = F)
head(chip_atlas)

sentinel_atlas <- merge(sentinel_snps, chip_atlas, by.x = 'change',by.y = 'V4', all.x = T)
sentinel_atlas_kidney <- merge(sentinel_snps, chip_atlas_kidney, by.x = 'change',by.y = 'V4', all.x = T)
sentinel_atlas_neuronal <- merge(sentinel_snps, chip_atlas_neuronal, by.x = 'change',by.y = 'V4', all.x = T)

sentinel_atlas$TF <- gsub('.*=','', gsub('%.*', '', sentinel_atlas$V8))
sentinel_atlas_kidney$TF <- gsub('.*=','', gsub('%.*', '', sentinel_atlas_kidney$V8))
sentinel_atlas_neuronal$TF <- gsub('.*=','', gsub('%.*', '', sentinel_atlas_neuronal$V8))

sentinel_atlas$reg_CM <- 'not reg'
sentinel_atlas$reg_CM[sentinel_atlas$snp %in% CM_results_s$snp_info] <- 'reg'
table(sentinel_atlas$reg_CM)

sentinel_atlas$reg_VSMC <- 'not reg'
sentinel_atlas$reg_VSMC[sentinel_atlas$snp %in% VSMC_results_s$snp_info] <- 'reg'
table(sentinel_atlas$reg_VSMC)

length(unique(sentinel_atlas$snp[!is.na(sentinel_atlas$TF) & sentinel_atlas$reg_CM == 'reg']))
length(unique(sentinel_atlas$snp[!is.na(sentinel_atlas$TF) & sentinel_atlas$reg_VSMC == 'reg']))

library(VennDiagram)
venn.diagram(
  x = list(
    sentinel_atlas_table$Var1 , 
    sentinel_atlas_table_vsmc$Var1 
  ),
  category.names = c("TF in reg CMs" , "TF in reg VSMC"),
  filename = 'venn.Neurons.svg',
  output = TRUE ,
  imagetype="svg" ,
  height = 7 , 
  width = 7 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3))
)

venn.diagram(
  x = list(
    unique(sentinel_atlas_neuronal$snp[!is.na(sentinel_atlas_neuronal$TF) & sentinel_atlas_neuronal$snp %in% unique(sentinel_atlas$snp[sentinel_atlas$reg_VSMC=='reg'])]), 
    unique(sentinel_atlas_kidney$snp[!is.na(sentinel_atlas_kidney$TF) & sentinel_atlas_kidney$snp %in% unique(sentinel_atlas$snp[sentinel_atlas$reg_VSMC=='reg'])]),
    unique(sentinel_atlas$snp[!is.na(sentinel_atlas$TF) & sentinel_atlas$reg_VSMC=='reg'])
  ),
  category.names = c("Neural" , "Kidney","CVS"),
  filename = 'venn.TFVSMCs.svg',
  output = TRUE ,
  imagetype="svg" ,
  height = 7 , 
  width = 7 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  cat.cex = 2,
  cex=2,
  col=c("#440154ff", '#21908dff','#fde725ff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3),alpha('#fde725ff',0.3)))

## overlap TFs expresed 
chip_seq_info <- read.delim('~/Downloads/experimentList.tab', header = F, sep = '\t' 
                            # col.names = c('ID','Genome','Antigen','Cell_type_Class',
                            #               'cell_type','Desc','logs','logsBS','Title','name',
                            #               'cell-line','chip_ab','number')
                            )
head(chip_seq_info)
chip_seq_info <- chip_seq_info[chip_seq_info$V2=='hg38',]
chip_seq_info <- chip_seq_info[chip_seq_info$V3=='TFs and others',]
head(chip_seq_info)
neural <- chip_seq_info[chip_seq_info$V5=='Neural',]
Kidney <- chip_seq_info[chip_seq_info$V5=='Kidney',]
CVS <- chip_seq_info[chip_seq_info$V5=='Cardiovascular',]

neural_df <- as.data.frame(table(neural$V4))
Kidney_df <- as.data.frame(table(Kidney$V4))
CVS_df <- as.data.frame(table(CVS$V4))

venn.diagram(
  x = list(
    neural_df$Var1, 
    Kidney_df$Var1,
    CVS_df$Var1
  ),
  category.names = c("Neural" , "Kidney","CVS"),
  filename = 'venn.TFall.svg',
  output = TRUE ,
  imagetype="svg" ,
  height = 7 , 
  width = 7 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  cat.cex = 2,
  cex=2,
  col=c("#440154ff", '#21908dff','#fde725ff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3),alpha('#fde725ff',0.3))
)

  