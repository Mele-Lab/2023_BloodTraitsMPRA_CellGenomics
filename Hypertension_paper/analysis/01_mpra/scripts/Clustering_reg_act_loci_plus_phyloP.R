##### 
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
#mital_name <- unique(index$dupe_info[index$name %in% mital_names$name])

active_snps_cm <- unique(vals_significance_cm$snp_info[vals_significance_cm$CM_padj<0.05 & !vals_significance_cm$full_name %in% mital_names$name])
active_snps_vsmc <- unique(vals_significance_vsmc$snp_info[vals_significance_vsmc$VSMC_padj<0.05])

CM_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_ordered.CM.005.new_back.txt', sep='\t', header = T)
VSMC_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_ordered.VSMC.005.new_back.txt', sep='\t', header = T)

diff_snps <- unique(CM_results$snp[CM_results$fdr_comp<0.05 & CM_results$tile_type == 'WILDTYPE_SNP_INDIV'])
diff_snps_vsmc <- unique(VSMC_results$snp[VSMC_results$fdr_comp<0.05 & VSMC_results$tile_type == 'WILDTYPE_SNP_INDIV'])

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

#### read length block -----
length <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/sentinel_snps_length_block.txt',
                     sep = '\t', header = T)
#1,000,000

number_snps <- as.data.frame(table(sentinel_snps$sentinel[sentinel_snps$activeDiffCM == 'ActiveDiff']))
number_snps_all <- as.data.frame(table(sentinel_snps$sentinel))
number_snps_all <- merge(number_snps_all, length, by.x='Var1', by.y='sentinel')
number_snps_all$clust <- number_snps_all$Freq/(number_snps_all$length/1000000)

number_snps <- merge(number_snps, length, by.x='Var1', by.y='sentinel')
number_snps <- merge(number_snps, number_snps_all, by='Var1')
#number_snps$clust <- number_snps$Freq/(number_snps$length/1000000)

ggplot(number_snps[number_snps$length>1,], aes(x=clust)) + 
  geom_density(adjust = 1/3, color='darkgrey')+
  xlab('N of LD SNPs')+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'))

number_snps$clust <- as.numeric(number_snps$clust)
number_snps$length <- as.numeric(number_snps$length)
vector <- number_snps_all$clust[number_snps_all$length>1]
hist(log10(vector), col = 'steelblue', breaks = 60, main = 'Nº of LD SNPs / Mb', xlab = 'log10(NºSNPs/Mbase)') + 
  abline(v=median(log10(vector)), col='red')

ggplot(number_snps[number_snps$length.y>300,], aes(x=clust, y = Freq.x)) + 
  geom_point(color='darkgrey')+
  xlab('Nº SNPs / Mb')+
  ylab('Nº of Active + Differential SNPs')+
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=FALSE) +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50')) +
  theme_minimal()+
  geom_point() + 
  geom_text_repel(
    data = number_snps[number_snps$length.y>300,],
    aes(label = Var1),
    size = 4)
  geom_text(
    label=number_snps$Var1[number_snps$length.y>300], 
    nudge_x = 0.25, nudge_y = 10, 
    check_overlap = T
  )
  
#### plot ratio ------
number_snps <- as.data.frame(table(sentinel_snps$sentinel[sentinel_snps$activeDiffCM == 'ActiveDiff']))
number_snps <- merge(number_snps, number_snps_all, by='Var1')
  
number_snps$ratio <- (number_snps$Freq.x/number_snps$clust)
vector <- number_snps$Freq.x/number_snps$clust
hist((vector), col = 'steelblue', breaks = 60, main = 'Nº of Act SNPs / Mb', xlab = 'Active SNPs / clustering') + 
  abline(v=median((vector)), col='red')

write.table(number_snps[order(number_snps$ratio, decreasing = T),], 
            '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/ratio_regAct_CM_ordered.txt',
            sep = '\t', col.names = T, row.names = F, quote = F)

#### phylop score -----
  phylop_scores <- read.csv('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/04_PhyloP/PhyloP_scores_summary_tile.txt',
                            sep = '\t', header = T)
  
  phylop_scores <- read.csv('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/04_PhyloP/PhyloP_scores_summary_snp.txt',
                            sep = '\t', header = T)
  phylop_scores <- read.csv('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/04_PhyloP/PhyloP_scores_summary_tile_46_placenta_window_a.txt',
                            sep = '\t', header = T)
  phylop_scores_B <- read.csv('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/04_PhyloP/PhyloP_scores_summary_tile_46_placenta_window_b.txt',
                            sep = '\t', header = T)
  
  phylop_window <- rbind(phylop_scores, phylop_scores_B)
  
phylop_scores_g <- phylop_window %>% 
  group_by(snp, sentinel) %>% 
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  as.data.frame()  

head(phylop_scores_g)
write.table(phylop_scores_g, '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/04_PhyloP/PhyloP_scores_summary_tile_46_placenta_window_a_b.txt',
            sep = '\t', col.names = T, row.names = F, quote = F)

#### take one random oligo -----
phylop_scores_g_r <- lapply(c(1:500), function(x) {
  phylop_scores_g %>% 
  group_by(sentinel) %>% 
  sample_n(1) %>%
  as.data.frame()  
  })

phylop_scores_g_r_df <- do.call('rbind.data.frame',phylop_scores_g_r)
phylop_scores_g_r <- phylop_scores_g_r_df %>% 
  group_by(sentinel) %>% 
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  as.data.frame()  
  
number_snps_phylop <- merge(phylop_scores_g_r, number_snps_all, by.x='sentinel',by.y='Var1')

ggplot(number_snps_phylop[number_snps_phylop$length>300,], aes(x=clust, y = mean)) + 
  geom_point(color='darkgrey')+
  xlab('Nº SNPs / Mb')+
  ylab('mean PhyloP score')+
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=FALSE) +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50')) +
  theme_minimal()+
  geom_point() + 
  geom_text_repel(
    data = number_snps_phylop[number_snps_phylop$length>300,],
    aes(label = sentinel),
    size = 4)
geom_text(
  label=number_snps_phylop$sentinel[number_snps_phylop$length>300], 
  nudge_x = 0.25, nudge_y = 10, 
  check_overlap = T
)    

library("ggpubr")
ggscatter(number_snps_phylop[number_snps_phylop$length>300,], x = "clust", y = "mean", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Clustering", ylab = "mean(PhyloP)")

#### phyloP per type -----
phylop_type <- merge(phylop_scores, sentinel_snps, by.x='snp',by.y='snp')

p <- ggboxplot(phylop_type, x = "DiffCM", y = "mean",
               color = "DiffCM", palette = "jco")
p + stat_compare_means(aes(group = DiffCM) , label = "p.format")

library(ggridges)
library(viridis)
library(hrbrthemes)

ggplot(phylop_type, aes(x = mean, y = sig_both, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  labs(title = '') +
  theme_minimal() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )

library(plyr)
phylop_type$sig_both <- 'no sig'
phylop_type$sig_both[phylop_type$DiffCM == 'Diff'] <- 'sig CM'
phylop_type$sig_both[phylop_type$DiffVSMC == 'Diff'] <- 'sig VSMC'
phylop_type$sig_both[phylop_type$DiffCM == 'Diff' & phylop_type$DiffVSMC == 'Diff'] <- 'sig both'


mu <- ddply(phylop_type, "sig_both", summarise, grp.mean=mean(mean, na.rm = TRUE))
head(mu)

p<-ggplot(phylop_type, aes(x=mean, color=sig_both)) +
  geom_density()+
  xlim(c(-2,2))+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=sig_both),
             linetype="dashed") + theme_minimal()
p


p <- ggboxplot(phylop_type, x = "sig_both", y = "mean",
               color = "sig_both", palette = "jco")
p +   stat_compare_means(method = "anova", label.y = 6.5)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "no sig")  

### clustering snps ------
phylop_scores_freq_g <- phylop_scores_freq %>% 
  group_by(sentinel,diffVSMC) %>% 
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  as.data.frame()
number_snps_phylop <- merge(phylop_scores_g, number_snps, by.x='sentinel',by.y='Var1')
number_snps_phylop$clust_d <- number_snps$Freq.x / number_snps$Freq.y

ggplot(number_snps_phylop[number_snps_phylop$length.x>300,], aes(x=clust_d, y = mean)) + 
  geom_point(color='darkgrey')+
  xlab('Nº SNPs / Mb')+
  ylab('mean PhyloP score')+
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=FALSE) +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50')) +
  theme_minimal()+
  geom_point() + 
  geom_text_repel(
    data = number_snps_phylop[number_snps_phylop$length.x>300,],
    aes(label = sentinel),
    size = 4)

library("ggpubr")
ggscatter(number_snps_phylop[number_snps_phylop$length>300,], x = "clust", y = "mean", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Clustering", ylab = "mean(PhyloP)")

#### merge phylop and snps/mb --------
phylop_scores <- read.csv('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/04_PhyloP/PhyloP_scores_summary_tile.txt',
                          sep = '\t', header = T)

#phylop_scores <- read.csv('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/04_PhyloP/PhyloP_scores_summary_snp.txt',
#                          sep = '\t', header = T)

phylop_scores_g <- phylop_scores %>% 
  group_by(sentinel) %>% 
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  as.data.frame()  

number_snps <- as.data.frame(table(sentinel_snps$sentinel[sentinel_snps$DiffVSMC == 'Diff']))

number_snps <- merge(number_snps, length, by.x='Var1', by.y='sentinel')
number_snps <- merge(number_snps, number_snps_all, by='Var1')
#number_snps$clust <- number_snps$Freq/(number_snps$length/1000000)

number_snps_phylop <- merge(phylop_scores_g, number_snps, by.x='sentinel',by.y='Var1')

ggplot(number_snps_phylop[number_snps_phylop$Freq.x < 20 &number_snps_phylop$length.x >300 ,], aes(x=Freq.x, y = mean, color=clust)) + 
  geom_point()+
  scale_color_gradient2(low="black", mid='black', high="red", midpoint = 100)+
  xlab('Nº SNPs / Mb')+
  ylab('mean PhyloP score')+
  geom_smooth(method=lm , color="red", se=FALSE) +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50')) +
  theme_minimal()+
  geom_text_repel(
    data = number_snps_phylop[number_snps_phylop$Freq.x < 20 & number_snps_phylop$length.x >300,],
    aes(label = sentinel),
    size = 4)

ggscatter(number_snps_phylop[number_snps_phylop$Freq.x < 20 & number_snps_phylop$length.x >300,], x = "Freq.x", y = "mean", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Clustering", ylab = "mean(PhyloP)")

#### number of genes per LD block -----
library(GenomicRanges)
library(AnnotationDbi)
library(AnnotationHub)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

info <- as.data.frame(table(sentinel_snps$sentinel))
info <- merge(info, length, by.x='Var1', by.y='sentinel')
info <- merge(info, sentinel_snps[,c('sentinel','Chr')], by.x='Var1', by.y='sentinel')
info <- info %>% distinct()
info$clust <- info$Freq/(info$length/1000000)

colnames(info) <- c('sentinel','Freq','start','end','width','chr','clustering')

granges_trait <- makeGRangesFromDataFrame(
  info[,c('chr','start','end','width','sentinel','Freq','clustering')],
  seqnames.field = 'chr',
  start.field = 'start',
  end.field = 'end',
  keep.extra.columns = TRUE)
seqlevelsStyle(granges_trait) <- "UCSC"

transcriptsByGene = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
overlaps = findOverlaps(granges_trait, transcriptsByGene)

queryHits(overlaps)
subjectHits(overlaps)

number_genes <- as.data.frame(table(queryHits(overlaps)))
info$Var1 <- rownames(info)
number_genes <- merge(number_genes, info, by='Var1')

ggplot(number_genes[number_genes$width >300,], aes(x=clustering, y = Freq.x)) + 
  geom_point()+
  xlab('Nº SNPs / Mb')+
  ylab('Nº genes')+
  geom_smooth(method=lm , color="red", se=FALSE) +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50')) +
  theme_minimal()+
  geom_text_repel(
    data = number_genes[number_genes$width >300,],
    aes(label = sentinel),
    size = 4)

write.table(number_genes, '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/number_genes_per_sentinel.txt',
            sep = '\t', quote = F, row.names = F, col.names = T)


#### genes name per sentinel ------
library(org.Hs.eg.db)

closest_gene <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/clostest_gene_snps.txt', sep='\t', header=F)
head(sentinel_snps)
closest_gene$V4 <- gsub('\\|.*','', closest_gene$V4)

merge_info <- merge(sentinel_snps, closest_gene, by.x='change', by.y='V4', all.x=TRUE)
head(merge_info)

entrez_DE <- mapIds(org.Hs.eg.db,keys=gsub('\\..*','',merge_info$V7),column="SYMBOL", keytype="ENSEMBL", multiVals="first")
merge_info$hgnc <- entrez_DE

grouped_snps <- merge_info[!is.na(merge_info$hgnc),] %>% 
  group_by(sentinel) %>% 
  summarise(genes = paste(unique(hgnc), collapse = ',')) %>%
  as.data.frame()  

ratio_ordered <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/ratio_active_CM_ordered.txt',
                            sep = '\t', header = T)
ratio_ordered <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/ratio_active_VSMC_ordered.txt',
                            sep = '\t', header = T)
ratio_ordered <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/ratio_reg_CM_ordered.txt',
                            sep = '\t', header = T)
ratio_ordered <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/ratio_reg_VSMC_ordered.txt',
                            sep = '\t', header = T)
ratio_ordered <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/ratio_regAct_CM_ordered.txt',
                            sep = '\t', header = T)
ratio_ordered <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/ratio_regAct_VSMC_ordered.txt',
                            sep = '\t', header = T)

ratio <- merge(ratio_ordered, grouped_snps, all.x=T, by.x='Var1', by.y='sentinel')
colnames(ratio) <- c('sentinel', 'N.sig','SNPs in LD','start block','end block','length block','SNPs/MB','ratio','nearest genes')
write.table(ratio[order(ratio$ratio, decreasing = T),], '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/ratio_Act_CM_ordered.genes.txt',
            sep = '\t', quote = F, col.names = T, row.names = F)

