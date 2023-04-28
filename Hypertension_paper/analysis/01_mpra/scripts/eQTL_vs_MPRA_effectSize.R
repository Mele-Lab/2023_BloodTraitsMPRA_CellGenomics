#### Concordance MPRA - eQTL analysis ####
###gtex snps ----
### GTEx eQTLs ####
#### eQTL data downloaded from GTEx portal ######
gtex_path <- '~/marenostrum/Projects/GTEx_v8/Manuscript/'
TissueInfoFile <- "TissuesInfo.rds" 

TissueInfo <- readRDS(paste0(gtex_path,TissueInfoFile))
tissues <- as.character(TissueInfo$Tissue_id)
subset <- tissues[grep('Arter|Heart',tissues)]

eQTLs_path <- '~/marenostrum_scratch/bsc83535/GTEx/v8/cis_QTLs/cis_eQTLs/GTEx_Analysis_v8_eQTL/'

eQTLs <- lapply(subset, function(tissue) read.table(paste0(eQTLs_path,tissue,".v8.signif_variant_gene_pairs.txt.gz"), sep='\t', header=T) )
names(eQTLs) <- subset

eQTLs <- lapply(subset, function(tissue) eQTLs[[tissue]] %>%
                  separate(variant_id, c("chrom", "start","ref","alt",'genome'), "_"))
names(eQTLs) <- subset
eQTLs_n <- lapply(subset, function(tissue) eQTLs[[tissue]]$variant_id <- paste0(eQTLs[[tissue]]$chrom,
                                                                              '_',eQTLs[[tissue]]$start,
                                                                              '_',eQTLs[[tissue]]$ref,
                                                                              '_',eQTLs[[tissue]]$alt,
                                                                              '_',eQTLs[[tissue]]$genome))
names(eQTLs_n) <- subset

liftover_coord <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/snp_coord_master.hg38.simple.chr.sorted.bed', sep='\t', header=F)
colnames(liftover_coord) <- c('chr_b38', 'start_b38', 'end_b38','info_snp')

sentinel_snps = read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/sentinel_trait_ld_snps.txt', sep = '\t', header = T)
head(sentinel_snps)
sentinel_snps = sentinel_snps[sentinel_snps$Trait %in% c('DBP','PP','SBP'),]
liftover_coord <- merge(liftover_coord, sentinel_snps, by.x='info_snp', by.y='change')
liftover_coord <- liftover_coord %>% distinct()

liftover_coord$variant_id <- paste0(liftover_coord$chr_b38,'_',as.numeric(liftover_coord$start_b38)+1,'_',gsub(':.*','',liftover_coord$E),'_',gsub('.*:','',liftover_coord$E),'_b38')
head(liftover_coord)

## slope > 0 --> ALT more active
## slope < 0 --> REF more active 
liftover_coord$start_b38 <-  liftover_coord$start_b38+1
eQTLs_mpra <- lapply(subset, function(tissue) as.data.frame(merge(liftover_coord, eQTLs[[tissue]][,c("chrom", "start","ref","alt","slope","gene_id")], 
                                                                  by.x=c('chr_b38', 'start_b38'), 
                                                                  by.y=c("chrom", "start"))))
names(eQTLs_mpra) <- subset
eQTLs_mpra_S <- lapply(subset, function(tissue) as.data.frame(eQTLs_mpra[[tissue]][eQTLs_mpra[[tissue]]$variant_id %in% eQTLs_n[[tissue]],]))
names(eQTLs_mpra_S) <- subset

#### merge with mpra results #####
CM_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_with_significance.CM.new_back.005.txt', sep='\t', header = T)
VSMC_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_with_significance.VSMC.new_back.005.txt', sep='\t', header = T)

index <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/Hypertension__pooled.index.txt', sep = '\t', header = T)
head(index)

index['dupe_info'] <- gsub('\\..*','',index$tile_id)

CM_results = merge(CM_results,index[index$tile_type != 'RANDOM',c("name","dupe_info")], by="dupe_info")
VSMC_results = merge(VSMC_results,index[index$tile_type != 'RANDOM',c("name","dupe_info")], by="dupe_info")

CM_results <- CM_results %>% distinct()
VSMC_results <- VSMC_results %>% distinct()

CM_results <- CM_results %>%
  separate(name, c("info","snp","coord"), "__")
VSMC_results <- VSMC_results %>%
  separate(name, c("info","snp","coord"), "__")

## FC > 0 --> REF more active
## FC < 0 --> ALT more active

eQTLs_mpra_CM <- lapply(subset, function(tissue) as.data.frame(merge(eQTLs_mpra_S[[tissue]], CM_results[,c('snp','logFC_comp','coord','fdr_comp')], 
                                                                      by.x=c('snp'), 
                                                                      by.y=c('snp'))))
names(eQTLs_mpra_CM) <- subset
eQTLs_mpra_VSMC <- lapply(subset, function(tissue) as.data.frame(merge(eQTLs_mpra_S[[tissue]], VSMC_results[,c('snp','logFC_comp','coord','fdr_comp')], 
                                                                     by.x=c('snp'), 
                                                                     by.y=c('snp'))))
names(eQTLs_mpra_VSMC) <- subset

library(dplyr) #basic data management & %>%
library(likert) 

CM_clust <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/CM_clustering_info.txt')
VSMC_clust <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/VSMC_clustering_info.txt')

for (tissue in subset) {
  data <- eQTLs_mpra_CM[[tissue]]
  data$concordance <- 'discordant'
  data$concordance[sign(data$slope) != sign(data$logFC_comp)] <- 'concordant'
  data$regulatory <- 'not regulatory'
  data$regulatory[data$fdr_comp < 0.05] <- 'regulatory'
  
  ggplot(data, aes(x=slope)) + geom_density() +
    ggtitle(paste0(tissue,' betas eQTLs all'))
  
  pdf(file = paste0("~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/plots/",tissue,"_cumulative_mpra_reg_conc.pdf"),   # The directory you want to save the file in
      width = 6, # The width of the plot in inches
      height = 4)
  print(ggplot(data, aes(abs(logFC_comp), color=interaction(regulatory, concordance), 
                   group=interaction(regulatory, concordance))) + 
    stat_ecdf(geom = "step")+
    labs(title=paste0(tissue, " MPRA [Effect Size]"),
         y = "Empirical Cumulative density", x="[MPRA Effect Size]")+
    theme_classic())
  dev.off()
  
  data$EffectSize <- 'Negative'
  data$EffectSize[sign(data$logFC_comp)<0] <- 'Positive'
  data$strong <- 'Low'
  data$strong[abs(data$logFC_comp)>0.8 & data$fdr_comp < 0.05] <- 'strong'
  table(data$strong)
  
  pdf(file = paste0("~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/plots/",tissue,"_eqTL_beta_strong_reg.pdf"),   # The directory you want to save the file in
      width = 4, # The width of the plot in inches
      height = 3)
  print(ggplot(data[data$regulatory == 'regulatory' & data$strong == 'strong',], aes(x=slope, fill = EffectSize)) + 
    geom_density() +
    theme_classic() +
    ggtitle(paste0(tissue,' betas eQTLs strong regulatory')))
  dev.off()
  
  mat <- table(data$concordance, data$strong)
  mat <- mat[,c(2,1)]
  print(mat)
  mat[is.na(mat)] <- 0
  f <- fisher.test(mat)
  print(f)

  sentinel_concord <- as.data.frame(table(data$sentinel[data$regulatory=='regulatory'], 
        data$concordance[data$regulatory=='regulatory']))
  
  pdf(file = paste0("~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/plots/",tissue,"_sentinel_concord.pdf"),   # The directory you want to save the file in
      width = 4, # The width of the plot in inches
      height = 6)
  #plot(likert(sentinel_concord))

  print(ggplot(sentinel_concord[sentinel_concord$Var1 %in% CM_clust$sentinel[CM_clust$clustering_group=='Top density'],], aes(x=Freq, y=Var1, fill = Var2)) + 
    geom_bar(position="dodge", stat="identity") +
    theme_classic() +
    ylab('')+ xlab('nº reg variants')+
    ggtitle(paste0(tissue,' betas concordant reg variants')))
  dev.off()
  
  n_pos <- nrow(data[data$concordance == 'concordant' & data$sentinel %in% CM_clust$sentinel[CM_clust$clustering_group=='Top density'] & data$regulatory == 'regulatory',])
  total <- nrow(data[data$sentinel %in% CM_clust$sentinel[CM_clust$clustering_group=='Top density'] & data$regulatory == 'regulatory',])
  
  binomial <- binom.test(n_pos, total, 0.5)$p.value
  print(binomial)
  
}

for (tissue in subset) {
  data <- eQTLs_mpra_CM[[tissue]]
  data$concordance <- 'discordant'
  data$concordance[sign(data$slope) != sign(data$logFC_comp)] <- 'concordant'
  data$regulatory <- 'not regulatory'
  data$regulatory[data$fdr_comp < 0.05] <- 'regulatory'
  
  sentinel_concord <- as.data.frame(table(data$gene_id[data$regulatory=='regulatory'], 
                                          data$concordance[data$regulatory=='regulatory']))
  
  symbol_DE <- mapIds(org.Hs.eg.db,keys=gsub('\\..*','',sentinel_concord$Var1),column="SYMBOL", keytype="ENSEMBL", multiVals="first")
  sentinel_concord$gene <- symbol_DE
  
  sentinel_concord$gene[is.na(sentinel_concord$gene)] <- as.character(sentinel_concord$Var1[is.na(sentinel_concord$gene)]) 
  
  pdf(file = paste0("~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/plots/",tissue,"_topgenes_concord.pdf"),   # The directory you want to save the file in
      width = 4, # The width of the plot in inches
      height = 6)
  #plot(likert(sentinel_concord))
  
  print(ggplot(sentinel_concord[sentinel_concord$Freq > 8,], aes(x=Freq, y=gene, fill = Var2)) + 
          geom_bar(position="dodge", stat="identity") +
          theme_classic() +
          ylab('')+ xlab('nº reg variants')+
          ggtitle(paste0(tissue,' betas concordant reg variants')))
  dev.off()

}
